#include "pfile_catalog_storage.hpp"
#include "pfile_catalog_matrix.hpp"

#include "duckdb/parser/parser.hpp"
#include "duckdb/storage/storage_extension.hpp"
#include "duckdb/parser/parsed_data/attach_info.hpp"
#include "duckdb/parser/parsed_data/create_view_info.hpp"
#include "duckdb/catalog/duck_catalog.hpp"
#include "duckdb/transaction/duck_transaction_manager.hpp"
#include "duckdb/catalog/default/default_views.hpp"
#include "duckdb/catalog/catalog_entry/duck_schema_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/main/attached_database.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/exception.hpp"

namespace duckdb {

//! A default-view generator seeded at attach time with (view name -> SELECT sql).
//! Modelled on DuckDB's OpenFileDefaultGenerator: views are created lazily, on first
//! reference, so ATTACH itself opens no data files.
class PfileCatalogViewGenerator : public DefaultGenerator {
public:
	PfileCatalogViewGenerator(Catalog &catalog, SchemaCatalogEntry &schema,
	                          vector<std::pair<Identifier, string>> views_p)
	    : DefaultGenerator(catalog), schema(schema), views(std::move(views_p)) {
		for (auto &v : views) {
			view_names.push_back(v.first);
		}
	}

	unique_ptr<CatalogEntry> CreateDefaultEntry(ClientContext &context, const Identifier &entry_name) override {
		for (auto &v : views) {
			if (entry_name == v.first) {
				auto result = make_uniq<CreateViewInfo>();
				result->SetQualifiedName(QualifiedName({Identifier::DefaultSchema()}, v.first));
				result->sql = v.second;
				auto parser = Parser::GetBuiltinParser();
				auto view_info = CreateViewInfo::FromSelect(parser, std::move(result));
				return make_uniq_base<CatalogEntry, ViewCatalogEntry>(catalog, schema, *view_info);
			}
		}
		return nullptr;
	}

	vector<Identifier> GetDefaultEntries() override {
		return view_names;
	}

private:
	SchemaCatalogEntry &schema;
	vector<std::pair<Identifier, string>> views;
	vector<Identifier> view_names;
};

//! A SQL list literal of single-quoted, escaped paths: ['a', 'b', ...].
static string SqlPathList(const vector<string> &paths) {
	string out = "[";
	for (idx_t i = 0; i < paths.size(); i++) {
		if (i > 0) {
			out += ", ";
		}
		out += SQLString(paths[i]);
	}
	out += "]";
	return out;
}

//! A double-quoted SQL identifier (internal quotes doubled).
static string QuoteIdent(const string &name) {
	return "\"" + StringUtil::Replace(name, "\"", "\"\"") + "\"";
}

//! Per-axis block tiling: DISTINCT (axis key values, span) over the layout reader
//! (which includes holes, so the tiling is complete even if a block's cells are all holes).
static string BuildBlockViewSql(const AxisDef &axis, const string &path_sql, const char *start_col,
                                const char *stop_col) {
	string sql = "SELECT DISTINCT ";
	for (auto &key : axis.keys) {
		sql += "keys[" + SQLString(key.name) + "] AS " + QuoteIdent(key.name) + ", ";
	}
	sql += string(start_col) + ", " + string(stop_col) + " FROM read_pfile_catalog_layout(" + path_sql + ")";
	return sql;
}

//! Append `s` to `v` only if non-empty and not already present (small N, order preserved).
static void PushUnique(vector<string> &v, const string &s) {
	if (s.empty()) {
		return;
	}
	for (auto &e : v) {
		if (e == s) {
			return;
		}
	}
	v.push_back(s);
}

//! c.samples: read_psam is single-file, so union the distinct sample-block .psam BY NAME.
static string BuildPsamUnionSql(const vector<string> &psams) {
	string sql;
	for (idx_t i = 0; i < psams.size(); i++) {
		if (i > 0) {
			sql += " UNION ALL BY NAME ";
		}
		sql += "SELECT * FROM read_psam(" + SQLString::ToString(psams[i]) + ")";
	}
	return sql;
}

static string ReadWholeFile(FileSystem &fs, const string &path) {
	string content;
	auto handle = fs.OpenFile(path, FileFlags::FILE_FLAGS_READ);
	idx_t size = handle->GetFileSize();
	content.resize(size);
	if (size > 0) {
		handle->Read(reinterpret_cast<void *>(&content[0]), size);
	}
	return content;
}

static unique_ptr<Catalog> PfileCatalogAttach(optional_ptr<StorageExtensionInfo> storage_info, ClientContext &context,
                                              AttachedDatabase &db, const string &name, AttachInfo &info,
                                              AttachOptions &options) {
	// The catalog is a read-only view over the shipped readers. access_mode is already
	// fixed before attach() runs, so force the read-only type here (rejects CREATE/INSERT/etc).
	db.SetReadOnlyDatabase();

	auto &fs = FileSystem::GetFileSystem(context);
	string catalog_path = info.path;
	string content = ReadWholeFile(fs, catalog_path);

	// Parse + resolve the surviving cell grid (opens no data files).
	auto matrix = ParseMatrixCatalog(fs, catalog_path, content);
	auto cells = BuildCells(matrix, fs, /*include_holes=*/false);

	vector<string> cell_prefixes;
	vector<string> pvar_paths;
	vector<string> psam_paths;
	for (auto &cell : cells) {
		if (cell.present) {
			cell_prefixes.push_back(cell.prefix);
			PushUnique(pvar_paths, cell.pvar);
			PushUnique(psam_paths, cell.psam);
		}
	}

	// An all-holes catalog has no data to read. The shipped read_pfile/read_pvar reject an empty
	// file list, so rather than fabricate reader schemas we fail loud here. (Consistency note:
	// read_pfile_catalog returns empty for all-holes; giving read_pfile an empty-list path is a
	// separate reader change — flagged in the ledger.)
	if (cell_prefixes.empty()) {
		throw InvalidInputException("ATTACH pfile_catalog '%s': the catalog has no present cells "
		                            "(all cells are holes); nothing to attach",
		                            catalog_path);
	}

	// Bake the per-relation SQL over the shipped readers.
	string path_sql = SQLString::ToString(catalog_path);
	vector<std::pair<Identifier, string>> views;
	views.emplace_back(Identifier("genotypes"), "SELECT * FROM read_pfile(" + SqlPathList(cell_prefixes) + ")");
	views.emplace_back(Identifier("cells"), "SELECT * FROM read_pfile_catalog(" + path_sql + ")");
	views.emplace_back(Identifier("variant_blocks"),
	                   BuildBlockViewSql(matrix.variants, path_sql, "variant_start", "variant_stop"));
	views.emplace_back(Identifier("sample_blocks"),
	                   BuildBlockViewSql(matrix.samples, path_sql, "sample_start", "sample_stop"));
	views.emplace_back(Identifier("variants"), "SELECT * FROM read_pvar(" + SqlPathList(pvar_paths) + ")");
	views.emplace_back(Identifier("samples"), BuildPsamUnionSql(psam_paths));
	views.emplace_back(Identifier("genotypes_by_sample"),
	                   "SELECT * FROM read_pfile(" + SqlPathList(cell_prefixes) + ", orient := 'sample')");

	// Build a native in-memory catalog seeded with those views.
	info.path = ":memory:";
	auto catalog = make_uniq<DuckCatalog>(db);
	catalog->Initialize(false);

	auto system_transaction = CatalogTransaction::GetSystemTransaction(db.GetDatabase());
	auto &schema = catalog->GetSchema(system_transaction, Identifier::DefaultSchema());
	auto &duck_schema = schema.Cast<DuckSchemaEntry>();
	auto &view_set = duck_schema.GetCatalogSet(CatalogType::VIEW_ENTRY);
	view_set.SetDefaultGenerator(make_uniq<PfileCatalogViewGenerator>(*catalog, schema, std::move(views)));

	return std::move(catalog);
}

static unique_ptr<TransactionManager>
PfileCatalogTransactionManager(optional_ptr<StorageExtensionInfo> storage_info, AttachedDatabase &db, Catalog &catalog) {
	return make_uniq<DuckTransactionManager>(db);
}

void RegisterPfileCatalogStorage(ExtensionLoader &loader) {
	auto extension = make_shared_ptr<StorageExtension>();
	extension->attach = PfileCatalogAttach;
	extension->create_transaction_manager = PfileCatalogTransactionManager;

	auto &config = DBConfig::GetConfig(loader.GetDatabaseInstance());
	StorageExtension::Register(config, "pfile_catalog", std::move(extension));
}

} // namespace duckdb
