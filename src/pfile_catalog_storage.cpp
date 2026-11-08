#include "pfile_catalog_storage.hpp"
#include "pfile_catalog_matrix.hpp"
#include "pfile_catalog.hpp"
#include "pfile_catalog_util.hpp"

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
#include "duckdb/common/enums/access_mode.hpp"

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

//! v1 block tiling. v1's manifest carries one concat axis of `partition_start/stop` + `concat_dim`;
//! the other axis is shared (one block, extent not declared). `axis` is "variant" or "sample".
static string BuildV1BlockSql(const string &path_sql, const string &concat_dim, const string &axis) {
	bool axis_is_concat = (concat_dim == axis + "s") || (concat_dim.empty() && axis == "variant");
	if (axis_is_concat) {
		return "SELECT DISTINCT partition_start AS " + axis + "_start, partition_stop AS " + axis +
		       "_stop FROM read_pfile_catalog(" + path_sql + ") ORDER BY 1";
	}
	// Shared axis: one block; the v1 manifest does not declare the shared-axis extent.
	return "SELECT NULL::BIGINT AS " + axis + "_start, NULL::BIGINT AS " + axis + "_stop";
}

//! c.genotypes / c.genotypes_by_sample.
//! A multi-file cohort that shares one .psam (variant-sharded) must pass it explicitly — read_pfile
//! applies it to the first source as the shared sample metadata. A single fileset passes NOTHING:
//! read_pfile discovers its own companions (incl. .pvar.parquet / .bim / .fam). A pvar override is
//! never passed — read_pfile rejects it on a multi-file list, and each cell/shard already carries
//! its own co-located .pvar.
static string GenotypesSql(const vector<string> &prefixes, const string &shared_psam, bool by_sample) {
	string sql = "SELECT * FROM read_pfile(" + SqlPathList(prefixes);
	if (prefixes.size() > 1 && !shared_psam.empty()) {
		sql += ", psam := " + SQLString::ToString(shared_psam);
	}
	if (by_sample) {
		sql += ", orient := 'sample'";
	}
	sql += ")";
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
	string raw_path = info.path;

	// Resolve the ATTACH target: an explicit '.pfile.json', a prefix ('<p>.pfile.json'), or a bare
	// fileset ('<p>.pgen' with no sidecar → synthesize a one-cell catalog, P13-005 D1).
	string catalog_path; // the .pfile.json to parse; empty for a bare fileset
	string bare_prefix;  // the fileset prefix; empty unless bare
	if (StringUtil::EndsWith(raw_path, ".pfile.json")) {
		if (!fs.FileExists(raw_path)) {
			throw InvalidInputException("ATTACH pfile_catalog: catalog file not found: '%s'", raw_path);
		}
		catalog_path = raw_path;
	} else if (fs.FileExists(raw_path + ".pfile.json")) {
		catalog_path = raw_path + ".pfile.json";
	} else if (fs.FileExists(raw_path + ".pgen")) {
		bare_prefix = raw_path;
	} else {
		throw InvalidInputException("ATTACH pfile_catalog: nothing found at '%s' (looked for it, "
		                            "'%s.pfile.json', and '%s.pgen')",
		                            raw_path, raw_path, raw_path);
	}

	// Resolve the cell list (prefix/pvar/psam per surviving cell) and the per-kind structural-view
	// SQL. Absolute paths so the lazy views resolve the same regardless of the cwd at query time.
	vector<string> cell_prefixes;
	vector<string> pvar_paths;
	vector<string> psam_paths;
	string variant_blocks_sql;
	string sample_blocks_sql;
	string cells_sql;

	auto take_cell = [&](const string &prefix, const string &pvar, const string &psam) {
		// pvar_paths/psam_paths feed c.variants (read_pvar) and c.samples (read_psam), which need an
		// explicit path. A catalog that omits a companion falls back to the conventional
		// <prefix>.pvar / <prefix>.psam. (A .bim/.fam/.parquet-only fileset is not covered for those
		// two relations; c.genotypes still discovers it, since GenotypesSql passes no per-cell path.)
		cell_prefixes.push_back(prefix);
		PushUnique(pvar_paths, pvar.empty() ? prefix + ".pvar" : pvar);
		PushUnique(psam_paths, psam.empty() ? prefix + ".psam" : psam);
	};

	if (!bare_prefix.empty()) {
		// Bare fileset: one synthetic cell, no catalog file. read_pfile discovers its own companions.
		string abs_prefix = AbsPath(fs, bare_prefix);
		take_cell(abs_prefix, "", "");
		variant_blocks_sql = "SELECT NULL::BIGINT AS variant_start, NULL::BIGINT AS variant_stop";
		sample_blocks_sql = "SELECT NULL::BIGINT AS sample_start, NULL::BIGINT AS sample_stop";
		cells_sql = "SELECT " + SQLString::ToString(abs_prefix) + " AS prefix, " +
		            SQLString::ToString(abs_prefix + ".pvar") + " AS pvar, " +
		            SQLString::ToString(abs_prefix + ".psam") + " AS psam";
	} else {
		string content = ReadWholeFile(fs, catalog_path);
		string path_sql = SQLString::ToString(AbsPath(fs, catalog_path));
		cells_sql = "SELECT * FROM read_pfile_catalog(" + path_sql + ")";
		if (IsMatrixCatalog(content)) {
			auto matrix = ParseMatrixCatalog(fs, catalog_path, content);
			for (auto &cell : BuildCells(matrix, fs, /*include_holes=*/false)) {
				if (cell.present) {
					take_cell(cell.prefix, cell.pvar, cell.psam);
				}
			}
			variant_blocks_sql = BuildBlockViewSql(matrix.variants, path_sql, "variant_start", "variant_stop");
			sample_blocks_sql = BuildBlockViewSql(matrix.samples, path_sql, "sample_start", "sample_stop");
		} else {
			auto cat = ParseV1Catalog(fs, catalog_path, content);
			if (cat.concat_dim == "samples") {
				// Sample-sharded genotypes need read_pfile's sample-concat (CONCATENATE), which is not
				// implemented, and a pvar override on the shared .pvar is rejected on a multi-file list.
				// Fail loud with a kind-specific message rather than a cryptic reader error on query.
				throw InvalidInputException(
				    "ATTACH pfile_catalog '%s': sample-sharded cohorts are not supported yet — the "
				    "genotype relations need read_pfile's sample-concat, which is not implemented. "
				    "Variant-sharded, single-fileset, and matrix catalogs attach.",
				    raw_path);
			}
			for (auto &e : cat.files) {
				take_cell(e.prefix, e.pvar, e.psam);
			}
			variant_blocks_sql = BuildV1BlockSql(path_sql, cat.concat_dim, "variant");
			sample_blocks_sql = BuildV1BlockSql(path_sql, cat.concat_dim, "sample");
		}
	}

	// A catalog with no present cells has nothing to read (the shipped readers reject an empty list).
	if (cell_prefixes.empty()) {
		throw InvalidInputException("ATTACH pfile_catalog '%s': the catalog has no present cells; "
		                            "nothing to attach",
		                            raw_path);
	}

	// Bake the per-relation SQL over the shipped readers (one builder per relation, so P13-006 can
	// swap c.genotypes's backing to the pushdown cascade in one place). A shared .psam (one distinct
	// psam across >1 cell = a variant-sharded cohort) is passed to read_pfile explicitly; a single
	// fileset lets read_pfile discover its own companions.
	string shared_psam = (cell_prefixes.size() > 1 && psam_paths.size() == 1) ? psam_paths[0] : "";
	vector<std::pair<Identifier, string>> views;
	views.emplace_back(Identifier("genotypes"), GenotypesSql(cell_prefixes, shared_psam, false));
	views.emplace_back(Identifier("genotypes_by_sample"), GenotypesSql(cell_prefixes, shared_psam, true));
	views.emplace_back(Identifier("variants"), "SELECT * FROM read_pvar(" + SqlPathList(pvar_paths) + ")");
	views.emplace_back(Identifier("samples"), BuildPsamUnionSql(psam_paths));
	views.emplace_back(Identifier("cells"), cells_sql);
	views.emplace_back(Identifier("variant_blocks"), variant_blocks_sql);
	views.emplace_back(Identifier("sample_blocks"), sample_blocks_sql);

	// Build a native in-memory catalog seeded with those views. The backing storage is :memory:,
	// which cannot be opened read-only — so force READ_WRITE here even if the user passed READ_ONLY;
	// the catalog's read-only guarantee comes from SetReadOnlyDatabase() above, not the storage mode.
	info.path = ":memory:";
	options.access_mode = AccessMode::READ_WRITE;
	auto catalog = make_uniq<DuckCatalog>(db);
	catalog->Initialize(false);

	auto system_transaction = CatalogTransaction::GetSystemTransaction(db.GetDatabase());
	auto &schema = catalog->GetSchema(system_transaction, Identifier::DefaultSchema());
	auto &duck_schema = schema.Cast<DuckSchemaEntry>();
	auto &view_set = duck_schema.GetCatalogSet(CatalogType::VIEW_ENTRY);
	view_set.SetDefaultGenerator(make_uniq<PfileCatalogViewGenerator>(*catalog, schema, std::move(views)));

	return std::move(catalog);
}

static unique_ptr<TransactionManager> PfileCatalogTransactionManager(optional_ptr<StorageExtensionInfo> storage_info,
                                                                     AttachedDatabase &db, Catalog &catalog) {
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
