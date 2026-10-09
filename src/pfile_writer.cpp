#include "pfile_writer.hpp"

#include "duckdb/common/file_system.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/function/copy_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <set>

namespace duckdb {

// ===========================================================================
// FORMAT pvar — text `.pvar` writer (P9 Phase A, first increment)
//
// Writes a `.pvar` text file from a query, re-readable by `read_pvar` to the
// same decoded values. Variant-metadata only; genotype/sample files are later
// increments. See docs/planning/P9-001-writers-design (Decisions 1, 2, 6, 9).
//
// Surface (DuckDB COPY options are `KEY value`, space-separated — NOT `:=`,
// which the COPY parser rejects; the P9 design docs show `:=` but that is the
// table-function argument operator, not COPY option syntax):
//   COPY (SELECT ... ORDER BY CHROM, POS) TO 'out.pvar' (FORMAT pvar
//     [, variant_fields {chrom:'CHROM', pos:'POS', id:'ID', ref:'REF', alt:'ALT'}]
//     [, pvar 'override/path.pvar']);
// variant_fields is a STRUCT literal: bareword keys (roles), string-literal
// values (input column names) — both constant-fold at bind.
//
// Defaults (reader-aligned, so `COPY (SELECT * FROM read_pvar('x')) TO 'y'`
// round-trips with no options): each required role maps to its same-named
// input column (CHROM/POS/ID/REF/ALT, matched case-insensitively), and every
// other input column is carried through as an extra `.pvar` column under its
// own name, in input order — mirroring `read_pvar`'s output schema.
// ===========================================================================

namespace {

// The five required variant roles, in canonical `.pvar` output order. The
// output header always uses these exact (uppercase) names for the roles, even
// when `variant_fields` remaps a role to a differently-named input column.
constexpr const char *kRoleKeys[] = {"chrom", "pos", "id", "ref", "alt"};
constexpr const char *kRoleHeaders[] = {"CHROM", "POS", "ID", "REF", "ALT"};
constexpr idx_t kRoleCount = 5;

bool IsIntegerType(const LogicalType &type) {
	switch (type.id()) {
	case LogicalTypeId::TINYINT:
	case LogicalTypeId::SMALLINT:
	case LogicalTypeId::INTEGER:
	case LogicalTypeId::BIGINT:
	case LogicalTypeId::HUGEINT:
	case LogicalTypeId::UTINYINT:
	case LogicalTypeId::USMALLINT:
	case LogicalTypeId::UINTEGER:
	case LogicalTypeId::UBIGINT:
	case LogicalTypeId::UHUGEINT:
		return true;
	default:
		return false;
	}
}

// Case-insensitive lookup of `target` in `names`; returns the index or throws
// an actionable error naming the available columns.
idx_t ResolveColumn(const vector<string> &names, const string &target, const string &role) {
	for (idx_t i = 0; i < names.size(); i++) {
		if (StringUtil::CIEquals(names[i], target)) {
			return i;
		}
	}
	throw BinderException("write pvar: required field '%s' maps to column '%s', which is not in the input "
	                      "(available columns: %s). Supply variant_fields {%s: '<col>'} to remap it.",
	                      role, target, StringUtil::Join(names, ", "), role);
}

struct PvarWriteBindData : public FunctionData {
	//! Output header names, in output order: the five canonical role headers
	//! followed by any carried-through extra columns.
	vector<string> out_headers;
	//! For each output column, the index into the input chunk to read from.
	vector<idx_t> out_source_idx;
	//! Input column index for the pos / chrom roles (for the monotonicity check).
	idx_t pos_input_idx = 0;
	idx_t chrom_input_idx = 0;
	//! File-path override from the `pvar` option; empty means use DuckDB's COPY target.
	string path_override;

	unique_ptr<FunctionData> Copy() const override {
		auto r = make_uniq<PvarWriteBindData>();
		r->out_headers = out_headers;
		r->out_source_idx = out_source_idx;
		r->pos_input_idx = pos_input_idx;
		r->chrom_input_idx = chrom_input_idx;
		r->path_override = path_override;
		return std::move(r);
	}
	bool Equals(const FunctionData &other_p) const override {
		auto &o = other_p.Cast<PvarWriteBindData>();
		return out_headers == o.out_headers && out_source_idx == o.out_source_idx &&
		       pos_input_idx == o.pos_input_idx && chrom_input_idx == o.chrom_input_idx &&
		       path_override == o.path_override;
	}
};

struct PvarWriteGlobalState : public GlobalFunctionData {
	unique_ptr<FileHandle> handle;
	string path; //! the path actually opened (for cleanup-on-failure)

	// (CHROM, POS) monotonicity tracking. The sink is single-threaded + ordered
	// (REGULAR_COPY_TO_FILE), so this global state sees rows in query order.
	bool has_last = false;
	string last_chrom;
	int64_t last_pos = 0;
	std::set<string> seen_chroms;

	//! Close the handle and remove the partial file. Used both on the failure
	//! path (sink throw) and so cleanup is idempotent. DuckDB only auto-removes
	//! its single COPY target; a `pvar` override would otherwise be orphaned.
	void CleanupPartial(FileSystem &fs) {
		if (handle) {
			handle.reset();
		}
		try {
			if (!path.empty() && fs.FileExists(path)) {
				fs.RemoveFile(path);
			}
		} catch (...) {
			// Best-effort cleanup; never mask the original error.
		}
	}
};

struct PvarWriteLocalState : public LocalFunctionData {};

// Convert one cell to its `.pvar` text form. NULL -> "." (round-trips with the
// reader, which maps "." -> NULL for every column type).
//
// A `.pvar` is strictly tab-delimited with newline-terminated records, and the
// reader drops every '\r'. A cell containing any of those bytes cannot be
// represented without an escaping scheme (which PLINK's format does not define),
// and writing it raw would silently shift/split fields on read-back. Reject it
// with an actionable error rather than corrupt the file.
string CellToField(const Value &v, const string &col_name) {
	if (v.IsNull()) {
		return ".";
	}
	auto s = v.ToString();
	if (s.find('\t') != string::npos || s.find('\n') != string::npos || s.find('\r') != string::npos) {
		throw InvalidInputException("write pvar: field '%s' value contains a tab, newline, or carriage return, "
		                            "which cannot be written to a tab-delimited .pvar",
		                            col_name);
	}
	return s;
}

unique_ptr<FunctionData> PvarCopyBind(ClientContext &, CopyFunctionBindInput &input, const vector<string> &names,
                                      const vector<LogicalType> &sql_types) {
	auto bind = make_uniq<PvarWriteBindData>();

	// --- parse options ------------------------------------------------------
	// role -> input column name; defaults map each role to its same-named column.
	vector<string> role_target(kRoleCount);
	for (idx_t r = 0; r < kRoleCount; r++) {
		role_target[r] = kRoleHeaders[r];
	}

	for (auto &opt : input.info.options) {
		auto key = StringUtil::Lower(opt.first);
		if (key == "variant_fields") {
			if (opt.second.empty()) {
				continue;
			}
			auto &val = opt.second[0];
			if (val.type().id() != LogicalTypeId::STRUCT) {
				throw BinderException("write pvar: variant_fields must be a STRUCT of role -> column name, e.g. "
				                      "{chrom: 'CHROM', pos: 'POS'}");
			}
			auto &child_types = StructType::GetChildTypes(val.type());
			auto &children = StructValue::GetChildren(val);
			for (idx_t c = 0; c < child_types.size(); c++) {
				auto role = StringUtil::Lower(child_types[c].first);
				idx_t r = kRoleCount;
				for (idx_t i = 0; i < kRoleCount; i++) {
					if (role == kRoleKeys[i]) {
						r = i;
						break;
					}
				}
				if (r == kRoleCount) {
					throw BinderException("write pvar: unknown variant_fields role '%s' (expected one of "
					                      "chrom, pos, id, ref, alt)",
					                      child_types[c].first);
				}
				if (children[c].IsNull()) {
					throw BinderException("write pvar: variant_fields role '%s' has a NULL column name", role);
				}
				role_target[r] = children[c].ToString();
			}
		} else if (key == "pvar") {
			if (!opt.second.empty() && !opt.second[0].IsNull()) {
				bind->path_override = opt.second[0].ToString();
			}
		}
		// Other options (sample_fields / genotype_fields / samples / path
		// overrides for other members) belong to later increments; ignore here.
	}

	// --- resolve roles against the input schema -----------------------------
	vector<idx_t> role_idx(kRoleCount);
	for (idx_t r = 0; r < kRoleCount; r++) {
		role_idx[r] = ResolveColumn(names, role_target[r], kRoleKeys[r]);
	}
	bind->chrom_input_idx = role_idx[0];
	bind->pos_input_idx = role_idx[1];

	// POS must be an integer type so the monotonicity check compares positions
	// numerically (the reader emits POS as INTEGER).
	if (!IsIntegerType(sql_types[role_idx[1]])) {
		throw BinderException("write pvar: pos field maps to column '%s' of type %s, but POS must be an integer "
		                      "type",
		                      names[role_idx[1]], sql_types[role_idx[1]].ToString());
	}

	// --- build the output column set -----------------------------------------
	// Canonical role headers first (sourced from their mapped input columns)...
	std::set<idx_t> consumed;
	for (idx_t r = 0; r < kRoleCount; r++) {
		bind->out_headers.emplace_back(kRoleHeaders[r]);
		bind->out_source_idx.push_back(role_idx[r]);
		consumed.insert(role_idx[r]);
	}
	// ...then every remaining input column, carried through under its own name,
	// in input order (mirrors read_pvar's "extra columns" output).
	for (idx_t i = 0; i < names.size(); i++) {
		if (consumed.count(i) == 0) {
			bind->out_headers.push_back(names[i]);
			bind->out_source_idx.push_back(i);
		}
	}

	// A role remapped off a column whose name equals a canonical header, while
	// that same column is also carried through as an extra, would emit a
	// duplicate output header (e.g. "#CHROM ... CHROM"), which read_pvar then
	// refuses. Catch it at bind with an actionable error (case-insensitive,
	// since DuckDB column binding is case-insensitive).
	for (idx_t a = 0; a < bind->out_headers.size(); a++) {
		for (idx_t b = a + 1; b < bind->out_headers.size(); b++) {
			if (StringUtil::CIEquals(bind->out_headers[a], bind->out_headers[b])) {
				throw BinderException("write pvar: output column '%s' appears twice in the .pvar header — a "
				                      "variant_fields role was remapped off a column that is also carried "
				                      "through as an extra. Drop or rename the duplicate input column.",
				                      bind->out_headers[b]);
			}
		}
	}

	return std::move(bind);
}

unique_ptr<GlobalFunctionData> PvarCopyInitGlobal(ClientContext &context, FunctionData &bind_data,
                                                  const string &file_path) {
	auto &bind = bind_data.Cast<PvarWriteBindData>();
	auto g = make_uniq<PvarWriteGlobalState>();

	auto &fs = FileSystem::GetFileSystem(context);
	g->path = bind.path_override.empty() ? file_path : bind.path_override;
	// FILE_CREATE_NEW maps to OverwriteExistingFile() -> O_CREAT | O_TRUNC in
	// LocalFileSystem (not O_EXCL), so re-running overwrites cleanly.
	g->handle = fs.OpenFile(g->path, FileFlags::FILE_FLAGS_WRITE | FileFlags::FILE_FLAGS_FILE_CREATE_NEW);

	// Header line: "#CHROM\tPOS\t...\n". Written once here (global init) so even
	// a zero-row query produces a valid header-only `.pvar`.
	string header = "#";
	header += StringUtil::Join(bind.out_headers, "\t");
	header += "\n";
	g->handle->Write(const_cast<char *>(header.data()), header.size());

	return std::move(g);
}

unique_ptr<LocalFunctionData> PvarCopyInitLocal(ExecutionContext &, FunctionData &) {
	return make_uniq<PvarWriteLocalState>();
}

void PvarCopySink(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate, LocalFunctionData &,
                  DataChunk &input) {
	auto &bind = bind_data.Cast<PvarWriteBindData>();
	auto &g = gstate.Cast<PvarWriteGlobalState>();
	auto &fs = FileSystem::GetFileSystem(context.client);

	try {
		input.Flatten();
		string buffer;
		for (idx_t row = 0; row < input.size(); row++) {
			// --- (CHROM, POS) monotonicity ----------------------------------
			// A CopyFunction cannot force ORDER BY on its input, and a round-trip
			// test is blind to disorder (the `.pvar` stays internally consistent).
			// PLINK tooling + this extension's `region` pushdown assume position-
			// sorted, contiguous-per-chromosome variants, so enforce it here.
			auto chrom_val = input.GetValue(bind.chrom_input_idx, row);
			auto pos_val = input.GetValue(bind.pos_input_idx, row);
			if (!chrom_val.IsNull() && !pos_val.IsNull()) {
				auto chrom = chrom_val.ToString();
				auto pos = pos_val.GetValue<int64_t>();
				if (g.has_last && chrom == g.last_chrom) {
					if (pos < g.last_pos) {
						throw InvalidInputException(
						    "write pvar: input is not position-sorted — POS %lld follows %lld on CHROM '%s'. "
						    "Add ORDER BY CHROM, POS to the query (a `.pvar` must be position-sorted so PLINK "
						    "tooling and region pushdown work).",
						    static_cast<long long>(pos), static_cast<long long>(g.last_pos), chrom);
					}
				} else if (g.has_last && chrom != g.last_chrom) {
					// Chromosome changed: it must not be one we already left
					// (interleaved blocks break the contiguous-per-chromosome
					// layout that region pushdown relies on).
					if (g.seen_chroms.count(chrom) != 0) {
						throw InvalidInputException(
						    "write pvar: CHROM '%s' reappears after another chromosome — variants must be grouped "
						    "contiguously per chromosome. Add ORDER BY CHROM, POS to the query.",
						    chrom);
					}
				}
				g.seen_chroms.insert(chrom);
				g.last_chrom = chrom;
				g.last_pos = pos;
				g.has_last = true;
			}

			// --- emit the row -----------------------------------------------
			for (idx_t c = 0; c < bind.out_source_idx.size(); c++) {
				if (c > 0) {
					buffer += '\t';
				}
				buffer += CellToField(input.GetValue(bind.out_source_idx[c], row), bind.out_headers[c]);
			}
			buffer += '\n';
		}
		if (!buffer.empty()) {
			g.handle->Write(const_cast<char *>(buffer.data()), buffer.size());
		}
	} catch (...) {
		// Remove the partial file we created before propagating — DuckDB only
		// auto-removes its single COPY target, not a `pvar` override.
		g.CleanupPartial(fs);
		throw;
	}
}

void PvarCopyCombine(ExecutionContext &, FunctionData &, GlobalFunctionData &, LocalFunctionData &) {
	// No per-thread buffering: REGULAR_COPY_TO_FILE writes straight to the global
	// handle in the sink. Hook is required (DuckDB calls it); nothing to do.
}

void PvarCopyFinalize(ClientContext &context, FunctionData &, GlobalFunctionData &gstate) {
	auto &g = gstate.Cast<PvarWriteGlobalState>();
	if (g.handle) {
		g.handle->Sync();
		g.handle.reset();
	}
	(void)context;
}

CopyFunctionExecutionMode PvarCopyExecMode(bool /*preserve_insertion_order*/, bool /*supports_batch_index*/) {
	// Single-threaded, ordered, no temp file — required for the ordered sink and
	// the (CHROM, POS) monotonicity check to see rows in query order.
	return CopyFunctionExecutionMode::REGULAR_COPY_TO_FILE;
}

} // namespace

void RegisterPfileWriter(ExtensionLoader &loader) {
	CopyFunction fn("pvar");
	fn.copy_to_bind = PvarCopyBind;
	fn.copy_to_initialize_local = PvarCopyInitLocal;
	fn.copy_to_initialize_global = PvarCopyInitGlobal;
	fn.copy_to_sink = PvarCopySink;
	fn.copy_to_combine = PvarCopyCombine;
	fn.copy_to_finalize = PvarCopyFinalize;
	fn.execution_mode = PvarCopyExecMode;
	fn.extension = "pvar";
	loader.RegisterFunction(fn);
}

} // namespace duckdb
