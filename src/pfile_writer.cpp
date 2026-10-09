#include "pfile_writer.hpp"
#include "plink_common.hpp" // AlignedBuffer + pgenlib read/misc headers

#include "duckdb/common/file_system.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/function/copy_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <pgenlib_write.h> // STPgenWriter write path (P9 Phase B, FORMAT pgen)

#include <cstring>
#include <limits>
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

// Core `.pvar` bind: resolve variant roles + carried-through extras against the
// input schema. `exclude` names input column indices that must NOT be carried
// through as extra `.pvar` columns (the bundle passes the genotypes column here
// so it is never stringified into the `.pvar`). Shared by the standalone
// `FORMAT pvar` CopyFunction and the `FORMAT pfile` bundle.
unique_ptr<PvarWriteBindData> BuildPvarBind(CopyFunctionBindInput &input, const vector<string> &names,
                                            const vector<LogicalType> &sql_types, const std::set<idx_t> &exclude) {
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
	// in input order (mirrors read_pvar's "extra columns" output). `exclude`
	// columns (e.g. the bundle's genotypes column) are never carried through.
	for (idx_t i = 0; i < names.size(); i++) {
		if (consumed.count(i) == 0 && exclude.count(i) == 0) {
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

unique_ptr<FunctionData> PvarCopyBind(ClientContext &, CopyFunctionBindInput &input, const vector<string> &names,
                                      const vector<LogicalType> &sql_types) {
	return BuildPvarBind(input, names, sql_types, /*exclude=*/{});
}

// Open the `.pvar` handle at `path` (overwrite) and write its header line. Shared
// by the standalone InitGlobal and the bundle. `path` is the resolved file path.
void PvarOpenAndWriteHeader(PvarWriteGlobalState &g, FileSystem &fs, const PvarWriteBindData &bind, const string &path) {
	g.path = path;
	// FILE_CREATE_NEW maps to OverwriteExistingFile() -> O_CREAT | O_TRUNC in
	// LocalFileSystem (not O_EXCL), so re-running overwrites cleanly.
	g.handle = fs.OpenFile(g.path, FileFlags::FILE_FLAGS_WRITE | FileFlags::FILE_FLAGS_FILE_CREATE_NEW);

	// Header line: "#CHROM\tPOS\t...\n". Written once here (global init) so even
	// a zero-row query produces a valid header-only `.pvar`.
	string header = "#";
	header += StringUtil::Join(bind.out_headers, "\t");
	header += "\n";
	g.handle->Write(const_cast<char *>(header.data()), header.size());
}

unique_ptr<GlobalFunctionData> PvarCopyInitGlobal(ClientContext &context, FunctionData &bind_data,
                                                  const string &file_path) {
	auto &bind = bind_data.Cast<PvarWriteBindData>();
	auto g = make_uniq<PvarWriteGlobalState>();
	auto &fs = FileSystem::GetFileSystem(context);
	PvarOpenAndWriteHeader(*g, fs, bind, bind.path_override.empty() ? file_path : bind.path_override);
	return std::move(g);
}

unique_ptr<LocalFunctionData> PvarCopyInitLocal(ExecutionContext &, FunctionData &) {
	return make_uniq<PvarWriteLocalState>();
}

// Write one chunk of variant rows to the already-open `.pvar` handle, enforcing
// the (CHROM, POS) monotonicity contract. Throws on disorder / unrepresentable
// cell; the CALLER owns cleanup-on-failure (so the bundle can remove all
// siblings, not just the `.pvar`). Shared by the standalone sink and the bundle.
void PvarWriteChunk(const PvarWriteBindData &bind, PvarWriteGlobalState &g, DataChunk &input) {
	{
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
	}
}

void PvarCopySink(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate, LocalFunctionData &,
                  DataChunk &input) {
	auto &bind = bind_data.Cast<PvarWriteBindData>();
	auto &g = gstate.Cast<PvarWriteGlobalState>();
	auto &fs = FileSystem::GetFileSystem(context.client);
	try {
		PvarWriteChunk(bind, g, input);
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

// ===========================================================================
// FORMAT pgen — binary `.pgen` writer (P9 Phase B, biallelic HARDCALLS only)
//
// Writes a PLINK 2 `.pgen` from a variant-oriented query using pgenlib's
// STPgenWriter append path, re-readable by `read_pgen`/`read_pfile` to the same
// decoded int8 hardcall matrix. Companions (`.pvar`/`.psam`) are NOT written
// here (that is FORMAT pvar / the pfile bundle) — only the `.pgen` + its
// genotype payload. See docs/planning/P9-001-writers-design (Decisions 3, 4, 6).
//
// Surface (space-separated `KEY value`, like FORMAT pvar):
//   COPY (SELECT ... genotypes ... FROM read_pfile('x') ORDER BY CHROM, POS)
//   TO 'out.pgen' (FORMAT pgen
//     [, genotype_fields 'genotypes']   -- input genotypes column (default: 'genotypes')
//     [, variant_count 1000000]         -- upper bound on #variants (see below)
//     [, pgen 'override/path.pgen']);   -- file-PATH override (reader symmetry)
//
// Input contract: one row per variant; the genotypes column is a LIST/ARRAY of
// per-sample int8 hardcalls (`0`/`1`/`2`, and `-9` or SQL NULL = missing),
// positional over samples. The sample count is the genotype list length and
// MUST be identical on every row. This is exactly what `read_pfile`'s default
// variant orient emits (`genotypes TINYINT[N]`, missing -> NULL).
//
// Scope: biallelic diploid hardcalls ONLY. Phased (`TINYINT[2][N]`), dosage
// (`DOUBLE[N]`), and multiallelic (comma in ALT, best-effort) inputs are
// rejected with an actionable error at bind/sink — they are later sub-phases.
//
// Write mode: kPgenWriteAndCopy. Rationale — it produces a single, fully
// self-contained `.pgen` with an *embedded* index (genotype body staged in a
// sibling `<path>.tmp` and folded in at SpgwFinish), matching the shape of the
// plink2-produced fixtures. The reader opens it with `PgfiInitPhase1(path,
// nullptr, …)` — no `.pgi` sidecar to plumb through. (kPgenWriteSeparateIndex
// would emit a `.pgi` the reader does not look for; kPgenWriteBackwardSeek needs
// the *exact* variant count up front, which a streaming COPY does not know.)
//
// Variant-count bound: SpgwInitPhase1 wants a variant-count *limit* up front and
// allocates an in-memory index proportional to it; AndCopy only needs an UPPER
// bound (actual count finalized at SpgwFinish). pgenlib does NOT bounds-check the
// per-variant append (`vrec_len_buf[vidx*…]` is written unconditionally, and the
// `assert` is compiled out in release), so overrunning the limit is a heap
// overflow, not an error. The sink therefore GUARDS every append against the
// limit and throws cleanly if exceeded, instructing the user to pass
// `variant_count N`. With that guard, a too-small default is a clean actionable
// error — never corruption — so the default below is purely a memory knob.
//
// Local-only: AndCopy writes via raw FILE* (bypasses DuckDB's VFS, which is
// read-only "rb" for pgen). The target MUST be a local path; a remote/scheme
// target is rejected at bind (P9-001 Decision 6 — a mechanism constraint).
// ===========================================================================

// pgenlib pre-allocates ~1.5 bytes per variant of `variant_ct_limit` for the
// on-the-fly index (the vrec-length + vrtype buffers). At the 2^24 default that
// is ~25 MB of address space, but only ~8–10 MB resident in practice (measured
// maxRSS at 4 and 250 samples). 2^24 ≈ 16.7M variants (exome / chip / single-
// chromosome scale); larger datasets pass `variant_count N`. The sink guard
// makes an exceeded bound a clean error, so this is safe to exceed, not a cap.
constexpr uint32_t kDefaultPgenVariantCtLimit = 1u << 24;

// Upper sanity bound on the `variant_count` option: a value above this is almost
// certainly a typo and would reserve multiple GB (~1.5 bytes/variant). 2e9 is
// comfortably above any real single-file whole-genome variant count (~1e9), and
// below UINT32_MAX so the cast to variant_ct_limit is safe.
constexpr uint32_t kMaxPgenVariantCountOption = 2000000000u;

// True iff `path` begins with a URI scheme ("scheme://" where the scheme starts
// with a letter) — e.g. s3://, http://, pathmacro://. A local filename that
// merely contains "://" somewhere (not as a leading scheme) is NOT flagged.
bool HasUriScheme(const string &path) {
	auto pos = path.find("://");
	if (pos == string::npos || pos == 0) {
		return false;
	}
	auto is_alpha = [](char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'); };
	if (!is_alpha(path[0])) {
		return false;
	}
	for (idx_t i = 1; i < pos; i++) {
		char c = path[i];
		if (!is_alpha(c) && !(c >= '0' && c <= '9') && c != '+' && c != '-' && c != '.') {
			return false;
		}
	}
	return true;
}

bool PathIsLocalWritable(const string &path) {
	// raw FILE* target: reject real remote schemes (s3/http/gcs/… via IsRemoteFile)
	// and any other leading scheme:// URL (e.g. pathmacro://). A bare local path —
	// even one whose name happens to contain "://" mid-string — stays writable.
	return !FileSystem::IsRemoteFile(path) && !HasUriScheme(path);
}

struct PgenWriteBindData : public FunctionData {
	idx_t geno_input_idx = 0;    //! input chunk column index of the genotypes LIST/ARRAY/STRUCT
	bool sample_ct_known = false; //! true iff genotypes is a fixed-size ARRAY or STRUCT (N in the type)
	uint32_t sample_ct = 0;      //! sample count when known at bind
	bool geno_is_struct = false; //! true iff genotypes column is a STRUCT (field NAMES are IIDs)
	vector<string> struct_iids;  //! struct field names, in field order (= IIDs) when geno_is_struct
	bool has_alt = false;        //! an ALT-role column is present (multiallelic guard)
	idx_t alt_input_idx = 0;
	uint32_t variant_ct_limit = kDefaultPgenVariantCtLimit;
	string path_override;        //! `pgen` option; empty -> use the COPY target

	unique_ptr<FunctionData> Copy() const override {
		auto r = make_uniq<PgenWriteBindData>();
		r->geno_input_idx = geno_input_idx;
		r->sample_ct_known = sample_ct_known;
		r->sample_ct = sample_ct;
		r->geno_is_struct = geno_is_struct;
		r->struct_iids = struct_iids;
		r->has_alt = has_alt;
		r->alt_input_idx = alt_input_idx;
		r->variant_ct_limit = variant_ct_limit;
		r->path_override = path_override;
		return std::move(r);
	}
	bool Equals(const FunctionData &other_p) const override {
		auto &o = other_p.Cast<PgenWriteBindData>();
		return geno_input_idx == o.geno_input_idx && sample_ct_known == o.sample_ct_known &&
		       sample_ct == o.sample_ct && geno_is_struct == o.geno_is_struct && struct_iids == o.struct_iids &&
		       has_alt == o.has_alt && alt_input_idx == o.alt_input_idx &&
		       variant_ct_limit == o.variant_ct_limit && path_override == o.path_override;
	}
};

struct PgenWriteGlobalState : public GlobalFunctionData {
	string path;                 //! path handed to InitGlobal (what pgenlib fopens); for cleanup
	plink2::STPgenWriter spgw;
	AlignedBuffer spgw_alloc;    //! pgenlib's working allocation
	AlignedBuffer genovec;       //! packed 2-bit genotype buffer (one variant)
	uintptr_t genovec_bytes = 0; //! allocated size of `genovec` (for per-row memset)
	uint32_t sample_ct = 0;
	bool sample_ct_known = false;
	uint32_t variant_ct_limit = 0;
	uint64_t appended = 0;
	bool writer_inited = false;  //! SpgwInitPhase1/2 done (files opened)
	bool finished = false;       //! SpgwFinish done (don't cleanup a completed file)

	//! Close pgenlib + remove the partial `.pgen` and its AndCopy `.tmp` sibling.
	//! Best-effort; never masks the original error.
	void CleanupPartial(FileSystem &fs) {
		if (writer_inited && !finished) {
			plink2::PglErr e = plink2::kPglRetSuccess;
			plink2::CleanupSpgw(&spgw, &e); // closes the FILE*s so the files can be removed
			finished = true;               // don't let the destructor CleanupSpgw a second time
		}
		if (path.empty()) {
			return;
		}
		for (const char *suffix : {"", ".tmp"}) {
			try {
				string p = path + suffix;
				if (fs.FileExists(p)) {
					fs.RemoveFile(p);
				}
			} catch (...) {
				// ignore
			}
		}
	}

	~PgenWriteGlobalState() override {
		// If the query threw past finalize without our catch running (shouldn't,
		// but be safe), close pgenlib's handles so they are not leaked.
		if (writer_inited && !finished) {
			plink2::PglErr e = plink2::kPglRetSuccess;
			plink2::CleanupSpgw(&spgw, &e);
		}
	}
};

struct PgenWriteLocalState : public LocalFunctionData {};

//! SpgwInitPhase1/2 + allocate buffers. `sample_ct` is from the ARRAY type (bind)
//! or the first row's list length (sink). Opens the output file(s) via raw FILE*.
void PgenEnsureWriterInited(PgenWriteGlobalState &g, uint32_t sample_ct) {
	if (g.writer_inited) {
		return;
	}
	if (sample_ct == 0) {
		throw InvalidInputException("write pgen: sample count is 0 — a .pgen must have at least one sample");
	}
	plink2::PreinitSpgw(&g.spgw);
	uintptr_t alloc_cacheline_ct = 0;
	uint32_t max_vrec_len = 0;
	plink2::PglErr e = plink2::SpgwInitPhase1(g.path.c_str(),
	                                          /*allele_idx_offsets=*/nullptr,
	                                          /*explicit_nonref_flags=*/nullptr,
	                                          /*variant_ct_limit=*/g.variant_ct_limit,
	                                          /*sample_ct=*/sample_ct,
	                                          /*allele_ct_upper_bound=*/0, plink2::kPgenWriteAndCopy,
	                                          /*phase_dosage_gflags=*/plink2::kfPgenGlobal0,
	                                          /*nonref_flags_storage=*/1, &g.spgw, &alloc_cacheline_ct, &max_vrec_len);
	if (e != plink2::kPglRetSuccess) {
		throw IOException("write pgen: SpgwInitPhase1 failed for '%s' (pgenlib error %u)", g.path,
		                  static_cast<uint32_t>(e));
	}
	g.writer_inited = true; // files are now open — a later throw must clean them up
	g.spgw_alloc.Allocate(alloc_cacheline_ct * plink2::kCacheline);
	plink2::SpgwInitPhase2(max_vrec_len, &g.spgw, g.spgw_alloc.As<unsigned char>());
	g.genovec_bytes = plink2::NypCtToVecCt(sample_ct) * plink2::kBytesPerVec;
	g.genovec.Allocate(g.genovec_bytes);
	g.sample_ct = sample_ct;
}

// Validate a genotypes child/element type is a biallelic-hardcall int: reject
// PHASED (nested LIST/ARRAY), DOSAGE (DOUBLE/FLOAT), and non-integer. Shared by
// the LIST/ARRAY and STRUCT classification paths.
void ValidateHardcallElem(const LogicalType &child, const string &col) {
	LogicalTypeId cid = child.id();
	if (cid == LogicalTypeId::ARRAY || cid == LogicalTypeId::LIST) {
		throw BinderException("write pgen: genotypes column '%s' looks PHASED (element type %s) — phased .pgen writing "
		                      "is not supported in this increment (biallelic hardcalls only)",
		                      col, child.ToString());
	}
	if (cid == LogicalTypeId::DOUBLE || cid == LogicalTypeId::FLOAT) {
		throw BinderException("write pgen: genotypes column '%s' looks like DOSAGES (element type %s) — dosage .pgen "
		                      "writing is not supported in this increment (biallelic hardcalls only)",
		                      col, child.ToString());
	}
	if (!IsIntegerType(child)) {
		throw BinderException("write pgen: genotypes column '%s' has element type %s, but biallelic hardcalls must be "
		                      "an integer (0/1/2/-9)",
		                      col, child.ToString());
	}
}

// Core `.pgen` bind: resolve + classify the genotypes column, parse variant_count
// / pgen path override, and enforce the local-only target. Shared by the
// standalone `FORMAT pgen` CopyFunction and the `FORMAT pfile` bundle.
unique_ptr<PgenWriteBindData> BuildPgenBind(CopyFunctionBindInput &input, const vector<string> &names,
                                            const vector<LogicalType> &sql_types) {
	auto bind = make_uniq<PgenWriteBindData>();

	string geno_col = "genotypes"; // reader-aligned default

	for (auto &opt : input.info.options) {
		auto key = StringUtil::Lower(opt.first);
		if (key == "genotype_fields") {
			if (opt.second.empty() || opt.second[0].IsNull()) {
				continue;
			}
			auto &val = opt.second[0];
			if (val.type().id() == LogicalTypeId::STRUCT) {
				throw BinderException("write pgen: genotype_fields as a STRUCT (multiple genotype tracks) is not "
				                      "supported in this increment — pass a single column name, e.g. "
				                      "genotype_fields 'genotypes' (biallelic hardcalls only)");
			}
			geno_col = val.ToString();
		} else if (key == "variant_count") {
			if (!opt.second.empty() && !opt.second[0].IsNull()) {
				int64_t vc = opt.second[0].GetValue<int64_t>();
				if (vc <= 0) {
					throw BinderException("write pgen: variant_count must be a positive integer");
				}
				// variant_count is an UPPER BOUND pgenlib pre-allocates against at ~1.5
				// bytes per variant. Reject an implausibly large value (a likely typo)
				// that would reserve multiple GB before a single row is seen.
				if (static_cast<uint64_t>(vc) > kMaxPgenVariantCountOption) {
					throw BinderException(
					    "write pgen: variant_count %lld is implausibly large (> %u) — this looks like a typo and would "
					    "pre-allocate multiple GB (~1.5 bytes per variant of the bound). Pass an upper bound at or below "
					    "%u (above any real single-file whole-genome variant count), or omit variant_count to stream "
					    "against the default bound.",
					    static_cast<long long>(vc), kMaxPgenVariantCountOption, kMaxPgenVariantCountOption);
				}
				bind->variant_ct_limit = static_cast<uint32_t>(vc);
			}
		} else if (key == "pgen") {
			if (!opt.second.empty() && !opt.second[0].IsNull()) {
				bind->path_override = opt.second[0].ToString();
			}
		}
		// Other options (pvar/psam/variant_fields/sample_fields/samples) belong to
		// the text writers / bundle; ignore here (a round-trip SELECT may carry them).
	}

	// --- resolve the genotypes column (case-insensitive) ---
	idx_t gidx = names.size();
	for (idx_t i = 0; i < names.size(); i++) {
		if (StringUtil::CIEquals(names[i], geno_col)) {
			gidx = i;
			break;
		}
	}
	if (gidx == names.size()) {
		throw BinderException("write pgen: genotypes column '%s' is not in the input (available columns: %s). Supply "
		                      "genotype_fields '<col>' to name it.",
		                      geno_col, StringUtil::Join(names, ", "));
	}
	bind->geno_input_idx = gidx;

	// --- classify the genotypes column type; reject non-hardcall shapes ---
	// Accepted shapes: LIST/ARRAY of per-sample int calls (sample order is list
	// order; IIDs out-of-band via `samples`), or a STRUCT whose field NAMES are the
	// sample IIDs and whose values (in field order) are the per-sample int calls
	// (this is read_pfile(..., genotypes:='struct') variant orient).
	const LogicalType &gtype = sql_types[gidx];
	LogicalTypeId tid = gtype.id();
	if (tid == LogicalTypeId::STRUCT) {
		auto &child_types = StructType::GetChildTypes(gtype);
		if (child_types.empty()) {
			throw BinderException("write pgen: genotypes STRUCT '%s' has no fields — a .pgen needs at least one sample",
			                      names[gidx]);
		}
		// Every field must be a hardcall int; phased ARRAY(TINYINT,2) / dosage DOUBLE
		// struct fields reject with the same actionable messages as the LIST path.
		for (auto &ct : child_types) {
			ValidateHardcallElem(ct.second, names[gidx]);
			bind->struct_iids.push_back(ct.first); // field name = IID
		}
		bind->geno_is_struct = true;
		bind->sample_ct = static_cast<uint32_t>(child_types.size());
		bind->sample_ct_known = true;
	} else if (tid == LogicalTypeId::ARRAY || tid == LogicalTypeId::LIST) {
		const LogicalType &child = (tid == LogicalTypeId::ARRAY) ? ArrayType::GetChildType(gtype)
		                                                         : ListType::GetChildType(gtype);
		ValidateHardcallElem(child, names[gidx]);
		if (tid == LogicalTypeId::ARRAY) {
			bind->sample_ct = static_cast<uint32_t>(ArrayType::GetSize(gtype));
			bind->sample_ct_known = true;
		}
	} else {
		throw BinderException("write pgen: genotypes column '%s' has type %s, but a biallelic-hardcall .pgen needs a "
		                      "LIST or ARRAY of per-sample calls (e.g. TINYINT[N]), or a STRUCT keyed by sample IID",
		                      names[gidx], gtype.ToString());
	}

	// --- optional ALT role for the best-effort multiallelic guard ---
	for (idx_t i = 0; i < names.size(); i++) {
		if (StringUtil::CIEquals(names[i], "ALT")) {
			bind->has_alt = true;
			bind->alt_input_idx = i;
			break;
		}
	}

	// --- local-only target (raw FILE*; no VFS) — reject remote/scheme at bind ---
	const string &target = bind->path_override.empty() ? input.info.file_path : bind->path_override;
	if (!target.empty() && !PathIsLocalWritable(target)) {
		throw BinderException("write pgen: '%s' is not a local path. A .pgen is written through a raw file handle that "
		                      "bypasses DuckDB's VFS, so the target must be local (no s3://, http://, or pathmacro:// "
		                      "targets). Write to a local path, or write the .pvar/.psam companions to the remote "
		                      "store separately.",
		                      target);
	}

	return std::move(bind);
}

unique_ptr<FunctionData> PgenCopyBind(ClientContext &, CopyFunctionBindInput &input, const vector<string> &names,
                                      const vector<LogicalType> &sql_types) {
	return BuildPgenBind(input, names, sql_types);
}

unique_ptr<GlobalFunctionData> PgenCopyInitGlobal(ClientContext &, FunctionData &bind_data, const string &file_path) {
	auto &bind = bind_data.Cast<PgenWriteBindData>();
	auto g = make_uniq<PgenWriteGlobalState>();
	g->path = bind.path_override.empty() ? file_path : bind.path_override;
	g->sample_ct_known = bind.sample_ct_known;
	g->sample_ct = bind.sample_ct;
	g->variant_ct_limit = bind.variant_ct_limit;

	// Writer init (SpgwInitPhase1/2, which fopens the .pgen + its .tmp sibling) is
	// DEFERRED to the first sink row for BOTH input shapes — even when sample_ct
	// is already known from a fixed-size ARRAY type. Two reasons:
	//   1. A zero-row query then never opens a file, and Finalize fails cleanly
	//      (pgenlib cannot finish a 0-variant .pgen — SpgwFinish's PwcFinish
	//      asserts variant_ct, which is compiled out in release and then SIGSEGVs).
	//   2. An init failure happens INSIDE the sink's try-guard, so CleanupPartial
	//      removes the partial .pgen/.tmp instead of orphaning them (InitGlobal has
	//      no try/catch).
	// sample_ct_known/sample_ct are carried in the global state and used to size
	// the writer on the first row (so a wrong first-row length is still caught).
	return std::move(g);
}

unique_ptr<LocalFunctionData> PgenCopyInitLocal(ExecutionContext &, FunctionData &) {
	return make_uniq<PgenWriteLocalState>();
}

// Append one chunk of variant rows to the `.pgen` (deferred writer init on the
// first row). Throws on a bad shape/value/overrun; the CALLER owns cleanup (so
// the bundle removes all siblings). Shared by the standalone sink and the bundle.
void PgenWriteChunk(const PgenWriteBindData &bind, PgenWriteGlobalState &g, DataChunk &input) {
	{
		input.Flatten();
		for (idx_t row = 0; row < input.size(); row++) {
			// --- best-effort multiallelic guard (only when ALT is in the stream) ---
			// A genotypes-only input carries NO multiallelic signal (PgrGet collapses
			// a multiallelic site to biallelic 0/1/2/missing). When the natural
			// round-trip query carries ALT, a comma reveals it; reject those.
			if (bind.has_alt) {
				auto alt_val = input.GetValue(bind.alt_input_idx, row);
				if (!alt_val.IsNull() && alt_val.ToString().find(',') != string::npos) {
					throw InvalidInputException(
					    "write pgen: variant has a multiallelic ALT ('%s') — multiallelic .pgen writing is not "
					    "supported in this increment (biallelic hardcalls only)",
					    alt_val.ToString());
				}
			}

			// --- fetch the per-variant genotype list ---
			auto geno_val = input.GetValue(bind.geno_input_idx, row);
			if (geno_val.IsNull()) {
				throw InvalidInputException(
				    "write pgen: the genotypes cell is NULL for a variant — every variant needs a list of per-sample "
				    "calls (use -9 or NULL *inside* the list for a missing call)");
			}
			// LIST/ARRAY children, or STRUCT field values in field order (same
			// positional per-sample semantics — the struct's field NAMES carried the
			// IIDs to the .psam at bind, the VALUES are the per-sample calls here).
			LogicalTypeId gtid = geno_val.type().id();
			const vector<Value> &calls = (gtid == LogicalTypeId::ARRAY)    ? ArrayValue::GetChildren(geno_val)
			                             : (gtid == LogicalTypeId::STRUCT) ? StructValue::GetChildren(geno_val)
			                                                               : ListValue::GetChildren(geno_val);
			idx_t n = calls.size();

			// --- establish / validate the sample count ---
			// Use the bind-known sample_ct (fixed ARRAY width) when we have it, else
			// take it from this first row's list length (LIST input).
			if (!g.writer_inited) {
				uint32_t sc = g.sample_ct_known ? g.sample_ct : static_cast<uint32_t>(n);
				PgenEnsureWriterInited(g, sc);
			}
			if (n != g.sample_ct) {
				throw InvalidInputException(
				    "write pgen: variant %llu has %llu genotypes but the sample count is %u — every variant must list "
				    "the same number of samples",
				    static_cast<unsigned long long>(g.appended), static_cast<unsigned long long>(n), g.sample_ct);
			}

			// --- guard against exceeding the variant-count bound ---
			// pgenlib does NOT bounds-check the append; overrunning corrupts the heap.
			if (g.appended >= g.variant_ct_limit) {
				throw InvalidInputException(
				    "write pgen: the number of variants exceeds variant_count=%u. Pass `variant_count N` with an upper "
				    "bound >= the number of variants being written.",
				    g.variant_ct_limit);
			}

			// --- pack the int8 calls into a 2-bit genovec (sample 0 in low bits) ---
			//   0 -> 0b00, 1 -> 0b01, 2 -> 0b10, (-9 | NULL) -> 0b11 (missing).
			auto *genovec = g.genovec.As<uintptr_t>();
			std::memset(genovec, 0, g.genovec_bytes);
			for (idx_t s = 0; s < n; s++) {
				uintptr_t nyp;
				const Value &c = calls[s];
				if (c.IsNull()) {
					nyp = 3; // missing
				} else {
					int64_t code = c.GetValue<int64_t>();
					switch (code) {
					case 0:
						nyp = 0;
						break;
					case 1:
						nyp = 1;
						break;
					case 2:
						nyp = 2;
						break;
					case -9:
						nyp = 3;
						break;
					default:
						throw InvalidInputException(
						    "write pgen: genotype value %lld (variant %llu, sample %llu) is out of range — a biallelic "
						    "hardcall must be 0, 1, 2, or -9/NULL (missing)",
						    static_cast<long long>(code), static_cast<unsigned long long>(g.appended),
						    static_cast<unsigned long long>(s));
					}
				}
				uint32_t widx = static_cast<uint32_t>(s) / plink2::kBitsPerWordD2;
				uint32_t shift = 2 * (static_cast<uint32_t>(s) % plink2::kBitsPerWordD2);
				genovec[widx] |= nyp << shift;
			}
			plink2::ZeroTrailingBits(g.sample_ct * 2, genovec);

			plink2::PglErr e = plink2::SpgwAppendBiallelicGenovec(genovec, &g.spgw);
			if (e != plink2::kPglRetSuccess) {
				throw IOException("write pgen: SpgwAppendBiallelicGenovec failed at variant %llu (pgenlib error %u)",
				                  static_cast<unsigned long long>(g.appended), static_cast<uint32_t>(e));
			}
			g.appended++;
		}
	}
}

// Finalize the `.pgen` (SpgwFinish + CleanupSpgw). Throws on zero-variant input
// (pgenlib cannot finish a 0-variant .pgen) or an SpgwFinish error; the CALLER
// owns cleanup-on-failure. Sets g.finished so the destructor won't double-close.
void PgenFinish(PgenWriteGlobalState &g) {
	if (g.finished) {
		return;
	}
	if (!g.writer_inited || g.appended == 0) {
		// Zero variants. Init is deferred to the first sink row, so normally
		// !writer_inited here (no file opened). The appended==0 arm is defensive:
		// never feed vidx=0 to SpgwFinish — pgenlib asserts variant_ct>0 (compiled
		// out in release) and then SIGSEGVs.
		throw IOException("write pgen: cannot write a .pgen with zero variants — the input query produced no rows");
	}
	plink2::PglErr e = plink2::SpgwFinish(&g.spgw);
	if (e != plink2::kPglRetSuccess) {
		throw IOException("write pgen: SpgwFinish failed (pgenlib error %u)", static_cast<uint32_t>(e));
	}
	plink2::PglErr ce = plink2::kPglRetSuccess;
	plink2::CleanupSpgw(&g.spgw, &ce);
	g.finished = true;
}

void PgenCopySink(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate, LocalFunctionData &,
                  DataChunk &input) {
	auto &bind = bind_data.Cast<PgenWriteBindData>();
	auto &g = gstate.Cast<PgenWriteGlobalState>();
	auto &fs = FileSystem::GetFileSystem(context.client);
	try {
		PgenWriteChunk(bind, g, input);
	} catch (...) {
		g.CleanupPartial(fs);
		throw;
	}
}

void PgenCopyCombine(ExecutionContext &, FunctionData &, GlobalFunctionData &, LocalFunctionData &) {
	// REGULAR_COPY_TO_FILE: the sink writes straight to the global writer. Nothing
	// to merge; the hook is required but a no-op.
}

void PgenCopyFinalize(ClientContext &context, FunctionData &, GlobalFunctionData &gstate) {
	auto &g = gstate.Cast<PgenWriteGlobalState>();
	auto &fs = FileSystem::GetFileSystem(context);
	try {
		PgenFinish(g);
	} catch (...) {
		// Remove any partial .pgen/.tmp we opened (zero-variant case opened nothing,
		// but a mid-finish failure must not orphan the AndCopy temp).
		g.CleanupPartial(fs);
		throw;
	}
}

CopyFunctionExecutionMode PgenCopyExecMode(bool /*preserve_insertion_order*/, bool /*supports_batch_index*/) {
	// Single-threaded, ordered: .pgen records are strictly positional and the
	// single STPgenWriter appends in query row order.
	return CopyFunctionExecutionMode::REGULAR_COPY_TO_FILE;
}

// ===========================================================================
// FORMAT psam — text `.psam` writer + THE SAMPLE-IDENTITY CONTRACT (P9 Phase C)
//
// Writes a MINIMAL `#IID`-only `.psam` (one IID per sample), re-readable by
// `read_psam`/`read_pfile`. Per-sample metadata beyond IID (SEX/phenotype/
// covariates) is DEFERRED — a variant-oriented stream has no per-sample axis to
// carry it. See docs/planning/P9-001-writers-design (Decision 2) + the task's
// SAMPLE-IDENTITY CONTRACT.
//
// IID source, in priority order (NEVER fabricate IIDs):
//   1. a STRUCT genotypes column — its field NAMES are the sample IIDs
//      (read_pfile(..., genotypes:='struct') variant orient). OUT-OF-BAND: the
//      whole `.psam` is written once (at init), independent of the row stream.
//   2. a `samples ['id1', …]` option — an explicit IID list. OUT-OF-BAND.
//   3. (standalone psam only) a per-row IID column (default `IID`, remappable via
//      `sample_fields {iid: '<col>'}`) — ROW-STREAM: one `.psam` line per input row
//      (a sample-oriented input, e.g. `SELECT IID FROM read_psam('x')`).
// If NONE is present for the bundle, a clean BinderException (no files written).
//
// `sample_fields {sex: …, …}` (metadata beyond IID) → explicit "not yet
// supported (follow-on)" error; we NEVER silently drop requested metadata.
// ===========================================================================

// Reject an IID that cannot be represented in a tab-delimited `.psam` (or that
// read_psam would reject). Empty / tab / newline / CR are all unrepresentable.
void ValidatePsamIid(const string &iid) {
	if (iid.empty()) {
		throw InvalidInputException("write psam: an IID is empty — .psam IIDs must be non-empty");
	}
	if (iid.find('\t') != string::npos || iid.find('\n') != string::npos || iid.find('\r') != string::npos) {
		throw InvalidInputException("write psam: IID '%s' contains a tab, newline, or carriage return, which cannot be "
		                            "written to a tab-delimited .psam",
		                            iid);
	}
}

// Write a minimal `#IID`-only `.psam` from a fixed IID list (out-of-band identity).
// Validates: non-empty list, each IID representable, no duplicates (read_psam
// rejects duplicate IIDs, so a dup here would write an unreadable file). Syncs +
// closes so the file is durable before any later sibling failure.
void WritePsamFromIids(FileSystem &fs, const string &path, const vector<string> &iids) {
	if (iids.empty()) {
		throw InvalidInputException("write psam: no samples — a .psam needs at least one IID");
	}
	std::set<string> seen;
	string buf = "#IID\n";
	for (auto &iid : iids) {
		ValidatePsamIid(iid);
		if (!seen.insert(iid).second) {
			throw InvalidInputException("write psam: duplicate IID '%s' — .psam IIDs must be unique", iid);
		}
		buf += iid;
		buf += '\n';
	}
	auto h = fs.OpenFile(path, FileFlags::FILE_FLAGS_WRITE | FileFlags::FILE_FLAGS_FILE_CREATE_NEW);
	h->Write(const_cast<char *>(buf.data()), buf.size());
	h->Sync();
	h.reset();
}

// Parse the `samples ['a','b', …]` option. Returns true iff the option is present
// (even if empty); fills `out` with the IIDs. A LIST literal arrives as a single
// LIST Value; be tolerant of a scalar-per-entry spelling too.
bool TryParseSamplesOption(CopyFunctionBindInput &input, vector<string> &out, bool &present) {
	present = false;
	for (auto &opt : input.info.options) {
		if (StringUtil::Lower(opt.first) != "samples") {
			continue;
		}
		present = true;
		if (opt.second.empty()) {
			return true;
		}
		auto &v0 = opt.second[0];
		if (v0.type().id() == LogicalTypeId::LIST) {
			for (auto &c : ListValue::GetChildren(v0)) {
				if (c.IsNull()) {
					throw BinderException("write: the `samples` list contains a NULL IID");
				}
				out.push_back(c.ToString());
			}
		} else {
			for (auto &v : opt.second) {
				if (v.IsNull()) {
					throw BinderException("write: the `samples` list contains a NULL IID");
				}
				out.push_back(v.ToString());
			}
		}
		return true;
	}
	return false;
}

// Inspect a `sample_fields {iid: '<col>', …}` option. Sets `iid_col`/`has_iid` when
// an `iid` mapping is present; THROWS on any other key (metadata beyond IID is
// deferred — never silently dropped). No option present leaves outputs untouched.
void CheckSampleFieldsIid(CopyFunctionBindInput &input, string &iid_col, bool &has_iid) {
	for (auto &opt : input.info.options) {
		if (StringUtil::Lower(opt.first) != "sample_fields") {
			continue;
		}
		if (opt.second.empty() || opt.second[0].IsNull()) {
			continue;
		}
		auto &val = opt.second[0];
		if (val.type().id() != LogicalTypeId::STRUCT) {
			throw BinderException("write: sample_fields must be a STRUCT of role -> column name, e.g. {iid: 'IID'}");
		}
		auto &child_types = StructType::GetChildTypes(val.type());
		auto &children = StructValue::GetChildren(val);
		for (idx_t c = 0; c < child_types.size(); c++) {
			auto role = StringUtil::Lower(child_types[c].first);
			if (role == "iid") {
				if (children[c].IsNull()) {
					throw BinderException("write: sample_fields role 'iid' has a NULL column name");
				}
				iid_col = children[c].ToString();
				has_iid = true;
			} else {
				throw BinderException("write: sample metadata beyond IID is not yet supported (follow-on): "
				                      "sample_fields role '%s'. For now .psam is IID-only.",
				                      child_types[c].first);
			}
		}
	}
}

// Find a STRUCT genotypes column (default name `genotypes`, or the
// `genotype_fields '<col>'` override) in the input schema; returns its field names
// (the IIDs) if found. Used to source out-of-band identity from a
// genotypes:='struct' stream. Returns false if the column is absent or not a STRUCT.
bool TryStructIidsFromSchema(CopyFunctionBindInput &input, const vector<string> &names,
                             const vector<LogicalType> &sql_types, vector<string> &out) {
	string geno_col = "genotypes";
	for (auto &opt : input.info.options) {
		if (StringUtil::Lower(opt.first) == "genotype_fields" && !opt.second.empty() && !opt.second[0].IsNull() &&
		    opt.second[0].type().id() != LogicalTypeId::STRUCT) {
			geno_col = opt.second[0].ToString();
		}
	}
	for (idx_t i = 0; i < names.size(); i++) {
		if (StringUtil::CIEquals(names[i], geno_col) && sql_types[i].id() == LogicalTypeId::STRUCT) {
			for (auto &ct : StructType::GetChildTypes(sql_types[i])) {
				out.push_back(ct.first);
			}
			return true;
		}
	}
	return false;
}

// Resolve the out-of-band IID list for the pfile bundle (and the out-of-band psam
// modes), enforcing the SAMPLE-IDENTITY CONTRACT. `struct_iids` is non-empty iff
// the genotypes column is a STRUCT. Precedence: struct field names, then a
// `samples` list; when BOTH are present they must AGREE (so an explicit `samples`
// is never silently ignored). Absent both → clean BinderException.
vector<string> ResolveIdentity(CopyFunctionBindInput &input, const vector<string> &struct_iids, bool is_struct) {
	vector<string> samples;
	bool samples_present = false;
	TryParseSamplesOption(input, samples, samples_present);

	vector<string> iids;
	if (is_struct) {
		iids = struct_iids;
		if (samples_present) {
			if (samples != struct_iids) {
				throw BinderException(
				    "write pfile: the `samples` list does not match the STRUCT genotypes field names (IIDs). "
				    "Omit `samples` to use the struct field names, or make them agree (same IIDs, same order).");
			}
		}
	} else if (samples_present) {
		iids = samples;
	} else {
		throw BinderException(
		    "FORMAT pfile needs sample identity to write .psam: read with genotypes:='struct' (the struct field "
		    "names become the IIDs), or pass samples ['id1', …]. (A plain TINYINT[N] genotypes LIST carries no "
		    "sample identity, so no .psam can be written without one of these.)");
	}
	// Validate up front so a bad IID fails at bind, before any file is opened.
	std::set<string> seen;
	for (auto &iid : iids) {
		ValidatePsamIid(iid);
		if (!seen.insert(iid).second) {
			throw BinderException("write pfile: duplicate IID '%s' — sample IIDs must be unique", iid);
		}
	}
	return iids;
}

// --- standalone FORMAT psam ------------------------------------------------

enum class PsamWriteMode { OUT_OF_BAND, ROW_STREAM };

struct PsamWriteBindData : public FunctionData {
	PsamWriteMode mode = PsamWriteMode::ROW_STREAM;
	vector<string> iids;      //! OUT_OF_BAND: the full IID list (written once at init)
	idx_t iid_input_idx = 0;  //! ROW_STREAM: input column index of the per-row IID
	string path_override;     //! `psam` option; empty -> use the COPY target

	unique_ptr<FunctionData> Copy() const override {
		auto r = make_uniq<PsamWriteBindData>();
		r->mode = mode;
		r->iids = iids;
		r->iid_input_idx = iid_input_idx;
		r->path_override = path_override;
		return std::move(r);
	}
	bool Equals(const FunctionData &other_p) const override {
		auto &o = other_p.Cast<PsamWriteBindData>();
		return mode == o.mode && iids == o.iids && iid_input_idx == o.iid_input_idx &&
		       path_override == o.path_override;
	}
};

struct PsamWriteGlobalState : public GlobalFunctionData {
	unique_ptr<FileHandle> handle;
	string path;
	std::set<string> seen_iids; //! ROW_STREAM dup-IID guard (read_psam rejects dups)

	void CleanupPartial(FileSystem &fs) {
		if (handle) {
			handle.reset();
		}
		try {
			if (!path.empty() && fs.FileExists(path)) {
				fs.RemoveFile(path);
			}
		} catch (...) {
		}
	}
};

struct PsamWriteLocalState : public LocalFunctionData {};

unique_ptr<FunctionData> PsamCopyBind(ClientContext &, CopyFunctionBindInput &input, const vector<string> &names,
                                      const vector<LogicalType> &sql_types) {
	auto bind = make_uniq<PsamWriteBindData>();

	for (auto &opt : input.info.options) {
		if (StringUtil::Lower(opt.first) == "psam" && !opt.second.empty() && !opt.second[0].IsNull()) {
			bind->path_override = opt.second[0].ToString();
		}
	}

	// IID column remap (standalone row-stream) + reject metadata-beyond-IID.
	string iid_col = "IID";
	bool has_iid_map = false;
	CheckSampleFieldsIid(input, iid_col, has_iid_map);

	// Priority: STRUCT genotypes field names > `samples` list > per-row IID column.
	vector<string> struct_iids;
	bool is_struct = TryStructIidsFromSchema(input, names, sql_types, struct_iids);
	vector<string> samples;
	bool samples_present = false;
	TryParseSamplesOption(input, samples, samples_present);

	if (is_struct || samples_present) {
		bind->mode = PsamWriteMode::OUT_OF_BAND;
		bind->iids = ResolveIdentity(input, struct_iids, is_struct);
	} else {
		// ROW_STREAM: resolve the IID column (default IID, or sample_fields {iid:'col'}).
		bind->mode = PsamWriteMode::ROW_STREAM;
		idx_t idx = names.size();
		for (idx_t i = 0; i < names.size(); i++) {
			if (StringUtil::CIEquals(names[i], iid_col)) {
				idx = i;
				break;
			}
		}
		if (idx == names.size()) {
			throw BinderException("write psam: no IID source — the input has no '%s' column, and neither a STRUCT "
			                      "genotypes column nor a `samples [...]` option was given (available columns: %s). "
			                      "Pass samples ['id1', …], read with genotypes:='struct', or provide an IID column "
			                      "(remap with sample_fields {iid: '<col>'}).",
			                      iid_col, StringUtil::Join(names, ", "));
		}
		bind->iid_input_idx = idx;
	}
	return std::move(bind);
}

unique_ptr<GlobalFunctionData> PsamCopyInitGlobal(ClientContext &context, FunctionData &bind_data,
                                                  const string &file_path) {
	auto &bind = bind_data.Cast<PsamWriteBindData>();
	auto g = make_uniq<PsamWriteGlobalState>();
	auto &fs = FileSystem::GetFileSystem(context);
	g->path = bind.path_override.empty() ? file_path : bind.path_override;

	if (bind.mode == PsamWriteMode::OUT_OF_BAND) {
		// Identity is fully known at bind — write the whole .psam now; the sink is a
		// no-op. (A zero-row variant stream still produces the correct .psam.)
		WritePsamFromIids(fs, g->path, bind.iids);
	} else {
		// ROW_STREAM: open + header now; one line per input row in the sink.
		g->handle = fs.OpenFile(g->path, FileFlags::FILE_FLAGS_WRITE | FileFlags::FILE_FLAGS_FILE_CREATE_NEW);
		const char *hdr = "#IID\n";
		g->handle->Write(const_cast<char *>(hdr), 5);
	}
	return std::move(g);
}

unique_ptr<LocalFunctionData> PsamCopyInitLocal(ExecutionContext &, FunctionData &) {
	return make_uniq<PsamWriteLocalState>();
}

void PsamCopySink(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate, LocalFunctionData &,
                  DataChunk &input) {
	auto &bind = bind_data.Cast<PsamWriteBindData>();
	auto &g = gstate.Cast<PsamWriteGlobalState>();
	auto &fs = FileSystem::GetFileSystem(context.client);
	if (bind.mode == PsamWriteMode::OUT_OF_BAND) {
		return; // .psam already fully written at init
	}
	try {
		input.Flatten();
		string buffer;
		for (idx_t row = 0; row < input.size(); row++) {
			auto v = input.GetValue(bind.iid_input_idx, row);
			if (v.IsNull()) {
				throw InvalidInputException("write psam: the IID column is NULL for a row — every sample needs an IID");
			}
			auto iid = v.ToString();
			ValidatePsamIid(iid);
			if (!g.seen_iids.insert(iid).second) {
				throw InvalidInputException("write psam: duplicate IID '%s' — .psam IIDs must be unique", iid);
			}
			buffer += iid;
			buffer += '\n';
		}
		if (!buffer.empty()) {
			g.handle->Write(const_cast<char *>(buffer.data()), buffer.size());
		}
	} catch (...) {
		g.CleanupPartial(fs);
		throw;
	}
}

void PsamCopyCombine(ExecutionContext &, FunctionData &, GlobalFunctionData &, LocalFunctionData &) {
}

void PsamCopyFinalize(ClientContext &context, FunctionData &, GlobalFunctionData &gstate) {
	auto &g = gstate.Cast<PsamWriteGlobalState>();
	if (g.handle) {
		g.handle->Sync();
		g.handle.reset();
	}
	(void)context;
}

CopyFunctionExecutionMode PsamCopyExecMode(bool, bool) {
	return CopyFunctionExecutionMode::REGULAR_COPY_TO_FILE;
}

// ===========================================================================
// FORMAT pfile — the bundle: one COPY → prefix.pvar + prefix.pgen + prefix.psam
//
//   COPY (SELECT … ORDER BY CHROM, POS) TO 'prefix' (FORMAT pfile
//     [, variant_fields {…}] [, genotype_fields '<col>']
//     [, samples ['id1', …]] [, variant_count N]
//     [, pvar '…'] [, pgen '…'] [, psam '…']);   -- path overrides (reader symmetry)
//
// Fans out to the three writers above: pvar (text), pgen (pgenlib, local-only),
// psam (IID-only, out-of-band identity). A single shared bind validates the input
// once. Because `.pgen` is involved the whole target must be LOCAL. On ANY failure
// ALL siblings (prefix.pvar / prefix.pgen / prefix.pgen.tmp / prefix.psam) are
// removed — DuckDB only auto-removes its single COPY-target path.
// ===========================================================================

struct PfileWriteBindData : public FunctionData {
	unique_ptr<PvarWriteBindData> pvar_bind;
	unique_ptr<PgenWriteBindData> pgen_bind;
	vector<string> iids;  //! out-of-band sample identity (struct field names or `samples`)
	string pvar_override; //! `pvar` path option (empty -> prefix + ".pvar")
	string pgen_override; //! `pgen` path option (empty -> prefix + ".pgen")
	string psam_override; //! `psam` path option (empty -> prefix + ".psam")

	unique_ptr<FunctionData> Copy() const override {
		auto r = make_uniq<PfileWriteBindData>();
		r->pvar_bind = unique_ptr_cast<FunctionData, PvarWriteBindData>(pvar_bind->Copy());
		r->pgen_bind = unique_ptr_cast<FunctionData, PgenWriteBindData>(pgen_bind->Copy());
		r->iids = iids;
		r->pvar_override = pvar_override;
		r->pgen_override = pgen_override;
		r->psam_override = psam_override;
		return std::move(r);
	}
	bool Equals(const FunctionData &other_p) const override {
		auto &o = other_p.Cast<PfileWriteBindData>();
		return pvar_bind->Equals(*o.pvar_bind) && pgen_bind->Equals(*o.pgen_bind) && iids == o.iids &&
		       pvar_override == o.pvar_override && pgen_override == o.pgen_override && psam_override == o.psam_override;
	}
};

struct PfileWriteGlobalState : public GlobalFunctionData {
	// ATOMIC STAGING. All three writers target `.staging` paths in the SAME directory
	// as the final files; InitGlobal therefore NEVER touches an existing fileset. Only
	// a SUCCESSFUL finalize renames staging -> final (overwriting intra-filesystem via
	// rename(2)). So an upstream-query error that bypasses the sink/finalize cleanup
	// cannot clobber a prior fileset — it only leaves staging files, removed by the
	// destructor. (`pvar_g.path` / `pgen_g.path` / `psam_path` all hold STAGING paths;
	// pgenlib owns `<pgen.staging>.tmp` underneath, hence the distinct `.staging`
	// suffix — never `.tmp`, which is pgenlib's.)
	PvarWriteGlobalState pvar_g; //! pvar_g.path = <pvar_final>.staging
	PgenWriteGlobalState pgen_g; //! pgen_g.path = <pgen_final>.staging
	string psam_path;            //! staging .psam path
	string pvar_final;
	string pgen_final;
	string psam_final;
	FileSystem *fs = nullptr; //! captured in InitGlobal so the destructor can clean staging
	bool cleaned = false;
	bool finished = false; //! staging fully promoted to final (don't clean in the destructor)

	// Remove the STAGING siblings we created, in a safe order: close the pvar handle,
	// then close pgenlib + remove .pgen.staging/.pgen.staging.tmp, then .psam.staging,
	// then .pvar.staging. NEVER touches the final files. Idempotent.
	void CleanupPartial(FileSystem &sysfs) {
		if (cleaned) {
			return;
		}
		cleaned = true;
		if (pvar_g.handle) {
			pvar_g.handle.reset();
		}
		pgen_g.CleanupPartial(sysfs); // closes STPgenWriter + removes .pgen.staging(+.tmp)
		for (const string *p : {&psam_path, &pvar_g.path}) {
			try {
				if (!p->empty() && sysfs.FileExists(*p)) {
					sysfs.RemoveFile(*p);
				}
			} catch (...) {
			}
		}
	}

	// On a fully-written fileset, rename each staging file onto its final path. The
	// `.pgen` is renamed LAST on purpose: its absence makes `read_pfile` fail cleanly
	// (no .pgen), so until it lands the fileset reads as "not there" rather than
	// half-new. The three renames are NOT atomic as a GROUP: a failure between them (these
	// MoveFile calls are not individually guarded — a throw unwinds into PfileCopyFinalize's
	// catch, whose CleanupPartial removes only staging, never the finals) can leave a NEW
	// .pvar/.psam beside an OLD .pgen. When the new fileset has a DIFFERENT variant count
	// than the old .pgen, `read_pfile` catches that as a variant-count mismatch — a clean
	// error, never a silent mis-serve. The ONE edge it does not cover: a `.pgen`-rename
	// failure on a SAME-variant-count overwrite (new .pvar/.psam already promoted, old .pgen
	// left in place) would serve the new metadata over the old genotypes undetected (the
	// counts match, so the reader's guard stays silent). That requires an intra-directory
	// rename(2) to fail on a just-fsync'd file — effectively impossible on a local fs — and
	// is still strictly better than the pre-staging truncate-in-place, which corrupted the
	// real targets on ANY mid-write failure. We do not claim the group-atomicity we lack.
	void PromoteStaging(FileSystem &sysfs) {
		sysfs.MoveFile(pvar_g.path, pvar_final);
		sysfs.MoveFile(psam_path, psam_final);
		sysfs.MoveFile(pgen_g.path, pgen_final); // LAST
		finished = true;
	}

	// Catches an UPSTREAM-query error (e.g. a cast failure) that errors the pipeline
	// WITHOUT running the sink/finalize cleanup: remove the staging files so they are
	// never orphaned. Only staging — the final fileset is only ever written by a
	// successful PromoteStaging, after which finished == true and this is a no-op.
	~PfileWriteGlobalState() override {
		if (!finished && fs) {
			CleanupPartial(*fs);
		}
	}
};

struct PfileWriteLocalState : public LocalFunctionData {};

// Strip a trailing ".pfile" from the COPY target so both `TO 'out'` and
// `TO 'out.pfile'` yield the companion prefix "out".
string PfilePrefix(const string &file_path) {
	const string ext = ".pfile";
	if (file_path.size() > ext.size() && StringUtil::Lower(file_path.substr(file_path.size() - ext.size())) == ext) {
		return file_path.substr(0, file_path.size() - ext.size());
	}
	return file_path;
}

unique_ptr<FunctionData> PfileCopyBind(ClientContext &context, CopyFunctionBindInput &input, const vector<string> &names,
                                       const vector<LogicalType> &sql_types) {
	auto bind = make_uniq<PfileWriteBindData>();

	// --- path overrides (reader symmetry), parsed first so the local-only check
	//     below sees the EFFECTIVE target of each of the three files. ---
	for (auto &opt : input.info.options) {
		auto key = StringUtil::Lower(opt.first);
		if (key == "pvar" && !opt.second.empty() && !opt.second[0].IsNull()) {
			bind->pvar_override = opt.second[0].ToString();
		} else if (key == "psam" && !opt.second.empty() && !opt.second[0].IsNull()) {
			bind->psam_override = opt.second[0].ToString();
		} else if (key == "pgen" && !opt.second.empty() && !opt.second[0].IsNull()) {
			bind->pgen_override = opt.second[0].ToString();
		}
	}

	// --- local-only (the whole fileset): .pgen is written through a raw FILE* that
	//     bypasses DuckDB's VFS, so EVERY one of the three files must be local. Check
	//     the EFFECTIVE target of each: the prefix-derived path when no override, the
	//     override when present. This rejects a remote PREFIX regardless of any
	//     per-file override (a remote prefix + local `pgen` override must NOT bind and
	//     then write .pvar/.psam to the remote store). Done BEFORE BuildPgenBind so the
	//     message is named for the pfile bundle, not for `write pgen`. ---
	string prefix = PfilePrefix(input.info.file_path);
	string eff_pvar = bind->pvar_override.empty() ? prefix + ".pvar" : bind->pvar_override;
	string eff_pgen = bind->pgen_override.empty() ? prefix + ".pgen" : bind->pgen_override;
	string eff_psam = bind->psam_override.empty() ? prefix + ".psam" : bind->psam_override;
	for (const string *t : {&eff_pvar, &eff_pgen, &eff_psam}) {
		if (!t->empty() && !PathIsLocalWritable(*t)) {
			throw BinderException(
			    "write pfile: '%s' is not a local path. The pfile bundle writes a .pgen through a raw file handle that "
			    "bypasses DuckDB's VFS, so the whole fileset (.pvar/.pgen/.psam) must be written to a local path (no "
			    "s3://, http://, or pathmacro:// targets).",
			    *t);
		}
	}

	// Genotypes / pgen roles (also detects STRUCT genotypes). Its own local-check is
	// now redundant — every effective target was just confirmed local above.
	bind->pgen_bind = BuildPgenBind(input, names, sql_types);

	// pvar roles + carried-through extras, EXCLUDING the genotypes column (so it is
	// never stringified into the .pvar as a bogus extra column).
	std::set<idx_t> exclude = {bind->pgen_bind->geno_input_idx};
	bind->pvar_bind = BuildPvarBind(input, names, sql_types, exclude);

	// Sample identity (the contract): struct field names or a `samples` list; also
	// rejects sample_fields metadata beyond IID.
	string iid_col = "IID";
	bool has_iid_map = false;
	CheckSampleFieldsIid(input, iid_col, has_iid_map);
	bind->iids = ResolveIdentity(input, bind->pgen_bind->struct_iids, bind->pgen_bind->geno_is_struct);

	// When the genotypes column FIXES the sample count at bind (a fixed-size ARRAY,
	// or a STRUCT), a `samples` list of a different length is a bind-time error —
	// clearer than deferring to the first-row length check. (For a STRUCT the two
	// are equal by construction; for a plain LIST sample_ct_known is false and the
	// check is deferred to the first row.)
	if (bind->pgen_bind->sample_ct_known &&
	    bind->pgen_bind->sample_ct != static_cast<uint32_t>(bind->iids.size())) {
		throw BinderException("write pfile: the `samples` list has %llu IIDs but the genotypes column has %u "
		                      "per-sample slots — they must match",
		                      static_cast<unsigned long long>(bind->iids.size()), bind->pgen_bind->sample_ct);
	}
	// The IID count is authoritative for the .pgen sample count: this makes the
	// pgen writer validate EVERY variant's genotype length against it — so a
	// `samples` list whose length != the genotype count fails cleanly (at the first
	// row for a LIST input, where sample_ct_known is false).
	bind->pgen_bind->sample_ct = static_cast<uint32_t>(bind->iids.size());
	bind->pgen_bind->sample_ct_known = true;

	(void)context;
	return std::move(bind);
}

unique_ptr<GlobalFunctionData> PfileCopyInitGlobal(ClientContext &context, FunctionData &bind_data,
                                                   const string &file_path) {
	auto &bind = bind_data.Cast<PfileWriteBindData>();
	auto &fs = FileSystem::GetFileSystem(context);
	auto g = make_uniq<PfileWriteGlobalState>();

	// Final destinations, and the `.staging` siblings we actually write first (same
	// directory → intra-filesystem rename at promote; never truncate the real files).
	string prefix = PfilePrefix(file_path);
	g->pvar_final = bind.pvar_override.empty() ? prefix + ".pvar" : bind.pvar_override;
	g->pgen_final = bind.pgen_override.empty() ? prefix + ".pgen" : bind.pgen_override;
	g->psam_final = bind.psam_override.empty() ? prefix + ".psam" : bind.psam_override;
	string pvar_staging = g->pvar_final + ".staging";
	string pgen_staging = g->pgen_final + ".staging";
	string psam_staging = g->psam_final + ".staging";
	g->psam_path = psam_staging;
	g->fs = &fs; // so the destructor can remove staging on an upstream-query error

	// pgen writer init is DEFERRED to the first sink row (same as standalone pgen):
	// carry the staging path + sample_ct + bound in its global state now. pgenlib owns
	// `<pgen_staging>.tmp` under this.
	g->pgen_g.path = pgen_staging;
	g->pgen_g.sample_ct_known = bind.pgen_bind->sample_ct_known;
	g->pgen_g.sample_ct = bind.pgen_bind->sample_ct;
	g->pgen_g.variant_ct_limit = bind.pgen_bind->variant_ct_limit;

	try {
		// pvar: open staging + header. psam: written in full to staging now (out-of-band
		// identity). Neither touches the existing final fileset.
		PvarOpenAndWriteHeader(g->pvar_g, fs, *bind.pvar_bind, pvar_staging);
		WritePsamFromIids(fs, psam_staging, bind.iids);
	} catch (...) {
		g->CleanupPartial(fs);
		throw;
	}
	return std::move(g);
}

unique_ptr<LocalFunctionData> PfileCopyInitLocal(ExecutionContext &, FunctionData &) {
	return make_uniq<PfileWriteLocalState>();
}

void PfileCopySink(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate, LocalFunctionData &,
                   DataChunk &input) {
	auto &bind = bind_data.Cast<PfileWriteBindData>();
	auto &g = gstate.Cast<PfileWriteGlobalState>();
	auto &fs = FileSystem::GetFileSystem(context.client);
	try {
		// pvar first: its (CHROM, POS) monotonicity check fires here, so unsorted
		// input throws before any .pgen append — and the catch removes all siblings.
		PvarWriteChunk(*bind.pvar_bind, g.pvar_g, input);
		PgenWriteChunk(*bind.pgen_bind, g.pgen_g, input);
	} catch (...) {
		g.CleanupPartial(fs);
		throw;
	}
}

void PfileCopyCombine(ExecutionContext &, FunctionData &, GlobalFunctionData &, LocalFunctionData &) {
}

void PfileCopyFinalize(ClientContext &context, FunctionData &, GlobalFunctionData &gstate) {
	auto &g = gstate.Cast<PfileWriteGlobalState>();
	auto &fs = FileSystem::GetFileSystem(context);
	try {
		// pgen finalize throws on zero variants (pgenlib cannot finish a 0-variant
		// .pgen) → cleanup removes the staging .pvar/.psam too (finals untouched).
		PgenFinish(g.pgen_g);
		if (g.pvar_g.handle) {
			g.pvar_g.handle->Sync();
			g.pvar_g.handle.reset();
		}
		// All three staging files are complete — promote them to the final paths
		// (rename, overwriting any prior fileset intra-filesystem). ONLY here does an
		// existing fileset change; an earlier failure never reached this point.
		g.PromoteStaging(fs);
	} catch (...) {
		g.CleanupPartial(fs);
		throw;
	}
}

CopyFunctionExecutionMode PfileCopyExecMode(bool, bool) {
	// Single-threaded, ordered: both the .pgen append and the .pvar monotonicity
	// check require rows in query order.
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

	// FORMAT pgen — binary .pgen writer (P9 Phase B, biallelic hardcalls only).
	CopyFunction pgen_fn("pgen");
	pgen_fn.copy_to_bind = PgenCopyBind;
	pgen_fn.copy_to_initialize_local = PgenCopyInitLocal;
	pgen_fn.copy_to_initialize_global = PgenCopyInitGlobal;
	pgen_fn.copy_to_sink = PgenCopySink;
	pgen_fn.copy_to_combine = PgenCopyCombine;
	pgen_fn.copy_to_finalize = PgenCopyFinalize;
	pgen_fn.execution_mode = PgenCopyExecMode;
	pgen_fn.extension = "pgen";
	loader.RegisterFunction(pgen_fn);

	// FORMAT psam — text .psam writer (P9 Phase C, IID-only).
	CopyFunction psam_fn("psam");
	psam_fn.copy_to_bind = PsamCopyBind;
	psam_fn.copy_to_initialize_local = PsamCopyInitLocal;
	psam_fn.copy_to_initialize_global = PsamCopyInitGlobal;
	psam_fn.copy_to_sink = PsamCopySink;
	psam_fn.copy_to_combine = PsamCopyCombine;
	psam_fn.copy_to_finalize = PsamCopyFinalize;
	psam_fn.execution_mode = PsamCopyExecMode;
	psam_fn.extension = "psam";
	loader.RegisterFunction(psam_fn);

	// FORMAT pfile — the bundle: one COPY fans out .pvar + .pgen + .psam (P9 Phase C).
	CopyFunction pfile_fn("pfile");
	pfile_fn.copy_to_bind = PfileCopyBind;
	pfile_fn.copy_to_initialize_local = PfileCopyInitLocal;
	pfile_fn.copy_to_initialize_global = PfileCopyInitGlobal;
	pfile_fn.copy_to_sink = PfileCopySink;
	pfile_fn.copy_to_combine = PfileCopyCombine;
	pfile_fn.copy_to_finalize = PfileCopyFinalize;
	pfile_fn.execution_mode = PfileCopyExecMode;
	pfile_fn.extension = "pfile";
	loader.RegisterFunction(pfile_fn);
}

} // namespace duckdb
