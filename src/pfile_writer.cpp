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

// pgenlib allocates ~a few bytes per variant in `variant_ct_limit` for the
// on-the-fly index. 2^24 ≈ 16.7M variants (exome / chip / single-chromosome
// scale); larger datasets pass `variant_count N`. The sink guard makes an
// exceeded bound a clean error, so this is safe to exceed, not a silent cap.
constexpr uint32_t kDefaultPgenVariantCtLimit = 1u << 24;

bool PathIsLocalWritable(const string &path) {
	// raw FILE* target: reject remote (s3/http/gcs/…) and any scheme:// URL
	// (e.g. pathmacro://). A bare local path (absolute/relative) is writable.
	if (FileSystem::IsRemoteFile(path)) {
		return false;
	}
	if (path.find("://") != string::npos) {
		return false;
	}
	return true;
}

struct PgenWriteBindData : public FunctionData {
	idx_t geno_input_idx = 0;    //! input chunk column index of the genotypes LIST/ARRAY
	bool sample_ct_known = false; //! true iff genotypes is a fixed-size ARRAY (N in the type)
	uint32_t sample_ct = 0;      //! sample count when known at bind
	bool has_alt = false;        //! an ALT-role column is present (multiallelic guard)
	idx_t alt_input_idx = 0;
	uint32_t variant_ct_limit = kDefaultPgenVariantCtLimit;
	string path_override;        //! `pgen` option; empty -> use the COPY target

	unique_ptr<FunctionData> Copy() const override {
		auto r = make_uniq<PgenWriteBindData>();
		r->geno_input_idx = geno_input_idx;
		r->sample_ct_known = sample_ct_known;
		r->sample_ct = sample_ct;
		r->has_alt = has_alt;
		r->alt_input_idx = alt_input_idx;
		r->variant_ct_limit = variant_ct_limit;
		r->path_override = path_override;
		return std::move(r);
	}
	bool Equals(const FunctionData &other_p) const override {
		auto &o = other_p.Cast<PgenWriteBindData>();
		return geno_input_idx == o.geno_input_idx && sample_ct_known == o.sample_ct_known &&
		       sample_ct == o.sample_ct && has_alt == o.has_alt && alt_input_idx == o.alt_input_idx &&
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

unique_ptr<FunctionData> PgenCopyBind(ClientContext &, CopyFunctionBindInput &input, const vector<string> &names,
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
				if (static_cast<uint64_t>(vc) > std::numeric_limits<uint32_t>::max()) {
					throw BinderException("write pgen: variant_count %lld exceeds the .pgen maximum",
					                      static_cast<long long>(vc));
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
	const LogicalType &gtype = sql_types[gidx];
	LogicalTypeId tid = gtype.id();
	if (tid != LogicalTypeId::ARRAY && tid != LogicalTypeId::LIST) {
		throw BinderException("write pgen: genotypes column '%s' has type %s, but a biallelic-hardcall .pgen needs a "
		                      "LIST or ARRAY of per-sample calls (e.g. TINYINT[N])",
		                      names[gidx], gtype.ToString());
	}
	const LogicalType &child = (tid == LogicalTypeId::ARRAY) ? ArrayType::GetChildType(gtype)
	                                                         : ListType::GetChildType(gtype);
	LogicalTypeId cid = child.id();
	if (cid == LogicalTypeId::ARRAY || cid == LogicalTypeId::LIST) {
		throw BinderException("write pgen: genotypes column '%s' looks PHASED (element type %s) — phased .pgen writing "
		                      "is not supported in this increment (biallelic hardcalls only)",
		                      names[gidx], child.ToString());
	}
	if (cid == LogicalTypeId::DOUBLE || cid == LogicalTypeId::FLOAT) {
		throw BinderException("write pgen: genotypes column '%s' looks like DOSAGES (element type %s) — dosage .pgen "
		                      "writing is not supported in this increment (biallelic hardcalls only)",
		                      names[gidx], child.ToString());
	}
	if (!IsIntegerType(child)) {
		throw BinderException("write pgen: genotypes column '%s' has element type %s, but biallelic hardcalls must be "
		                      "an integer (0/1/2/-9)",
		                      names[gidx], child.ToString());
	}
	if (tid == LogicalTypeId::ARRAY) {
		bind->sample_ct = static_cast<uint32_t>(ArrayType::GetSize(gtype));
		bind->sample_ct_known = true;
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

unique_ptr<GlobalFunctionData> PgenCopyInitGlobal(ClientContext &, FunctionData &bind_data, const string &file_path) {
	auto &bind = bind_data.Cast<PgenWriteBindData>();
	auto g = make_uniq<PgenWriteGlobalState>();
	g->path = bind.path_override.empty() ? file_path : bind.path_override;
	g->sample_ct_known = bind.sample_ct_known;
	g->sample_ct = bind.sample_ct;
	g->variant_ct_limit = bind.variant_ct_limit;

	// If sample_ct is known from a fixed-size ARRAY type, init the writer now so
	// even a zero-row query yields a valid empty .pgen. For a LIST input the
	// sample count is only known from the first row, so defer to the sink.
	if (bind.sample_ct_known) {
		PgenEnsureWriterInited(*g, bind.sample_ct);
	}
	return std::move(g);
}

unique_ptr<LocalFunctionData> PgenCopyInitLocal(ExecutionContext &, FunctionData &) {
	return make_uniq<PgenWriteLocalState>();
}

void PgenCopySink(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate, LocalFunctionData &,
                  DataChunk &input) {
	auto &bind = bind_data.Cast<PgenWriteBindData>();
	auto &g = gstate.Cast<PgenWriteGlobalState>();
	auto &fs = FileSystem::GetFileSystem(context.client);

	try {
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
			const vector<Value> &calls = (geno_val.type().id() == LogicalTypeId::ARRAY)
			                                  ? ArrayValue::GetChildren(geno_val)
			                                  : ListValue::GetChildren(geno_val);
			idx_t n = calls.size();

			// --- establish / validate the sample count ---
			if (!g.writer_inited) {
				PgenEnsureWriterInited(g, static_cast<uint32_t>(n));
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

	if (g.finished) {
		return;
	}
	if (!g.writer_inited) {
		// Zero rows AND sample_ct never determined (a LIST-typed zero-row query).
		// No file was opened; a fixed-ARRAY input would have inited in InitGlobal.
		throw IOException("write pgen: no variants were written and the sample count could not be determined from the "
		                  "genotypes type — provide a fixed-size genotypes ARRAY (TINYINT[N]) or at least one row");
	}
	try {
		plink2::PglErr e = plink2::SpgwFinish(&g.spgw);
		if (e != plink2::kPglRetSuccess) {
			throw IOException("write pgen: SpgwFinish failed (pgenlib error %u)", static_cast<uint32_t>(e));
		}
		plink2::PglErr ce = plink2::kPglRetSuccess;
		plink2::CleanupSpgw(&g.spgw, &ce);
		g.finished = true;
	} catch (...) {
		g.CleanupPartial(fs);
		throw;
	}
}

CopyFunctionExecutionMode PgenCopyExecMode(bool /*preserve_insertion_order*/, bool /*supports_batch_index*/) {
	// Single-threaded, ordered: .pgen records are strictly positional and the
	// single STPgenWriter appends in query row order.
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
}

} // namespace duckdb
