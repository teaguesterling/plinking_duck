#include "howard_reader.hpp"
#include "howard_catalog.hpp"
#include "plink_common.hpp"
#include "register_helper.hpp"

#include "duckdb/common/file_system.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include "yyjson.hpp"

#include <algorithm>
#include <cstring>
#include <limits>

using namespace duckdb_yyjson; // NOLINT

namespace duckdb {

// ===========================================================================
// read_howard_manifest — the parse + verify + prune surface plinking_duck owns for
// `.howard` READING (the sharded pfile writer owns the WRITE, src/pfile_writer.cpp).
//
//   read_howard_manifest(path VARCHAR
//       [, region := VARCHAR, variant_index := ANY, sample_index := ANY, verify_all := BOOLEAN])
//
// Parses a `.howard.json` with yyjson, resolves every path relative to the catalog's own
// directory (§7.2), folds the `shared` companion into each shard, prunes shards a predicate
// cannot touch (§6/§10 overlap), verifies the survivors' `.pgen` header counts against their
// companions (and their declared `range` against the real `.pvar` extent), and returns one row
// per SURVIVING shard. Phase A: this is the manifest only — NOT the `read_howard` scan or ATTACH.
//
// ERROR TAXONOMY (distinct message prefixes so the Python consumer can classify — spec §8/§10.1):
//   "read_howard_manifest: catalog parse error: ..."                 malformed JSON / schema  (InvalidInputException)
//   "read_howard_manifest: file not found: ..."                      a referenced file missing (IOException)
//   "read_howard_manifest: verification failed (count mismatch): ..." .pgen header vs companion (InvalidInputException)
//   "read_howard_manifest: verification failed (range under-claim): ..." declared range < real  (InvalidInputException)
//   "read_howard_manifest: invalid predicate: ..."                   wrong-axis / malformed arg (InvalidInputException)
// ===========================================================================

namespace {

constexpr int64_t HOWARD_VERSION = 1;

// ---- yyjson helpers -------------------------------------------------------

//! RAII for a parsed yyjson document, so every throw path frees it.
struct YyDocGuard {
	yyjson_doc *doc = nullptr;
	~YyDocGuard() {
		if (doc) {
			yyjson_doc_free(doc);
		}
	}
};

//! Read a JSON integer (yyjson stores non-negative ints as UINT, negatives as SINT; the safe
//! typed accessors each return 0 for the other subtype, so branch explicitly). An unsigned value
//! that will not fit in int64 is REJECTED (returns false) rather than silently wrapping negative
//! — the caller turns that into a `catalog parse error` (M3).
bool JsonGetInt64(yyjson_val *v, int64_t &out) {
	if (!v) {
		return false;
	}
	if (yyjson_is_uint(v)) {
		uint64_t u = yyjson_get_uint(v);
		if (u > static_cast<uint64_t>(std::numeric_limits<int64_t>::max())) {
			return false; // out of int64 range -> parse error upstream
		}
		out = static_cast<int64_t>(u);
		return true;
	}
	if (yyjson_is_sint(v)) {
		out = yyjson_get_sint(v);
		return true;
	}
	return false;
}

// ---- path resolution (§7.2) ----------------------------------------------

//! Collapse `.`/`..`/duplicate separators in a POSIX-style path (os.path.normpath semantics, enough
//! for catalog paths). Preserves a leading `/` (absolute) and a trailing-less form.
string NormPath(const string &path) {
	if (path.empty()) {
		return ".";
	}
	bool absolute = path[0] == '/';
	vector<string> out;
	size_t i = 0;
	while (i < path.size()) {
		size_t j = path.find('/', i);
		if (j == string::npos) {
			j = path.size();
		}
		string seg = path.substr(i, j - i);
		i = j + 1;
		if (seg.empty() || seg == ".") {
			continue;
		}
		if (seg == "..") {
			if (!out.empty() && out.back() != "..") {
				out.pop_back();
			} else if (!absolute) {
				out.push_back(seg); // keep leading .. on a relative path
			}
			continue;
		}
		out.push_back(seg);
	}
	string joined;
	for (size_t k = 0; k < out.size(); k++) {
		if (k != 0) {
			joined += "/";
		}
		joined += out[k];
	}
	if (absolute) {
		return "/" + joined;
	}
	return joined.empty() ? "." : joined;
}

//! Absolute-ize `path` against the process CWD (os.path.abspath), then normalize.
string AbsPath(FileSystem &fs, const string &path) {
	if (fs.IsPathAbsolute(path)) {
		return NormPath(path);
	}
	return NormPath(FileSystem::GetWorkingDirectory() + "/" + path);
}

//! Directory portion of a path. "." when the path has no separator.
string DirName(const string &path) {
	auto slash = path.find_last_of("/\\");
	if (slash == string::npos) {
		return ".";
	}
	if (slash == 0) {
		return path.substr(0, 1); // root "/"
	}
	return path.substr(0, slash);
}

//! Resolve a catalog path relative to the catalog's directory by default (§7.2); an absolute path
//! is honored as-is, and `path_is_relative=false` opts one path out. The result is always ABSOLUTE
//! and NORMALIZED (CWD-independent, `..` collapsed), so a consumer can use it as a chunk source path
//! regardless of its own working directory (I5 / catalog.py `_resolve`).
string ResolvePath(FileSystem &fs, const string &base_dir, const string &path, bool path_is_relative) {
	if (!path_is_relative || fs.IsPathAbsolute(path)) {
		return NormPath(path);
	}
	return NormPath(fs.JoinPath(base_dir, path));
}

// ---- SQL quoting ----------------------------------------------------------

string SqlStr(const string &s) {
	return "'" + StringUtil::Replace(s, "'", "''") + "'";
}

// ---- region string parsing (mirrors ParseRegion strictness, string half only) --------------

struct RegionPredicate {
	string contig;
	int64_t pos_lo = std::numeric_limits<int64_t>::lowest();
	int64_t pos_hi = std::numeric_limits<int64_t>::max();
};

//! Parse `"contig:lo-hi"` or a bare `"contig"` (whole-contig). Reject anything else. A manifest
//! prune must NOT open any `.pvar`, so this resolves only the string — contig + inclusive POS
//! window — never variant indices (unlike plink_common's ParseRegion).
RegionPredicate ParseRegionString(const string &region) {
	RegionPredicate out;
	auto colon = region.find(':');
	if (colon == string::npos) {
		if (region.empty()) {
			throw InvalidInputException("read_howard_manifest: invalid predicate: region is empty "
			                            "(expected 'contig:lo-hi' or a bare 'contig')");
		}
		out.contig = region; // bare contig: whole-contig window
		return out;
	}
	if (colon == 0) {
		throw InvalidInputException("read_howard_manifest: invalid predicate: region '%s' has an empty "
		                            "contig (expected 'contig:lo-hi')",
		                            region);
	}
	out.contig = region.substr(0, colon);
	string range_part = region.substr(colon + 1);
	auto dash = range_part.find('-');
	if (dash == string::npos) {
		throw InvalidInputException("read_howard_manifest: invalid predicate: region '%s' is not "
		                            "'contig:lo-hi'",
		                            region);
	}
	string lo_s = range_part.substr(0, dash);
	string hi_s = range_part.substr(dash + 1);
	try {
		out.pos_lo = std::stoll(lo_s);
		out.pos_hi = std::stoll(hi_s);
	} catch (std::exception &) {
		throw InvalidInputException("read_howard_manifest: invalid predicate: region '%s' has "
		                            "non-integer positions",
		                            region);
	}
	if (out.pos_hi < out.pos_lo) {
		throw InvalidInputException("read_howard_manifest: invalid predicate: region '%s' has inverted "
		                            "POS bounds (lo=%lld > hi=%lld)",
		                            region, (long long)out.pos_lo, (long long)out.pos_hi);
	}
	return out;
}

//! A 2-element inclusive index predicate `[lo, hi]` from an ANY value (a 2-element LIST).
void ParseIndexPredicate(const Value &val, const string &param, int64_t &lo, int64_t &hi) {
	if (val.type().id() != LogicalTypeId::LIST) {
		throw InvalidInputException("read_howard_manifest: invalid predicate: %s must be a 2-element "
		                            "list [lo, hi] (inclusive)",
		                            param);
	}
	auto &children = ListValue::GetChildren(val);
	if (children.size() != 2) {
		throw InvalidInputException("read_howard_manifest: invalid predicate: %s must have exactly 2 "
		                            "elements [lo, hi], got %llu",
		                            param, (unsigned long long)children.size());
	}
	lo = children[0].GetValue<int64_t>();
	hi = children[1].GetValue<int64_t>();
	if (hi < lo) {
		throw InvalidInputException("read_howard_manifest: invalid predicate: %s has inverted bounds "
		                            "(lo=%lld > hi=%lld)",
		                            param, (long long)lo, (long long)hi);
	}
}

// ---- companion value resolution ------------------------------------------

//! A companion value is a string path or an object {path, path_is_relative?, columns?} (§7.2).
//! `columns` is accepted and ignored this phase. Returns the resolved path.
string ResolveCompanion(FileSystem &fs, const string &base_dir, yyjson_val *v, const char *where) {
	if (yyjson_is_str(v)) {
		return ResolvePath(fs, base_dir, yyjson_get_str(v), true);
	}
	if (yyjson_is_obj(v)) {
		yyjson_val *p = yyjson_obj_get(v, "path");
		if (!p || !yyjson_is_str(p)) {
			throw InvalidInputException("read_howard_manifest: catalog parse error: %s companion object "
			                            "needs a string 'path'",
			                            where);
		}
		bool rel = true;
		yyjson_val *r = yyjson_obj_get(v, "path_is_relative");
		if (r && yyjson_is_bool(r)) {
			rel = yyjson_get_bool(r);
		}
		return ResolvePath(fs, base_dir, yyjson_get_str(p), rel);
	}
	throw InvalidInputException("read_howard_manifest: catalog parse error: %s companion must be a "
	                            "string path or an object",
	                            where);
}

// ---- catalog parse (yyjson) -----------------------------------------------

//! Parse + validate a multi-file `.howard.json`. Single-fileset (`pfile`) catalogs are out of
//! scope for the manifest reader this phase and raise a clear parse error.
HowardCatalog ParseCatalog(FileSystem &fs, const string &path, const string &content) {
	YyDocGuard guard;
	guard.doc = yyjson_read(content.c_str(), content.size(), 0);
	if (!guard.doc) {
		throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' is not valid JSON", path);
	}
	yyjson_val *root = yyjson_doc_get_root(guard.doc);
	if (!yyjson_is_obj(root)) {
		throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' top level must be a "
		                            "JSON object",
		                            path);
	}

	HowardCatalog cat;
	// Paths resolve relative to the catalog's own directory, ABSOLUTE-ized against CWD (§7.2 /
	// catalog.py: base_dir = dirname(abspath(path))) so every resolved shard path is CWD-independent.
	cat.base_dir = DirName(AbsPath(fs, path));

	yyjson_val *ver = yyjson_obj_get(root, "howard");
	int64_t version = 0;
	if (!ver || !JsonGetInt64(ver, version)) {
		throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' missing required "
		                            "integer version key 'howard'",
		                            path);
	}
	if (version != HOWARD_VERSION) {
		throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' unsupported 'howard' "
		                            "version %lld (this reader supports %lld)",
		                            path, (long long)version, (long long)HOWARD_VERSION);
	}
	cat.version = version;

	yyjson_val *ds = yyjson_obj_get(root, "dataset");
	if (!ds || !yyjson_is_str(ds)) {
		throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' missing required "
		                            "string field 'dataset'",
		                            path);
	}
	cat.dataset = yyjson_get_str(ds);

	yyjson_val *pfile = yyjson_obj_get(root, "pfile");
	bool has_pfile = pfile != nullptr;
	bool has_multi = yyjson_obj_get(root, "files") != nullptr || yyjson_obj_get(root, "concat_dim") != nullptr ||
	                 yyjson_obj_get(root, "shared") != nullptr;
	if (has_pfile && has_multi) {
		throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' has both a single-fileset "
		                            "'pfile' and multi-file ('files'/'concat_dim'/'shared') keys; a catalog is "
		                            "one or the other",
		                            path);
	}
	if (!has_pfile && !has_multi) {
		throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' must declare either a "
		                            "'pfile' (single fileset) or 'concat_dim'+'shared'+'files' (multi-file)",
		                            path);
	}

	// ---- single-fileset ('pfile') form (catalog.py `_parse_single`) ----
	// One fileset → one row; `concat_dim` is empty (NULL downstream — the signal howard routes on).
	// The partition/chrom/pos are not declared here; bind fills them by probing the single `.pvar`.
	if (has_pfile) {
		cat.is_single = true;
		if (!yyjson_is_obj(pfile)) {
			throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' 'pfile' must be an "
			                            "object with a 'prefix'",
			                            path);
		}
		yyjson_val *pfx = yyjson_obj_get(pfile, "prefix");
		if (!pfx || !yyjson_is_str(pfx)) {
			throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' 'pfile' needs a string "
			                            "'prefix'",
			                            path);
		}
		// Optional `chunks` {dim, chunk_size} — validated like catalog.py `_parse_single` (only variant
		// chunking this phase; a positive chunk_size). The value itself is not surfaced in the row (the
		// coordinator-fixed shape has no chunk_size column), but a malformed one must still be rejected
		// here so it does not reach the loader unchallenged.
		yyjson_val *chunks = yyjson_obj_get(root, "chunks");
		if (chunks) {
			if (!yyjson_is_obj(chunks)) {
				throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' 'chunks' must be "
				                            "an object {dim, chunk_size}",
				                            path);
			}
			yyjson_val *cdim = yyjson_obj_get(chunks, "dim");
			if (cdim && (!yyjson_is_str(cdim) || string(yyjson_get_str(cdim)) != "variants")) {
				throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' chunks.dim must "
				                            "be 'variants' (only variant chunking is supported this phase)",
				                            path);
			}
			yyjson_val *csize = yyjson_obj_get(chunks, "chunk_size");
			if (csize) {
				int64_t cs;
				if (!JsonGetInt64(csize, cs) || cs <= 0) {
					throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' "
					                            "chunks.chunk_size must be a positive integer",
					                            path);
				}
				cat.has_chunk_size = true;
				cat.chunk_size = cs; // surfaced in the row for howard's _register_single
			}
		}

		ShardEntry e;
		e.prefix = ResolvePath(fs, cat.base_dir, yyjson_get_str(pfx), true);
		yyjson_val *own_pvar = yyjson_obj_get(pfile, "pvar");
		yyjson_val *own_psam = yyjson_obj_get(pfile, "psam");
		// An omitted companion stays EMPTY (→ NULL column, matching catalog.py's `None`): the probe
		// then lets read_pfile discover it (incl. the `.pvar.parquet`/`.psam.parquet` sidecar and
		// legacy `.bim`/`.fam`), rather than hardcoding a text path that may not exist.
		e.pvar = own_pvar ? ResolveCompanion(fs, cat.base_dir, own_pvar, "pfile.pvar") : "";
		e.psam = own_psam ? ResolveCompanion(fs, cat.base_dir, own_psam, "pfile.psam") : "";
		// partition/range left at defaults; bind probes the .pvar to fill start=0/stop=count/chrom/pos.
		cat.files.push_back(std::move(e));
		return cat;
	}

	yyjson_val *cd = yyjson_obj_get(root, "concat_dim");
	if (!cd || !yyjson_is_str(cd)) {
		throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' needs a string "
		                            "'concat_dim' ('variants' or 'samples')",
		                            path);
	}
	cat.concat_dim = yyjson_get_str(cd);
	bool is_variants = cat.concat_dim == "variants";
	bool is_samples = cat.concat_dim == "samples";
	if (!is_variants && !is_samples) {
		throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' concat_dim must be "
		                            "'variants' or 'samples', got '%s'",
		                            path, cat.concat_dim);
	}

	// §5 invariant: shared companion must match the concat axis and cannot carry the forbidden one.
	yyjson_val *shared = yyjson_obj_get(root, "shared");
	if (!yyjson_is_obj(shared)) {
		throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' multi-file catalog "
		                            "needs a 'shared' companion object",
		                            path);
	}
	const char *expected_shared = is_variants ? "psam" : "pvar";
	const char *forbidden_shared = is_variants ? "pvar" : "psam";
	if (yyjson_obj_get(shared, forbidden_shared)) {
		throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' concat_dim='%s' "
		                            "requires a shared '%s' but 'shared' carries a '%s' — a contradiction "
		                            "(spec §5)",
		                            path, cat.concat_dim, expected_shared, forbidden_shared);
	}
	yyjson_val *sc = yyjson_obj_get(shared, expected_shared);
	if (!sc) {
		throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' concat_dim='%s' "
		                            "requires 'shared.%s'",
		                            path, cat.concat_dim, expected_shared);
	}
	string shared_resolved = ResolveCompanion(fs, cat.base_dir, sc, expected_shared);
	if (is_variants) {
		cat.shared_psam = shared_resolved;
	} else {
		cat.shared_pvar = shared_resolved;
	}

	yyjson_val *files = yyjson_obj_get(root, "files");
	if (!yyjson_is_arr(files) || yyjson_arr_size(files) == 0) {
		throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' multi-file catalog "
		                            "needs a non-empty 'files' array",
		                            path);
	}

	size_t idx, max;
	yyjson_val *fd;
	yyjson_arr_foreach(files, idx, max, fd) {
		string where = StringUtil::Format("files[%llu]", (unsigned long long)idx);
		if (!yyjson_is_obj(fd)) {
			throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' %s must be an "
			                            "object with a 'prefix'",
			                            path, where);
		}
		yyjson_val *pfx = yyjson_obj_get(fd, "prefix");
		if (!pfx || !yyjson_is_str(pfx)) {
			throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' %s needs a string "
			                            "'prefix'",
			                            path, where);
		}
		// A per-file companion of the SHARED type is a contradiction (spec §5/§7).
		if (yyjson_obj_get(fd, expected_shared)) {
			throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' %s declares a "
			                            "per-file '%s', but the catalog shares one (spec §5/§7)",
			                            path, where, expected_shared);
		}

		ShardEntry e;
		e.prefix = ResolvePath(fs, cat.base_dir, yyjson_get_str(pfx), true);

		// Each shard's own non-shared companion; an omitted one stays EMPTY (→ NULL column, matching
		// catalog.py's `None`) so the probe lets read_pfile discover it (incl. a `.parquet` sidecar /
		// legacy `.bim`/`.fam`). The SHARED companion (always declared) is folded into the other axis.
		yyjson_val *own_pvar = yyjson_obj_get(fd, "pvar");
		yyjson_val *own_psam = yyjson_obj_get(fd, "psam");
		if (is_variants) {
			e.pvar = own_pvar ? ResolveCompanion(fs, cat.base_dir, own_pvar, (where + ".pvar").c_str()) : "";
			e.psam = cat.shared_psam; // fold the shared companion in
		} else {
			e.psam = own_psam ? ResolveCompanion(fs, cat.base_dir, own_psam, (where + ".psam").c_str()) : "";
			e.pvar = cat.shared_pvar;
		}

		// partition {start, stop, source_ordinal} — required, half-open [start, stop) (§6).
		yyjson_val *part = yyjson_obj_get(fd, "partition");
		if (!yyjson_is_obj(part)) {
			throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' %s needs a "
			                            "'partition' object {start, stop, source_ordinal}",
			                            path, where);
		}
		if (!JsonGetInt64(yyjson_obj_get(part, "start"), e.partition_start) ||
		    !JsonGetInt64(yyjson_obj_get(part, "stop"), e.partition_stop) ||
		    !JsonGetInt64(yyjson_obj_get(part, "source_ordinal"), e.source_ordinal)) {
			throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' %s partition needs "
			                            "integer start, stop, source_ordinal",
			                            path, where);
		}
		if (e.partition_stop <= e.partition_start) {
			throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' %s partition "
			                            "[%lld, %lld) is empty or inverted (stop must be > start, §6)",
			                            path, where, (long long)e.partition_start, (long long)e.partition_stop);
		}

		// Optional range {CHROM:[...], POS:{min,max}} — the genomic prune stat (§10).
		yyjson_val *range = yyjson_obj_get(fd, "range");
		if (yyjson_is_obj(range)) {
			e.has_range = true;
			yyjson_val *chrom = yyjson_obj_get(range, "CHROM");
			if (yyjson_is_arr(chrom)) {
				size_t ci, cmax;
				yyjson_val *cv;
				yyjson_arr_foreach(chrom, ci, cmax, cv) {
					// Coerce every element to string (catalog.py `str(c)`): a numeric
					// `"CHROM":[1,2,3]` must participate in prune-drop AND under-claim detection, not
					// be silently dropped to an empty list (I4).
					if (yyjson_is_str(cv)) {
						e.chrom.emplace_back(yyjson_get_str(cv));
					} else {
						int64_t ci_val;
						if (JsonGetInt64(cv, ci_val)) {
							e.chrom.emplace_back(std::to_string(ci_val));
						}
					}
				}
			}
			yyjson_val *pos = yyjson_obj_get(range, "POS");
			if (yyjson_is_obj(pos)) {
				int64_t pmin, pmax;
				if (JsonGetInt64(yyjson_obj_get(pos, "min"), pmin) && JsonGetInt64(yyjson_obj_get(pos, "max"), pmax)) {
					e.has_pos = true;
					e.pos_min = pmin;
					e.pos_max = pmax;
				}
			}
		}

		cat.files.push_back(std::move(e));
	}

	// Canonicalize to ascending partition.start (catalog.py `parse_catalog`): the loader's canonical
	// order. Two shards sharing a start cannot both be valid — each must occupy a distinct range on
	// the concat axis — so that is a parse error. (The cross-shard [0, total) tiling/contiguity check
	// is deliberately NOT done here: it stays in howard, which owns rebasing and re-validates the
	// registered set — a design seam. This reader asserts only per-shard width + sort + dup-start.)
	std::sort(cat.files.begin(), cat.files.end(),
	          [](const ShardEntry &a, const ShardEntry &b) { return a.partition_start < b.partition_start; });
	for (size_t i = 1; i < cat.files.size(); i++) {
		if (cat.files[i].partition_start == cat.files[i - 1].partition_start) {
			throw InvalidInputException("read_howard_manifest: catalog parse error: '%s' two files declare the "
			                            "same partition start %lld — each shard must occupy a distinct range on "
			                            "the '%s' axis",
			                            path, (long long)cat.files[i].partition_start, cat.concat_dim);
		}
	}

	return cat;
}

// ---- verification (spec §8/§10.1) -----------------------------------------

//! A shard's non-shared companion is either DECLARED (explicit path) or OMITTED (read_pfile then
//! discovers it at bind). This resolves an omitted companion EXACTLY the way read_pfile does
//! (src/pfile_reader.cpp, "Discover .pvar/.bim" / ".psam/.fam"): the declared path wins
//! unconditionally; otherwise FindCompanionFileWithParquet returns the first existing default in
//! read_pfile's order — the `.pvar.parquet` / `.psam.parquet` sidecar FIRST when
//! `plinking_use_parquet_companions` is on, then the native text `.pvar`/`.bim` (or `.psam`/`.fam`).
//! Returns "" when nothing exists, so a caller can raise the correct (missing-file) error class.
//! Mirroring read_pfile's exact order keeps the manifest from ever disagreeing with what read_pfile
//! would actually open. (NB: read_pfile's extension lists are {".pvar",".bim"} / {".psam",".fam"} —
//! deliberately NO `.zst`, even though IsNativePlinkFormat recognizes `.pvar.zst`; we mirror that.)
enum class CompanionKind { PVAR, PSAM };

string ResolveExistingCompanion(ClientContext &context, FileSystem &fs, const string &prefix, const string &declared,
                                CompanionKind kind) {
	if (!declared.empty()) {
		return declared;
	}
	const vector<string> exts =
	    kind == CompanionKind::PVAR ? vector<string> {".pvar", ".bim"} : vector<string> {".psam", ".fam"};
	// `prefix + ".pgen"` is read_pfile's `eff_prefix + ".pgen"` base (FindCompanionFileWithParquet
	// strips the `.pgen` via ReplaceExtension before appending each candidate extension).
	return FindCompanionFileWithParquet(context, fs, prefix + ".pgen", exts);
}

//! Count verify: force read_pfile's bind-time `.pgen`-header-vs-companion check with a zero-row
//! probe and re-raise as a count-mismatch (same discipline as howard's `_verify_header`). Existence
//! is checked first so a surviving error is unambiguously a count lie, not a missing file.
void VerifyCounts(ClientContext &context, FileSystem &fs, const ShardEntry &e) {
	string pgen = e.prefix + ".pgen";
	if (!fs.FileExists(pgen)) {
		throw IOException("read_howard_manifest: file not found: .pgen '%s'", pgen);
	}
	if (!e.pvar.empty() && !fs.FileExists(e.pvar)) {
		throw IOException("read_howard_manifest: file not found: .pvar '%s'", e.pvar);
	}
	if (!e.psam.empty() && !fs.FileExists(e.psam)) {
		throw IOException("read_howard_manifest: file not found: .psam '%s'", e.psam);
	}
	// An OMITTED companion whose default file is simply ABSENT is a MISSING FILE (howard's reference
	// loader existence-checks the default companion → CatalogError), NOT an unreadable one. Left
	// unguarded, read_pfile's discovery fails with "cannot find .pvar" — a message WITHOUT "count
	// mismatch" — so the probe below misclassifies it as "companion unreadable" (→ AlignmentError in
	// the Python consumer). Resolve it the way read_pfile would and raise file-not-found (→ CatalogError)
	// when nothing exists at all.
	if (e.pvar.empty() && ResolveExistingCompanion(context, fs, e.prefix, "", CompanionKind::PVAR).empty()) {
		throw IOException("read_howard_manifest: file not found: .pvar '%s'", e.prefix + ".pvar");
	}
	if (e.psam.empty() && ResolveExistingCompanion(context, fs, e.prefix, "", CompanionKind::PSAM).empty()) {
		throw IOException("read_howard_manifest: file not found: .psam '%s'", e.prefix + ".psam");
	}
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	string probe = "SELECT 1 FROM read_pfile(" + SqlStr(e.prefix);
	if (!e.pvar.empty()) {
		probe += ", pvar := " + SqlStr(e.pvar);
	}
	if (!e.psam.empty()) {
		probe += ", psam := " + SqlStr(e.psam);
	}
	probe += ") LIMIT 0";
	auto res = conn.Query(probe);
	if (res->HasError()) {
		// Discriminate the two failure classes the probe can surface, so the caller is not sent down
		// the wrong branch: a header-vs-companion COUNT mismatch (the §8 trust check) vs an unreadable
		// / wrong-format companion (a corrupt .pgen header, an unparseable .pvar). read_pfile reports a
		// count mismatch with the literal "count mismatch"; anything else is a companion-readability
		// failure.
		string err = res->GetError();
		if (StringUtil::Contains(err, "count mismatch")) {
			throw InvalidInputException("read_howard_manifest: verification failed (count mismatch): shard "
			                            "'%s': %s",
			                            e.prefix, err);
		}
		throw InvalidInputException("read_howard_manifest: verification failed (companion unreadable): "
		                            "shard '%s': %s",
		                            e.prefix, err);
	}
}

//! Real stats from a shard's own `.pvar` (variants axis): row count, POS extent, and the distinct
//! CHROM set (coerced to string). One query serves BOTH the width check and the range under-claim
//! check (and single-file enrichment).
struct PvarStats {
	int64_t count = 0;
	bool has_pos = false;
	int64_t pos_min = 0;
	int64_t pos_max = 0;
	vector<string> chroms;
};

PvarStats QueryPvarStats(ClientContext &context, const string &pvar, const string &shard_prefix) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	auto res = conn.Query("SELECT count(*), min(POS), max(POS), "
	                      "list_sort(list_distinct(list(CHROM::VARCHAR))) FROM read_pvar(" +
	                      SqlStr(pvar) + ")");
	if (res->HasError()) {
		throw InvalidInputException("read_howard_manifest: verification failed (companion unreadable): shard "
		                            "'%s': cannot read .pvar '%s': %s",
		                            shard_prefix, pvar, res->GetError());
	}
	auto row = res->Fetch();
	PvarStats st;
	if (row && row->size() != 0) {
		st.count = row->GetValue(0, 0).GetValue<int64_t>();
		if (!row->GetValue(1, 0).IsNull()) {
			st.has_pos = true;
			st.pos_min = row->GetValue(1, 0).GetValue<int64_t>();
			st.pos_max = row->GetValue(2, 0).GetValue<int64_t>();
		}
		Value chrom_list = row->GetValue(3, 0);
		if (!chrom_list.IsNull()) {
			for (auto &c : ListValue::GetChildren(chrom_list)) {
				if (!c.IsNull()) {
					st.chroms.push_back(c.GetValue<string>());
				}
			}
		}
	}
	return st;
}

//! Count the rows of a `.psam` (samples axis) for the width check.
int64_t QueryPsamCount(ClientContext &context, const string &psam, const string &shard_prefix) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	auto res = conn.Query("SELECT count(*) FROM read_psam(" + SqlStr(psam) + ")");
	if (res->HasError()) {
		throw InvalidInputException("read_howard_manifest: verification failed (companion unreadable): shard "
		                            "'%s': cannot read .psam '%s': %s",
		                            shard_prefix, psam, res->GetError());
	}
	return res->Fetch()->GetValue(0, 0).GetValue<int64_t>();
}

//! Verify a shard's declared width against its real companion row count, and (variants only) its
//! declared `range` against the real `.pvar` extent; for a single-fileset catalog, instead ENRICH
//! the entry from the probed `.pvar` (partition 0..count, source_ordinal 0, chrom/POS).
//!
//!  * WIDTH (catalog.py `_verify_partition_tiling`): declared `partition_stop - partition_start`
//!    MUST equal the real companion count (pvar for variants, psam for samples) — a count mismatch.
//!    This is why `verified` cannot lie: a wrong declared width is caught here.
//!  * RANGE under-claim (spec §10): declared `range` must COVER the real `.pvar` extent (⊇). An
//!    under-claim raises; an over-claim is harmless. A shard with no `range` is skipped.
void VerifyWidthAndRange(ClientContext &context, ShardEntry &e, const string &concat_dim, bool is_single) {
	auto &fs = FileSystem::GetFileSystem(context);
	if (concat_dim == "samples") {
		// Read the SAME file read_pfile would: a declared path, else the discovered default — which may
		// be a `.psam.parquet` sidecar, not the hardcoded text `.psam`. read_psam reads parquet natively.
		// (VerifyCounts has already guaranteed something exists; the fallback is just belt-and-braces.)
		string psam = ResolveExistingCompanion(context, fs, e.prefix, e.psam, CompanionKind::PSAM);
		if (psam.empty()) {
			psam = e.prefix + ".psam";
		}
		int64_t real = QueryPsamCount(context, psam, e.prefix);
		int64_t width = e.partition_stop - e.partition_start;
		if (width != real) {
			throw InvalidInputException("read_howard_manifest: verification failed (count mismatch): shard "
			                            "'%s': partition declares %lld samples ([%lld, %lld)) but the .psam "
			                            "has %lld",
			                            e.prefix, (long long)width, (long long)e.partition_start,
			                            (long long)e.partition_stop, (long long)real);
		}
		return;
	}

	// variants axis (or single-fileset, counted on the variants axis).
	// Read the SAME file read_pfile would: a declared path, else the discovered default — which may be
	// a `.pvar.parquet` sidecar (read_pvar reads a non-native source natively), NOT the hardcoded text
	// `.pvar`. So a shard whose only variant companion is `<prefix>.pvar.parquet` passes the extent check
	// instead of failing to open a text `.pvar` that was never written. (VerifyCounts guards all-absent.)
	string pvar = ResolveExistingCompanion(context, fs, e.prefix, e.pvar, CompanionKind::PVAR);
	if (pvar.empty()) {
		pvar = e.prefix + ".pvar";
	}
	PvarStats st = QueryPvarStats(context, pvar, e.prefix);

	if (is_single) {
		// One fileset: the partition/range are not declared — derive them from the real .pvar.
		e.partition_start = 0;
		e.partition_stop = st.count;
		e.source_ordinal = 0;
		e.has_range = true;
		e.chrom = st.chroms;
		e.has_pos = st.has_pos;
		e.pos_min = st.pos_min;
		e.pos_max = st.pos_max;
		return;
	}

	int64_t width = e.partition_stop - e.partition_start;
	if (width != st.count) {
		throw InvalidInputException("read_howard_manifest: verification failed (count mismatch): shard '%s': "
		                            "partition declares %lld variants ([%lld, %lld)) but the .pvar has %lld",
		                            e.prefix, (long long)width, (long long)e.partition_start,
		                            (long long)e.partition_stop, (long long)st.count);
	}

	if (!e.has_range) {
		return;
	}
	if (!e.chrom.empty()) {
		for (auto &actual : st.chroms) {
			bool declared = false;
			for (auto &c : e.chrom) {
				if (c == actual) {
					declared = true;
					break;
				}
			}
			if (!declared) {
				throw InvalidInputException(
				    "read_howard_manifest: verification failed (range under-claim): shard '%s': declared "
				    "range.CHROM does not cover contig '%s' the .pvar actually holds — a stale/under-claiming "
				    "range would wrongly drop this shard (spec §10)",
				    e.prefix, actual);
			}
		}
	}
	if (e.has_pos && st.has_pos && (st.pos_min < e.pos_min || st.pos_max > e.pos_max)) {
		throw InvalidInputException(
		    "read_howard_manifest: verification failed (range under-claim): shard '%s': declared range.POS "
		    "[%lld, %lld] does not cover the real POS extent [%lld, %lld] — an under-claiming range would "
		    "wrongly drop this shard (spec §10)",
		    e.prefix, (long long)e.pos_min, (long long)e.pos_max, (long long)st.pos_min, (long long)st.pos_max);
	}
}

} // namespace

// ---------------------------------------------------------------------------
// Bind / scan
// ---------------------------------------------------------------------------

struct HowardManifestBindData : public TableFunctionData {
	vector<ShardEntry> survivors;
	string dataset;
	string concat_dim;
	bool is_single = false;      //! single-fileset catalog → emit NULL concat_dim (the signal howard routes on)
	bool has_chunk_size = false; //! surface chunks.chunk_size (single-fileset only; multi-shard → NULL)
	int64_t chunk_size = 0;
};

struct HowardManifestGlobalState : public GlobalTableFunctionState {
	idx_t row_idx = 0;
};

static unique_ptr<FunctionData> HowardManifestBind(ClientContext &context, TableFunctionBindInput &input,
                                                   vector<LogicalType> &return_types, vector<string> &names) {
	auto result = make_uniq<HowardManifestBindData>();
	auto &fs = FileSystem::GetFileSystem(context);

	string path = input.inputs[0].GetValue<string>();
	if (!fs.FileExists(path)) {
		throw IOException("read_howard_manifest: file not found: catalog '%s'", path);
	}
	string content;
	{
		auto handle = fs.OpenFile(path, FileFlags::FILE_FLAGS_READ);
		idx_t size = handle->GetFileSize();
		content.resize(size);
		if (size > 0) {
			handle->Read(reinterpret_cast<void *>(&content[0]), size);
		}
	}

	HowardCatalog cat = ParseCatalog(fs, path, content);
	result->dataset = cat.dataset;
	result->concat_dim = cat.concat_dim;
	result->is_single = cat.is_single;
	result->has_chunk_size = cat.has_chunk_size;
	result->chunk_size = cat.chunk_size;

	// ---- single-fileset: one fileset, always the whole candidate (catalog.py prune_files returns
	// [pfile] for ANY predicate — no pruning, no wrong-axis check). Verify + enrich the one row. ----
	if (cat.is_single) {
		auto &e = cat.files[0];
		VerifyCounts(context, fs, e);
		VerifyWidthAndRange(context, e, cat.concat_dim, /*is_single=*/true);
		result->survivors.push_back(e);
	} else {
		// ---- predicates ----
		bool has_region = false, has_variant_index = false, has_sample_index = false, verify_all = false;
		RegionPredicate region;
		int64_t idx_lo = 0, idx_hi = 0;
		bool use_index_prune = false;

		auto np_region = input.named_parameters.find("region");
		if (np_region != input.named_parameters.end() && !np_region->second.IsNull()) {
			has_region = true;
			region = ParseRegionString(np_region->second.GetValue<string>());
		}
		auto np_vi = input.named_parameters.find("variant_index");
		if (np_vi != input.named_parameters.end() && !np_vi->second.IsNull()) {
			has_variant_index = true;
		}
		auto np_si = input.named_parameters.find("sample_index");
		if (np_si != input.named_parameters.end() && !np_si->second.IsNull()) {
			has_sample_index = true;
		}
		auto np_va = input.named_parameters.find("verify_all");
		if (np_va != input.named_parameters.end() && !np_va->second.IsNull()) {
			verify_all = np_va->second.GetValue<bool>();
		}

		// Wrong-axis coherence (a silently-ignored off-axis predicate would quietly return the whole
		// cohort — the unsafe failure mode, spec §10.1). Each predicate is valid on exactly one axis.
		bool is_variants = cat.concat_dim == "variants";
		if (has_variant_index && !is_variants) {
			throw InvalidInputException("read_howard_manifest: invalid predicate: variant_index requires "
			                            "concat_dim='variants', but this catalog is sharded on '%s'",
			                            cat.concat_dim);
		}
		if (has_region && !is_variants) {
			throw InvalidInputException("read_howard_manifest: invalid predicate: region (genomic) requires "
			                            "concat_dim='variants', but this catalog is sharded on '%s' (variants "
			                            "are identical across sample-shards, so a genomic predicate prunes "
			                            "nothing, spec §10)",
			                            cat.concat_dim);
		}
		if (has_sample_index && is_variants) {
			throw InvalidInputException("read_howard_manifest: invalid predicate: sample_index requires "
			                            "concat_dim='samples', but this catalog is sharded on '%s'",
			                            cat.concat_dim);
		}
		if (has_variant_index) {
			ParseIndexPredicate(np_vi->second, "variant_index", idx_lo, idx_hi);
			use_index_prune = true;
		}
		if (has_sample_index) {
			ParseIndexPredicate(np_si->second, "sample_index", idx_lo, idx_hi);
			use_index_prune = true;
		}

		// ---- verify_all: verify EVERY shard before pruning (closes the dropped-shard hole, §10.1) ----
		if (verify_all) {
			for (auto &e : cat.files) {
				VerifyCounts(context, fs, e);
				VerifyWidthAndRange(context, e, cat.concat_dim, /*is_single=*/false);
			}
		}

		// ---- prune (pure, in C++, via the shared overlap predicates) ----
		for (auto &e : cat.files) {
			if (use_index_prune && !e.IndexOverlaps(idx_lo, idx_hi)) {
				continue;
			}
			if (has_region && !e.RangeOverlaps(region.contig, region.pos_lo, region.pos_hi)) {
				continue;
			}
			result->survivors.push_back(e);
		}

		// ---- verify survivors (default) ----
		if (!verify_all) {
			for (auto &e : result->survivors) {
				VerifyCounts(context, fs, e);
				VerifyWidthAndRange(context, e, cat.concat_dim, /*is_single=*/false);
			}
		}
	}

	names = {"dataset",        "concat_dim", "prefix",  "pvar",    "psam",     "partition_start", "partition_stop",
	         "source_ordinal", "chrom",      "pos_min", "pos_max", "verified", "chunk_size"};
	return_types = {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
	                LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::BIGINT,
	                LogicalType::BIGINT,  LogicalType::BIGINT,  LogicalType::LIST(LogicalType::VARCHAR),
	                LogicalType::BIGINT,  LogicalType::BIGINT,  LogicalType::BOOLEAN,
	                LogicalType::BIGINT};
	return std::move(result);
}

static unique_ptr<GlobalTableFunctionState> HowardManifestInitGlobal(ClientContext &context,
                                                                     TableFunctionInitInput &input) {
	return make_uniq<HowardManifestGlobalState>();
}

static void HowardManifestScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind = data_p.bind_data->Cast<HowardManifestBindData>();
	auto &gstate = data_p.global_state->Cast<HowardManifestGlobalState>();

	idx_t remaining = bind.survivors.size() - gstate.row_idx;
	idx_t count = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);
	for (idx_t i = 0; i < count; i++) {
		auto &e = bind.survivors[gstate.row_idx + i];
		output.data[0].SetValue(i, Value(bind.dataset));
		// A single-fileset catalog emits NULL concat_dim — the signal howard routes a whole-fileset
		// load on, vs a sharded cohort.
		output.data[1].SetValue(i, bind.is_single ? Value(LogicalType::VARCHAR) : Value(bind.concat_dim));
		output.data[2].SetValue(i, Value(e.prefix));
		output.data[3].SetValue(i, e.pvar.empty() ? Value(LogicalType::VARCHAR) : Value(e.pvar));
		output.data[4].SetValue(i, e.psam.empty() ? Value(LogicalType::VARCHAR) : Value(e.psam));
		output.data[5].SetValue(i, Value::BIGINT(e.partition_start));
		output.data[6].SetValue(i, Value::BIGINT(e.partition_stop));
		output.data[7].SetValue(i, Value::BIGINT(e.source_ordinal));
		if (e.has_range) {
			vector<Value> contigs;
			for (auto &c : e.chrom) {
				contigs.push_back(Value(c));
			}
			output.data[8].SetValue(i, Value::LIST(LogicalType::VARCHAR, std::move(contigs)));
		} else {
			output.data[8].SetValue(i, Value(LogicalType::LIST(LogicalType::VARCHAR)));
		}
		output.data[9].SetValue(i, e.has_pos ? Value::BIGINT(e.pos_min) : Value(LogicalType::BIGINT));
		output.data[10].SetValue(i, e.has_pos ? Value::BIGINT(e.pos_max) : Value(LogicalType::BIGINT));
		// Every returned row has been verified during this call (survivors, or all shards when
		// verify_all); a verification failure throws, so this is currently always true. The column
		// exists for a future degrade-gracefully mode that returns unverified shards.
		output.data[11].SetValue(i, Value::BOOLEAN(true));
		// chunk_size: a single-fileset catalog's declared chunks.chunk_size (consumed by howard's
		// _register_single); NULL for multi-shard rows (the shards ARE the chunks — attach_chunks,
		// registered with chunks=False — so there is no chunk_size).
		output.data[12].SetValue(i, bind.has_chunk_size ? Value::BIGINT(bind.chunk_size) : Value(LogicalType::BIGINT));
	}
	gstate.row_idx += count;
	output.SetCardinality(count);
}

void RegisterHowardReader(ExtensionLoader &loader) {
	TableFunction fn("read_howard_manifest", {LogicalType::VARCHAR}, HowardManifestScan, HowardManifestBind,
	                 HowardManifestInitGlobal);
	fn.named_parameters["region"] = LogicalType::VARCHAR;
	fn.named_parameters["variant_index"] = LogicalType::ANY;
	fn.named_parameters["sample_index"] = LogicalType::ANY;
	fn.named_parameters["verify_all"] = LogicalType::BOOLEAN;
	RegisterTableWithDesc(
	    loader, fn, {"path", "region", "variant_index", "sample_index", "verify_all"},
	    "Parse, verify, and prune a .howard.json shard catalog, returning one row per surviving shard. "
	    "The parse+verify+prune surface plinking_duck owns for .howard reading (write side = the sharded "
	    "pfile writer).",
	    {"SELECT * FROM read_howard_manifest('cohort.howard.json')",
	     "SELECT prefix FROM read_howard_manifest('cohort.howard.json', region := '2:1-5000000')"});
}

} // namespace duckdb
