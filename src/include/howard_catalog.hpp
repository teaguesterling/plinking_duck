#pragma once

#include "duckdb/common/common.hpp" // duckdb::string / duckdb::vector aliases

namespace duckdb {

// ---------------------------------------------------------------------------
// Shared `.howard` catalog model (spec:
//   ~/Projects/howard/docs/superpowers/specs/2026-10-09-howard-catalog-design.md §6/§7/§10).
//
// ONE struct for both sides of the `.howard` boundary plinking_duck owns:
//   * the sharded pfile WRITER (`src/pfile_writer.cpp`) fills these as it finalizes each
//     shard, then serializes them to `out.howard.json` (paths RELATIVE to the catalog dir,
//     `psam` left empty — the shared companion lives at the catalog root);
//   * the `read_howard_manifest` READER (`src/howard_reader.cpp`) parses them back, RESOLVES
//     every path against the catalog's own directory, folds the shared companion into each
//     shard, and prunes with the overlap predicates below.
//
// Boundary law (spec §6, stated explicitly so nothing is left to convention):
//   * `partition [start, stop)` is HALF-OPEN   — start included, stop excluded (concat-axis ordinal).
//   * `range.POS [min, max]`    is INCLUSIVE    — both endpoints are coordinates the shard holds.
//   * `range.CHROM`             is MEMBERSHIP   — a contig is covered or it is not.
// ---------------------------------------------------------------------------

//! One shard's catalog facts → one `files[]` entry in a `.howard.json`.
struct ShardEntry {
	string prefix;               //! fileset path, NO extension (writer: relative; reader: resolved)
	string pvar;                 //! this shard's own `.pvar` (writer: relative; reader: resolved)
	string psam;                 //! this shard's `.psam` (reader folds the shared companion here; writer: empty)
	int64_t partition_start = 0; //! global half-open [start, stop) ordinal on the concat axis (§6)
	int64_t partition_stop = 0;
	int64_t source_ordinal = 0;
	vector<string> chrom;   //! range.CHROM — the contigs this shard covers (membership)
	bool has_range = false; //! a `range` object was present (reader) / is emitted (writer)
	bool has_pos = false;   //! range.POS present (inclusive [pos_min, pos_max])
	int64_t pos_min = 0;
	int64_t pos_max = 0;

	//! INDEX prune (spec §6/§10.1). Does this shard's HALF-OPEN partition `[start, stop)` overlap
	//! the INCLUSIVE index predicate `[lo, hi]`? The predicate is inclusive to match SQL
	//! `WHERE variant_index BETWEEN lo AND hi`. A shard holds ordinals `start .. stop-1`, so it
	//! overlaps iff `start <= hi AND stop > lo`. (Boundary: `[1000,1000]` touches only the shard
	//! whose partition contains ordinal 1000, i.e. `[1000, 2000)`, not `[0, 1000)`.)
	bool IndexOverlaps(int64_t lo, int64_t hi) const {
		return partition_start <= hi && partition_stop > lo;
	}

	//! GENOMIC prune (spec §6/§10). CONSERVATIVE: returns true (KEEP) unless the declared `range`
	//! AFFIRMATIVELY PROVES the predicate cannot touch this shard. A shard with no `range` is always
	//! a candidate (never wrongly dropped). Drop only when: the contig is absent from a declared
	//! non-empty CHROM list, OR a declared inclusive POS extent is entirely disjoint from the
	//! predicate window `[pos_lo, pos_hi]` (also inclusive). Contig names are matched VERBATIM
	//! (no `chr`-prefix / numeric normalization — the caller spells the contig as the `.pvar` does).
	bool RangeOverlaps(const string &contig, int64_t pos_lo, int64_t pos_hi) const {
		if (!has_range) {
			return true; // nothing declared -> conservative keep
		}
		if (!chrom.empty()) {
			bool found = false;
			for (auto &c : chrom) {
				if (c == contig) {
					found = true;
					break;
				}
			}
			if (!found) {
				return false; // provably wrong contig -> drop
			}
		}
		if (has_pos && (pos_hi < pos_min || pos_lo > pos_max)) {
			return false; // predicate window disjoint from POS extent -> drop
		}
		return true;
	}
};

//! A parsed `.howard` catalog. `base_dir` is the (absolute) directory every path resolved against.
//! Both forms are supported: a single-fileset (`pfile`) catalog (`is_single`, empty `concat_dim`,
//! one entry) and the multi-file / sharded `files` form.
struct HowardCatalog {
	int64_t version = 0;
	string dataset;
	string base_dir;
	bool is_single = false; //! single-fileset ('pfile') catalog — `concat_dim` is empty / NULL downstream
	string concat_dim;      //! "variants" | "samples" (empty for a single-fileset catalog)
	string shared_pvar;     //! resolved (set when concat_dim == "samples")
	string shared_psam;     //! resolved (set when concat_dim == "variants")
	vector<ShardEntry> files;
};

} // namespace duckdb
