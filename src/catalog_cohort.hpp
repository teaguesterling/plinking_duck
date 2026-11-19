// CatalogCohort: a pfile catalog (v1 or v2; variant-sharded, sample-sharded or matrix) resolved for
// the per-variant analysis functions (plink_freq, plink_missing). The variant blocks are read in
// catalog order. Each block is read over all its sample blocks through a ConcatPgenReader, so a
// count adds across sample blocks. The cohort's samples are the sample blocks' .psam rows,
// concatenated in sample-block order.
#pragma once

#include "duckdb.hpp"
#include "plink_common.hpp"
#include "pgen_concat_reader.hpp"
#include "pgen_vfs_opener.hpp"

namespace duckdb {

//! The catalog that an analysis function's inputs name, or "" when they name one .pgen. The rules
//! follow read_pfile:
//!   * catalog := '<x>.pfile.json' names the catalog (with no positional argument);
//!   * a positional '<x>.pfile.json' is that catalog;
//!   * a positional local prefix '<p>' (not ending in '.pgen') is the catalog '<p>.pfile.json' when
//!     that file exists. An explicit '<p>.pgen' is always a file.
//! catalog := false turns discovery off. catalog := true requires a catalog at the prefix.
//! With a catalog, `rejected` parameters (the per-file companions, and any option the catalog read
//! cannot honor) are errors.
string ResolveAnalysisCatalog(ClientContext &context, const TableFunctionBindInput &input,
                              const vector<string> &rejected, const string &func_name);

struct CatalogCohort {
	//! One variant block: its cells' .pgen files in sample-block order, and its .pvar.
	struct Block {
		vector<PgenSampleShard> shards;
		VariantMetadataIndex variants;
		uint32_t variant_ct = 0;
	};
	//! One contiguous run of scanned variants inside one block: block-local [start, end), at scan
	//! position scan_start.
	struct Segment {
		idx_t block = 0;
		uint32_t start = 0;
		uint32_t end = 0;
		uint32_t scan_start = 0;
	};

	string path;
	vector<Block> blocks;
	//! The sample blocks' .psam rows in sample-block order. `sexes` is empty unless every block has
	//! SEX. `fids` is empty unless a block has FIDs (a block without them adds empty FIDs).
	SampleInfo samples;
	uint32_t sample_ct = 0;
	uint32_t variant_ct = 0;
	bool use_vfs = false;
	vector<PgenLocalizeGuard> localize_guards;

	//! The scanned variants, in scan order (filled by ResolveCohortSelection).
	vector<Segment> segments;
	uint32_t scan_ct = 0;

	//! The segment that holds scan position `pos` (pos < scan_ct).
	const Segment &SegmentAt(uint32_t pos) const;
	//! Variants in the scan with more than one ALT allele.
	uint64_t CountMultiallelic() const;
};

//! Load a catalog. Every cell's .pgen header is checked against its block's .pvar and .psam.
unique_ptr<CatalogCohort> LoadCatalogCohort(ClientContext &context, const string &catalog_path,
                                            const string &func_name);

//! Resolve the selection over a cohort, as ResolveSelection does over one fileset: `samples` over
//! the cohort's samples (integer indices are cohort-wide), `region` per block, and
//! `variant_range := [lo, hi)` over cohort-wide ordinals. Fills cohort.segments.
ResolvedSelection ResolveCohortSelection(const Selection &selection, CatalogCohort &cohort, const string &func_name);

//! Per-thread reader over a cohort: holds the ConcatPgenReader of the block it last read.
struct CatalogCohortReader {
	unique_ptr<ConcatPgenReader> reader;
	idx_t block = DConstants::INVALID_INDEX;

	//! The reader for `block`, opened over the selected samples. Reopens when the block changes.
	ConcatPgenReader &Open(ClientContext &context, const CatalogCohort &cohort, idx_t block, bool has_subset,
	                       const vector<uint32_t> &subset_sorted);
};

} // namespace duckdb
