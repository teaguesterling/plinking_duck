#pragma once

#include "duckdb.hpp"
#include "duckdb/common/file_system.hpp"
#include "psam_reader.hpp"

#include <pgenlib_read.h>
#include <pgenlib_ffi_support.h>
#include <pgenlib_misc.h>

#include <cstring>
#include <unordered_set>

namespace duckdb {

// ---------------------------------------------------------------------------
// Genotype output mode
// ---------------------------------------------------------------------------

//! Genotype column output mode for read_pgen / read_pfile default mode.
enum class GenotypeMode : uint8_t {
	ARRAY,   //!< ARRAY(TINYINT, N) — requires N <= MAX_ARRAY_SIZE
	LIST,    //!< LIST(TINYINT) — any sample count
	COLUMNS, //!< One scalar TINYINT column per sample (variant orient) or per variant (sample orient)
	STRUCT,  //!< STRUCT with named fields per sample/variant — like COLUMNS but single logical column
	COUNTS,  //!< STRUCT(hom_ref, het, hom_alt, missing) — fast genotype counting via PgrGetCounts
	STATS    //!< STRUCT with counts + derived stats (af, maf, missing_rate, carrier_count, het_rate)
};

//! Resolve a genotypes parameter string to a GenotypeMode.
//! Case-insensitive. "auto" → ARRAY if sample_ct ≤ MAX_ARRAY_SIZE, else LIST.
GenotypeMode ResolveGenotypeMode(const string &mode_str, uint32_t sample_ct, const string &func_name);

//! Build the STRUCT type for genotypes='counts': {hom_ref, het, hom_alt, missing}
LogicalType MakeGenotypeCountsType();

//! Build the STRUCT type for genotypes='stats': counts + derived fields
LogicalType MakeGenotypeStatsType();

//! Check if a genotype mode is an aggregate mode (counts/stats)
inline bool IsAggregateGenotypeMode(GenotypeMode mode) {
	return mode == GenotypeMode::COUNTS || mode == GenotypeMode::STATS;
}

// ---------------------------------------------------------------------------
// Orient mode
// ---------------------------------------------------------------------------

//! Orientation mode for read_pfile / read_pgen output.
enum class OrientMode : uint8_t {
	VARIANT,  //!< One row per variant (default)
	GENOTYPE, //!< One row per variant×sample
	SAMPLE    //!< One row per sample, genotypes across variants
};

//! Resolve orient parameter string to an OrientMode.
//! Errors if orient value is invalid.
OrientMode ResolveOrientMode(const string &orient_str, const string &func_name);

// ---------------------------------------------------------------------------
// RAII wrapper for cache-aligned allocations from pgenlib
// ---------------------------------------------------------------------------

//! Uses plink2::aligned_free() which expects the aligned_malloc header.
struct AlignedBuffer {
	void *ptr = nullptr;

	~AlignedBuffer() {
		if (ptr) {
			plink2::aligned_free(ptr);
		}
	}

	AlignedBuffer() = default;
	AlignedBuffer(const AlignedBuffer &) = delete;
	AlignedBuffer &operator=(const AlignedBuffer &) = delete;
	AlignedBuffer(AlignedBuffer &&other) noexcept : ptr(other.ptr) {
		other.ptr = nullptr;
	}
	AlignedBuffer &operator=(AlignedBuffer &&other) noexcept {
		if (this != &other) {
			if (ptr) {
				plink2::aligned_free(ptr);
			}
			ptr = other.ptr;
			other.ptr = nullptr;
		}
		return *this;
	}

	//! Allocate a cache-aligned buffer of the given size in bytes.
	void Allocate(uintptr_t size) {
		if (plink2::cachealigned_malloc(size, &ptr)) {
			throw IOException("failed to allocate %llu bytes of aligned memory", static_cast<unsigned long long>(size));
		}
	}

	template <typename T>
	T *As() {
		return static_cast<T *>(ptr);
	}

	template <typename T>
	const T *As() const {
		return static_cast<const T *>(ptr);
	}
};

// ---------------------------------------------------------------------------
// Offset-indexed variant metadata (memory-efficient Scan-time access)
// ---------------------------------------------------------------------------

//! Columnar variant metadata index.
//!
//! Two modes:
//!  * Dense  (vidx_map empty): vectors are indexed directly by file-row vidx.
//!  * Sparse (vidx_map non-empty): vectors hold a filtered subset; vidx_map
//!    maps file-row vidx → local index. `variant_ct` still reflects the
//!    total source row count (for count-mismatch validation), while
//!    `chroms.size()` is the loaded subset size.
//!
//! All accessors take a file-row vidx. Thread-safe for concurrent reads
//! once construction (or lazy EnsureAlleles/EnsureIds) has completed.
struct VariantMetadataIndex {
	//! Columnar typed backing
	vector<string> chroms;
	vector<int32_t> positions;
	vector<string> ids; // "." normalized to "" at load time
	vector<string> refs;
	vector<string> alts; // "." normalized to "" at load time

	//! Sparse mode: file-row vidx → local index in the vectors. Empty in dense mode.
	unordered_map<uint32_t, uint32_t> vidx_map;

	//! Fast path for region-pushdown subsets whose file-row variant indices are
	//! contiguous. Avoids an unordered_map lookup per GetChrom()/GetPos() call in
	//! dense LD window scans.
	bool has_contiguous_vidx_range = false;
	uint32_t contiguous_start_vidx = 0;

	//! Sparse mode: local index → file-row vidx (reverse of vidx_map).
	//! Empty in dense mode (where local index == file-row vidx).
	//! Keeping both directions lets Local()/GetChrom()/etc. stay O(1) AND lets
	//! scan/filter code look up "what file-row vidx is this local row?" in O(1)
	//! (used by ResolveByCpra and BuildVariantIdIndex in sparse mode).
	vector<uint32_t> local_to_vidx;

	//! Chromosome offset index for fast region lookups (dense mode only).
	//! chrom → [first_local_idx, past_end_local_idx). Assumes variants are
	//! (CHROM, POS)-sorted as required by the PLINK spec.
	unordered_map<string, std::pair<idx_t, idx_t>> chrom_offsets;

	//! Total variant count in the source file (not the loaded subset).
	idx_t variant_ct = 0;

	//! Whether the source was .bim format (affects column-name defaults only).
	bool is_bim = false;

	//! Whether ID/REF/ALT are loaded. Dense text load always sets both true.
	//! Parquet lazy paths may skip loading these if the query doesn't need them.
	bool has_ids = true;
	bool has_alleles = true;

	//! Local index for a file-row vidx.
	inline idx_t Local(idx_t vidx) const {
		if (has_contiguous_vidx_range) {
			if (vidx >= contiguous_start_vidx && vidx < contiguous_start_vidx + chroms.size()) {
				return vidx - contiguous_start_vidx;
			}
			throw InternalException("VariantMetadataIndex: vidx %llu not in contiguous loaded subset",
			                        static_cast<unsigned long long>(vidx));
		}
		if (vidx_map.empty()) {
			return vidx;
		}
		auto it = vidx_map.find(static_cast<uint32_t>(vidx));
		if (it == vidx_map.end()) {
			throw InternalException("VariantMetadataIndex: vidx %llu not in loaded subset",
			                        static_cast<unsigned long long>(vidx));
		}
		return it->second;
	}

	//! True iff this index is dense (fully loaded, indexed by vidx directly).
	inline bool IsDense() const {
		return vidx_map.empty() && !has_contiguous_vidx_range;
	}

	//! Returns [first_local_idx, past_end_local_idx) for a chromosome; {0,0} if absent.
	//! Dense mode only — sparse indexes have been pre-filtered.
	inline std::pair<idx_t, idx_t> ChromRange(const string &chrom) const {
		auto it = chrom_offsets.find(chrom);
		return it == chrom_offsets.end() ? std::make_pair(idx_t {0}, idx_t {0}) : it->second;
	}

	inline const string &GetChrom(idx_t vidx) const {
		return chroms[Local(vidx)];
	}
	inline int32_t GetPos(idx_t vidx) const {
		return positions[Local(vidx)];
	}
	inline const string &GetId(idx_t vidx) const {
		return ids[Local(vidx)];
	}
	inline const string &GetRef(idx_t vidx) const {
		return refs[Local(vidx)];
	}
	inline const string &GetAlt(idx_t vidx) const {
		return alts[Local(vidx)];
	}

	//! Map a local index (into chroms/positions/ids/refs/alts) back to a file-row vidx.
	//! O(1) in both modes. Use this when iterating the loaded subset (e.g. all of
	//! chroms.size()) and you need the pgenlib-compatible vidx for each row.
	inline uint32_t VidxForLocal(idx_t local) const {
		if (has_contiguous_vidx_range) {
			return contiguous_start_vidx + static_cast<uint32_t>(local);
		}
		if (local_to_vidx.empty()) {
			return static_cast<uint32_t>(local); // dense: local == vidx
		}
		return local_to_vidx[local];
	}
};

//! Build a map from variant ID (rsid) to file-row vidx. Handles both dense and
//! sparse variant indexes: in sparse mode only the loaded subset is indexed
//! (which is the correct semantic — user-supplied rsids can only resolve to
//! variants that are actually in the loaded subset, typically a region).
//! Empty IDs are skipped.
unordered_map<string, uint32_t> BuildVariantIdIndex(const VariantMetadataIndex &variants);

//! Build an offset-indexed metadata index from a .pvar/.bim file.
//! Reads the file once into a single buffer, builds line offsets, and parses
//! the header to determine column layout. Does NOT parse data lines.
VariantMetadataIndex LoadVariantMetadataIndex(ClientContext &context, const string &path, const string &func_name);

// ---------------------------------------------------------------------------
// File utilities
// ---------------------------------------------------------------------------

//! Read an entire file via DuckDB's VFS and split into lines.
vector<string> ReadFileLines(ClientContext &context, const string &path);

//! Split a line on tab characters.
vector<string> SplitTabLine(const string &line);

//! Split a line on whitespace (for .bim format).
vector<string> SplitWhitespaceLine(const string &line);

// ---------------------------------------------------------------------------
// Companion file discovery
// ---------------------------------------------------------------------------

//! Replace the extension of a file path.
string ReplaceExtension(const string &path, const string &new_ext);

//! Try to find a companion file by replacing the .pgen extension.
//! Returns the first existing path from candidates, or empty string.
string FindCompanionFile(FileSystem &fs, const string &pgen_path, const vector<string> &extensions);

//! Check if a file path refers to a parquet file.
bool IsParquetFile(const string &path);

//! Check if a path refers to a native PLINK format file.
//! Returns true for .pvar, .bim, .psam, .fam, .pvar.zst, .psam.zst extensions.
bool IsNativePlinkFormat(const string &path);

//! Read the plinking_use_parquet_companions config option.
bool GetUseParquetCompanions(ClientContext &context);

//! Find companion file with parquet priority when config is enabled.
//! For variant files: checks .pvar.parquet before .pvar/.bim
//! For sample files: checks .psam.parquet before .psam/.fam
string FindCompanionFileWithParquet(ClientContext &context, FileSystem &fs, const string &pgen_path,
                                    const vector<string> &extensions);

//! Load variant metadata from a parquet file via DuckDB's parquet reader.
//! Uses a separate Connection to avoid bind-phase reentrancy.
VariantMetadataIndex LoadVariantMetadataFromParquet(ClientContext &context, const string &path,
                                                    const string &func_name);

//! Region-pushdown parquet loader — materializes only variants in [pos_start, pos_end]
//! on the named chrom. Returns a sparse VariantMetadataIndex whose vidx_map keys are
//! file row numbers (pgenlib-compatible global vidx).
//!
//! `variant_ct_hint` sets `idx.variant_ct`. Callers should pass the expected total
//! (typically `pgfi.raw_variant_ct`) so the downstream count-mismatch check is a
//! sanity check against pgen — we deliberately do NOT re-scan the full parquet to
//! count rows, which would defeat the point of pushdown at 170M-variant scale.
VariantMetadataIndex LoadVariantMetadataFromParquetRegion(ClientContext &context, const string &path,
                                                          const string &chrom, int64_t pos_start, int64_t pos_end,
                                                          idx_t variant_ct_hint, const string &func_name);

//! Region-pushdown text loader for native .pvar/.bim companions. Scans only
//! until the sorted target interval has been passed and materializes only
//! matching rows. Returns a sparse VariantMetadataIndex keyed by global
//! pgenlib-compatible file-row variant indices.
VariantMetadataIndex LoadVariantMetadataFromTextRegion(ClientContext &context, const string &path,
                                                       const string &chrom, int64_t pos_start, int64_t pos_end,
                                                       idx_t variant_ct_hint, const string &func_name);

//! Row count from a parquet file via row-group metadata aggregation
//! (DuckDB-optimized `COUNT(*)`). Typically sub-ms but O(num_row_groups),
//! not literally O(1). Use for count-only metadata paths (psam fast path);
//! the region-pushdown path avoids this entirely by trusting pgen's count.
idx_t GetParquetRowCount(ClientContext &context, const string &path);

//! Load sample info from a parquet file via DuckDB's parquet reader.
SampleInfo LoadSampleInfoFromParquet(ClientContext &context, const string &path);

//! Load variant metadata from any DuckDB-readable source (CSV, table, view, etc.).
//! Executes SELECT * FROM '<source>' via a separate Connection, maps columns by name
//! (case-insensitive), synthesizes a pvar-format text buffer, and returns a standard
//! VariantMetadataIndex.
VariantMetadataIndex LoadVariantMetadataFromSource(ClientContext &context, const string &source,
                                                   const string &func_name);

//! Load sample info from any DuckDB-readable source (CSV, table, view, etc.).
//! Same approach: query via separate Connection, map columns, synthesize psam-format buffer,
//! return standard SampleInfo.
SampleInfo LoadSampleInfoFromSource(ClientContext &context, const string &source);

//! Load variant metadata, auto-dispatching between parquet, native text, and arbitrary sources.
VariantMetadataIndex LoadVariantMetadata(ClientContext &context, const string &path, const string &func_name);

//! Load sample info, auto-dispatching between parquet, native text, and arbitrary sources.
SampleInfo LoadSampleMetadata(ClientContext &context, const string &path);

//! Cheap count-only psam loader: returns SampleInfo with `sample_ct` set but
//! `iids`/`fids`/`iid_to_idx` empty. For parquet, uses footer metadata (O(1)).
//! For text/source, falls back to full load and discards strings to keep memory low.
//! Use when the query doesn't need IID strings (no VARCHAR sample filter, no
//! genotypes:='columns'/'struct', not orient='genotype'/'sample'). At biobank
//! scale this saves the ~600ms columnar string copy of 7M IIDs.
SampleInfo LoadSampleCount(ClientContext &context, const string &path);

// ---------------------------------------------------------------------------
// Sample subsetting (for PgrGetCounts / PgrGet)
// ---------------------------------------------------------------------------

//! Resolve the `samples` named parameter into 0-based sample indices.
//! Handles both LIST(INTEGER) and LIST(VARCHAR) dispatch.
vector<uint32_t> ResolveSampleIndices(const Value &samples_val, uint32_t raw_sample_ct, const SampleInfo *sample_info,
                                      const string &func_name);

//! Sample subsetting buffers shared across all scan threads.
//! Built once in bind, provides sample_include bitmask, interleaved vec
//! (for PgrGetCounts), and cumulative popcounts (for PgrSampleSubsetIndex).
struct SampleSubset {
	AlignedBuffer sample_include_buf;
	AlignedBuffer interleaved_vec_buf;
	AlignedBuffer cumulative_popcounts_buf;
	uint32_t subset_sample_ct = 0;
	uint32_t raw_sample_ct = 0;

	SampleSubset() = default;
	SampleSubset(SampleSubset &&) = default;
	SampleSubset &operator=(SampleSubset &&) = default;

	const uintptr_t *SampleInclude() const {
		return sample_include_buf.As<uintptr_t>();
	}
	const uintptr_t *InterleavedVec() const {
		return interleaved_vec_buf.As<uintptr_t>();
	}
	const uint32_t *CumulativePopcounts() const {
		return cumulative_popcounts_buf.As<uint32_t>();
	}
};

//! Build a SampleSubset from 0-based sample indices.
SampleSubset BuildSampleSubset(uint32_t raw_sample_ct, const vector<uint32_t> &sample_indices);

// ---------------------------------------------------------------------------
// Region filtering
// ---------------------------------------------------------------------------

//! Parsed variant range from a region string.
struct VariantRange {
	uint32_t start_idx = 0; //!< First matching variant index (inclusive)
	uint32_t end_idx = 0;   //!< Past-the-end matching variant index (exclusive)
	bool has_filter = false;
};

//! Parse "chr:start-end" region string and resolve to variant index range
//! using on-demand field parsing from VariantMetadataIndex.
//! Assumes variants are sorted by (CHROM, POS) as required by PLINK format.
//! Returns a range with start_idx == end_idx if no variants match.
VariantRange ParseRegion(const string &region_str, const VariantMetadataIndex &variants, const string &func_name);

// ---------------------------------------------------------------------------
// Count-based filtering (af_range, ac_range)
// ---------------------------------------------------------------------------

struct RangeFilter {
	double min = -std::numeric_limits<double>::infinity();
	double max = std::numeric_limits<double>::infinity();
	bool active = false;
	bool Passes(double value) const {
		return value >= min && value <= max;
	}
};

struct CountFilter {
	RangeFilter af_filter;
	RangeFilter ac_filter;
	bool HasFilter() const {
		return af_filter.active || ac_filter.active;
	}
};

//! Parse a STRUCT value into a RangeFilter with optional min/max fields.
//! valid_min/valid_max define the legal range for values (e.g. 0.0-1.0 for AF).
//! When include_missing_out is non-null, an "include_missing" BOOLEAN field is
//! also accepted and written out (used by genotype_range); otherwise that field
//! name is rejected like any other unknown field.
RangeFilter ParseRangeFilter(const Value &val, const string &param_name, double valid_min, double valid_max,
                             const string &func_name, bool *include_missing_out = nullptr);

//! Check if a variant passes count-based filters given its genotype counts.
//! genocounts: [hom_ref, het, hom_alt, missing] from PgrGetCounts.
bool VariantPassesCountFilter(const CountFilter &filter, const STD_ARRAY_REF(uint32_t, 4) genocounts,
                              uint32_t effective_sample_ct);

// ---------------------------------------------------------------------------
// Genotype range filtering (genotype_range)
// ---------------------------------------------------------------------------

struct GenotypeRangeResult {
	bool any_pass; // at least one non-missing sample has value in [min, max]
	bool all_pass; // every non-missing sample has value in [min, max]
};

//! Selects which hardcall genotype categories a genotype filter keeps.
//! allowed[0..2] gate hom_ref / het / hom_alt (ALT dosage 0/1/2); include_missing
//! gates missing (-9). Populated from either the `genotype_range` (numeric alias)
//! or `include_genotypes` (named-category) parameter — both compile to this mask,
//! so all engine sites use a single predicate.
struct GenotypeRangeFilter {
	bool active = false;
	bool allowed[3] = {false, false, false}; // hom_ref, het, hom_alt
	bool include_missing = false;

	//! True when a hardcall ALT dosage (0/1/2, or a phased diploid sum) is a
	//! selected category. Non-integer or out-of-domain values never match. Takes a
	//! double so it slots in where the old RangeFilter::Passes(double) was called.
	bool AllowsCall(double dosage) const {
		int v = static_cast<int>(dosage);
		if (v < 0 || v > 2 || static_cast<double>(v) != dosage) {
			return false;
		}
		return allowed[v];
	}

	//! Populate the mask from a numeric [min,max] range over {0,1,2} (the
	//! `genotype_range` alias). include_missing comes from the range's own field.
	void SetFromRange(const RangeFilter &range, bool inc_missing) {
		for (int g = 0; g <= 2; g++) {
			allowed[g] = range.Passes(static_cast<double>(g));
		}
		include_missing = inc_missing;
		active = range.active;
	}
};

GenotypeRangeResult CheckGenotypeRange(const GenotypeRangeFilter &filter, const STD_ARRAY_REF(uint32_t, 4) genocounts);

//! Parse an `include_genotypes` LIST of category labels into a GenotypeRangeFilter.
//! Accepts 'hom_ref', 'het', 'hom_alt', 'missing' (case-insensitive). Empty list
//! is a no-op (inactive). Throws on unknown labels.
void ParseIncludeGenotypes(const Value &val, GenotypeRangeFilter &out, const string &func_name);

//! Combined pre-decompression filter check: count filter + genotype range.
//! Returns: skip=true → skip variant entirely. skip=false → variant passes,
//! all_pass indicates whether per-element range check can be skipped.
struct PreDecompFilterResult {
	bool skip;     //!< variant should be skipped entirely
	bool all_pass; //!< all non-missing genotypes pass range filter (skip per-element check)
};

PreDecompFilterResult CheckPreDecompFilters(const CountFilter &count_filter, const GenotypeRangeFilter &genotype_filter,
                                            const STD_ARRAY_REF(uint32_t, 4) genocounts, uint32_t sample_ct);

// ---------------------------------------------------------------------------
// Genotype normalization for PCA (Price et al. 2006)
// ---------------------------------------------------------------------------

//! Per-variant normalization parameters for standardized genotype matrix.
//! p = ALT allele frequency (matches plink_freq ALT_FREQ convention).
struct VariantNorm {
	double center;    //!< 2 * allele_freq
	double inv_stdev; //!< 1.0 / sqrt(2 * af * (1 - af))
	bool skip;        //!< true if monomorphic (af == 0 or af == 1)
};

//! Compute normalization parameters from ALT allele frequency.
VariantNorm ComputeVariantNorm(double alt_freq);

//! Normalize int8 genotypes (from GenoarrToBytesMinus9) into centered, standardized doubles.
//! Missing values (-9) become 0.0 (mean imputation after centering).
//! output must be pre-allocated to sample_ct doubles.
void NormalizeGenotypes(const int8_t *genotypes, uint32_t sample_ct, const VariantNorm &norm, double *output);

// ---------------------------------------------------------------------------
// Phased genotype unpacking
// ---------------------------------------------------------------------------

//! Unpack genotype bytes + phase bitarrays into haplotype pairs.
//! Takes already-unpacked genotype bytes (from GenoarrToBytesMinus9) plus
//! phasepresent/phaseinfo bitarrays from PgrGetP. Writes sample_ct * 2
//! int8_t values to output_pairs with stride-2 layout:
//!   [allele1_s0, allele2_s0, allele1_s1, allele2_s1, ...]
//! Missing genotypes (-9) produce [-9, -9] pairs; caller handles NULL at vector level.
void UnpackPhasedGenotypes(const int8_t *genotype_bytes, const uintptr_t *phasepresent, const uintptr_t *phaseinfo,
                           uint32_t sample_ct, int8_t *output_pairs);

// ---------------------------------------------------------------------------
// Shared genotype output vector filling
// ---------------------------------------------------------------------------

//! Fill a genotype output vector for one row in ARRAY or LIST mode.
//! Handles unphased (genotype_bytes) and phased (phased_pairs) output.
//! For unphased: writes int8 genotype values, -9 → NULL.
//! For phased: writes ARRAY(TINYINT, 2) pairs, -9 → NULL pair.
//! Does NOT handle COLUMNS, STRUCT, COUNTS, STATS, or dosage modes.
void FillGenotypeVector(Vector &vec, idx_t row_idx, GenotypeMode mode, uint32_t output_sample_ct,
                        const int8_t *genotype_bytes, const int8_t *phased_pairs, bool include_phased);

// ---------------------------------------------------------------------------
// Unified variants parameter resolution
// ---------------------------------------------------------------------------

//! Resolve a variants parameter value into a sorted list of 0-based variant indices.
//! Handles all supported input types: single integer, single string (rsid or CPRA),
//! CPRA struct, range struct, and lists of any of those.
vector<uint32_t> ResolveVariantsParameter(const Value &val, const VariantMetadataIndex &variants,
                                          uint32_t raw_variant_ct, const string &func_name);

// ---------------------------------------------------------------------------
// Ploidy- and sex-aware statistics for sex/organelle chromosomes (chrX/Y/MT)
// ---------------------------------------------------------------------------

//! Ploidy class of a variant, determined by chromosome (and X position vs PAR).
enum class ChromPloidy : uint8_t {
	AUTOSOMAL, //!< Diploid for all samples (autosomes AND the chrX pseudo-autosomal region)
	CHR_X,     //!< chrX non-PAR: females diploid, males haploid, unknown-sex excluded
	CHR_Y,     //!< chrY: males haploid; females and unknown-sex excluded
	CHR_MT     //!< mitochondrial: haploid for everyone (sex-independent)
};

//! Pseudo-autosomal region boundaries for a genome build (1-based, inclusive).
//! chrX positions in [1, par1_end] (PAR1) or [par2_start, par2_end] (PAR2) are
//! diploid in both sexes and thus classified AUTOSOMAL.
struct ParBounds {
	int32_t par1_end = 0;
	int32_t par2_start = 0;
	int32_t par2_end = 0;
	bool active = false; //!< false → no coordinate-based PAR detection
};

//! Resolve a genome-build string to PAR boundaries. Recognizes GRCh38/hg38/b38,
//! GRCh37/hg19/b37, and 'none'/'' (disables coordinate PAR detection so PAR is
//! only recognized via explicit chromosome codes PAR1/PAR2/XY/25). Throws on an
//! unrecognized non-empty value.
ParBounds ResolveParBounds(const string &build, const string &func_name);

//! Classify a variant's ploidy from its chromosome label and position.
//! Chromosome matching is case-insensitive and ignores a leading "chr".
//! Explicit PAR codes (PAR1/PAR2/XY/25) → AUTOSOMAL. For X/23, positions within
//! `par` boundaries → AUTOSOMAL, otherwise CHR_X.
ChromPloidy ClassifyChromPloidy(const string &chrom, int32_t pos, const ParBounds &par);

//! Aggregated ploidy/sex-aware counts for one variant.
struct SexAwareCounts {
	//! Allele-based tallies for MAF/ALT_FREQ (haploid samples contribute 1 allele).
	uint32_t alt_allele_ct = 0;
	uint32_t obs_allele_ct = 0; //!< non-missing allele count (denominator)

	//! Diploid genotype counts of the HWE stratum: females for CHR_X. All zero and
	//! hwe_defined=false for CHR_Y/CHR_MT (HWE undefined on haploid data).
	uint32_t hwe_hom_ref = 0;
	uint32_t hwe_het = 0;
	uint32_t hwe_hom_alt = 0;
	bool hwe_defined = false;

	//! Genotype-distribution counts over included samples (haploid samples counted
	//! as homozygous). Feeds plink_freq's optional count columns.
	uint32_t geno_hom_ref = 0;
	uint32_t geno_het = 0;
	uint32_t geno_hom_alt = 0;
	uint32_t geno_missing = 0;

	//! True when the chromosome needs sex to stratify (CHR_X/CHR_Y) but no sex
	//! information was available → all derived stats should be emitted as NULL.
	bool sex_unavailable = false;
};

//! Compute ploidy/sex-aware counts from decoded hardcall bytes.
//! `geno_bytes[i]` ∈ {0,1,2} (ALT dosage) or -9 (missing), in effective-sample
//! order. `sex[i]` ∈ {0,1,2} aligned to the same order. `have_sex` indicates
//! whether `sex` is populated. On haploid strata a heterozygous hardcall
//! (value 1) is treated as missing (plink2 convention). Never called for
//! AUTOSOMAL — callers keep the fast PgrGetCounts path there.
SexAwareCounts ComputeSexAwareCounts(const int8_t *geno_bytes, uint32_t sample_ct, ChromPloidy ploidy,
                                     const uint8_t *sex, bool have_sex);

//! Build a sex array aligned to effective (post-subset) sample order.
//! Returns empty if the source has no sex data. `subset_sorted` lists original
//! sample indices in ascending order when a sample subset is active (matching
//! pgenlib's PgrGet ordering); pass nullptr for the full cohort.
vector<uint8_t> BuildAlignedSex(const SampleInfo &sample_info, const vector<uint32_t> *subset_sorted);

// ---------------------------------------------------------------------------
// Max threads config helper
// ---------------------------------------------------------------------------

//! Read the plinking_max_threads config option from the client context.
//! Returns 0 if unset/default (meaning "no override").
uint32_t GetPlinkingMaxThreads(ClientContext &context);

//! Apply the max threads cap to a computed thread count.
//! If config_max_threads > 0, returns min(computed, config_max_threads).
//! Otherwise returns min(computed, 16) (existing default behavior).
idx_t ApplyMaxThreadsCap(idx_t computed, uint32_t config_max_threads);

} // namespace duckdb
