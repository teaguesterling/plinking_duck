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

//! Offset-indexed variant metadata for memory-efficient Scan-time access.
//! Stores raw file content in a single buffer and line byte offsets; parses
//! fields on demand. Thread-safe for concurrent reads (all state is immutable
//! after construction).
struct VariantMetadataIndex {
	//! Raw file content (single allocation)
	string file_content;

	//! Byte offset of each data line's start within file_content.
	//! line_offsets[vidx] = byte offset of variant vidx's line.
	vector<uint64_t> line_offsets;

	//! Whether the source is .bim format (whitespace-delimited, different column order)
	bool is_bim = false;

	//! Physical field indices for each logical field (after header parsing)
	idx_t chrom_idx = 0;
	idx_t pos_idx = 0;
	idx_t id_idx = 0;
	idx_t ref_idx = 0;
	idx_t alt_idx = 0;

	//! Total variant count
	idx_t variant_ct = 0;

	//! Find the end of line vidx (exclusive, past trailing \r\n)
	size_t LineEnd(idx_t vidx) const;

	//! Extract the N-th delimited field from a line without allocating.
	//! For .pvar: tab-delimited. For .bim: whitespace-delimited.
	string GetField(idx_t vidx, idx_t field_idx) const;

	//! On-demand field access (thread-safe: const on file_content)
	string GetChrom(idx_t vidx) const;
	int32_t GetPos(idx_t vidx) const;
	string GetId(idx_t vidx) const;
	string GetRef(idx_t vidx) const;
	string GetAlt(idx_t vidx) const;
};

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

//! Load sample info from a parquet file via DuckDB's parquet reader.
SampleInfo LoadSampleInfoFromParquet(ClientContext &context, const string &path);

//! Load variant metadata from any DuckDB-readable source (CSV, table, view, etc.).
//! Executes SELECT * FROM '<source>' via a separate Connection, maps columns by name
//! (case-insensitive), synthesizes a pvar-format text buffer, and returns a standard
//! VariantMetadataIndex.
VariantMetadataIndex LoadVariantMetadataFromSource(ClientContext &context, const string &source, const string &func_name);

//! Load sample info from any DuckDB-readable source (CSV, table, view, etc.).
//! Same approach: query via separate Connection, map columns, synthesize psam-format buffer,
//! return standard SampleInfo.
SampleInfo LoadSampleInfoFromSource(ClientContext &context, const string &source);

//! Load variant metadata, auto-dispatching between parquet, native text, and arbitrary sources.
VariantMetadataIndex LoadVariantMetadata(ClientContext &context, const string &path, const string &func_name);

//! Load sample info, auto-dispatching between parquet, native text, and arbitrary sources.
SampleInfo LoadSampleMetadata(ClientContext &context, const string &path);

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
RangeFilter ParseRangeFilter(const Value &val, const string &param_name, double valid_min, double valid_max,
                             const string &func_name);

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

struct GenotypeRangeFilter {
	RangeFilter range; // reuses Phase 1 RangeFilter
	bool active = false;
};

GenotypeRangeResult CheckGenotypeRange(const RangeFilter &filter, const STD_ARRAY_REF(uint32_t, 4) genocounts);

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
// Unified variants parameter resolution
// ---------------------------------------------------------------------------

//! Resolve a variants parameter value into a sorted list of 0-based variant indices.
//! Handles all supported input types: single integer, single string (rsid or CPRA),
//! CPRA struct, range struct, and lists of any of those.
vector<uint32_t> ResolveVariantsParameter(const Value &val, const VariantMetadataIndex &variants,
                                           uint32_t raw_variant_ct, const string &func_name);

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
