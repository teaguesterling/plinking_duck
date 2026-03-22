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
	ARRAY,  //!< ARRAY(TINYINT, N) — requires N <= MAX_ARRAY_SIZE
	LIST,   //!< LIST(TINYINT) — any sample count
	COLUMNS //!< One scalar TINYINT column per sample (variant orient) or per variant (sample orient)
};

//! Resolve a genotypes parameter string to a GenotypeMode.
//! Case-insensitive. "auto" → ARRAY if sample_ct ≤ MAX_ARRAY_SIZE, else LIST.
GenotypeMode ResolveGenotypeMode(const string &mode_str, uint32_t sample_ct, const string &func_name);

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
// Parallel text file scanning (shared by pvar/psam/bim/fam readers)
// ---------------------------------------------------------------------------

//! Chunk size for parallel text file reading (~4MB per chunk).
//! Chosen to balance parallelism granularity against per-chunk overhead.
static constexpr uint64_t TEXT_FILE_CHUNK_SIZE = 4 * 1024 * 1024;

//! Global state for parallel text file scanning.
//! Manages byte-range chunk distribution across threads via atomic claiming.
struct TextFileGlobalState : public GlobalTableFunctionState {
	string file_path;
	uint64_t file_size = 0;
	uint64_t data_start_offset = 0;
	std::atomic<uint64_t> next_chunk_offset {0};
	vector<column_t> column_ids;

	idx_t MaxThreads() const override {
		uint64_t data_size = file_size > data_start_offset ? file_size - data_start_offset : 0;
		if (data_size < TEXT_FILE_CHUNK_SIZE) {
			return 1;
		}
		return std::min<idx_t>(data_size / TEXT_FILE_CHUNK_SIZE + 1, 16);
	}
};

//! Per-thread state for parallel text file scanning.
//! Each thread owns a file handle, reads chunks of the file, and extracts
//! complete lines using the byte-range boundary protocol:
//!   - If not the first chunk, skip to the next newline (previous thread owns partial line)
//!   - Read past chunk end to finish the current line
struct TextFileLocalState : public LocalTableFunctionState {
	unique_ptr<FileHandle> handle;

	//! Read buffer
	static constexpr size_t READ_BUF_SIZE = 64 * 1024;
	char read_buf[READ_BUF_SIZE];
	size_t buf_pos = 0;   //!< Current position within read_buf
	size_t buf_len = 0;   //!< Valid bytes in read_buf

	//! Current chunk boundaries (byte offsets in file)
	uint64_t chunk_start = 0;
	uint64_t chunk_end = 0;
	uint64_t file_pos = 0;   //!< Current byte offset in file (next byte to consume from buffer)
	bool chunk_exhausted = true;

	//! Initialize file handle for this thread
	void Init(ClientContext &context, const string &file_path) {
		auto &fs = FileSystem::GetFileSystem(context);
		handle = fs.OpenFile(file_path, FileFlags::FILE_FLAGS_READ);
	}

	//! Claim the next chunk from the global state. Returns false if no more work.
	bool ClaimChunk(TextFileGlobalState &gstate) {
		uint64_t start = gstate.next_chunk_offset.fetch_add(TEXT_FILE_CHUNK_SIZE);
		if (start >= gstate.file_size) {
			return false;
		}
		chunk_start = std::max(start, gstate.data_start_offset);
		chunk_end = std::min(start + TEXT_FILE_CHUNK_SIZE, gstate.file_size);

		// Seek to chunk start and reset buffer
		handle->Seek(chunk_start);
		file_pos = chunk_start;
		buf_pos = 0;
		buf_len = 0;
		chunk_exhausted = false;

		// If not the first data chunk, skip to the next newline
		// (the previous thread will finish the partial line at this boundary)
		if (chunk_start > gstate.data_start_offset) {
			string discard;
			ReadLine(discard, gstate.file_size);
		}
		return true;
	}

	//! Read one line from the buffered file handle, returning false at EOF.
	//! Handles \r\n and \n line endings (strips \r).
	//! file_size_limit prevents reading past end of file.
	bool ReadLine(string &line, uint64_t file_size_limit) {
		line.clear();
		while (true) {
			// Refill buffer if needed
			if (buf_pos >= buf_len) {
				if (file_pos >= file_size_limit) {
					return !line.empty();
				}
				size_t to_read = std::min<size_t>(READ_BUF_SIZE, file_size_limit - file_pos);
				buf_len = handle->Read(read_buf, to_read);
				buf_pos = 0;
				if (buf_len == 0) {
					return !line.empty();
				}
			}

			char c = read_buf[buf_pos++];
			file_pos++;
			if (c == '\n') {
				// Strip trailing \r
				if (!line.empty() && line.back() == '\r') {
					line.pop_back();
				}
				return true;
			}
			line += c;
		}
	}

	//! Read lines for this thread's current chunk, respecting boundaries.
	//! Returns false when the chunk (and any overflow line) is exhausted.
	//! When file_pos >= chunk_end, finishes the current line then stops.
	bool ReadChunkLine(string &line, uint64_t file_size) {
		if (chunk_exhausted) {
			return false;
		}
		// If we're past our chunk end, we're done (the previous call
		// already finished the line that straddled the boundary)
		if (file_pos >= chunk_end) {
			chunk_exhausted = true;
			return false;
		}
		// Read the line — allow reading past chunk_end to finish it
		bool got_line = ReadLine(line, file_size);
		if (!got_line) {
			chunk_exhausted = true;
		}
		return got_line;
	}
};

} // namespace duckdb
