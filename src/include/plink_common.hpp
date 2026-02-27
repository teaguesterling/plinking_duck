#pragma once

#include "duckdb.hpp"
#include "duckdb/common/file_system.hpp"
#include "pvar_reader.hpp"
#include "psam_reader.hpp"

#include <pgenlib_read.h>
#include <pgenlib_ffi_support.h>
#include <pgenlib_misc.h>

#include <atomic>
#include <cstring>
#include <unordered_set>

namespace duckdb {

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
			throw IOException("failed to allocate %llu bytes of aligned memory",
			                  static_cast<unsigned long long>(size));
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
// Variant metadata
// ---------------------------------------------------------------------------

//! Pre-loaded variant metadata from .pvar/.bim, in file order.
struct VariantMetadata {
	vector<string> chroms;
	vector<int32_t> positions;
	vector<string> ids;
	vector<string> refs;
	vector<string> alts;
	idx_t variant_ct = 0;
};

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
// Metadata loading
// ---------------------------------------------------------------------------

//! Load variant metadata from a .pvar or .bim file.
//! func_name is used in error messages (e.g., "read_pgen", "plink_freq").
VariantMetadata LoadVariantMetadata(ClientContext &context, const string &path, const string &func_name);

// ---------------------------------------------------------------------------
// Sample subsetting (for PgrGetCounts / PgrGet)
// ---------------------------------------------------------------------------

//! Resolve the `samples` named parameter into 0-based sample indices.
//! Handles both LIST(INTEGER) and LIST(VARCHAR) dispatch.
vector<uint32_t> ResolveSampleIndices(const Value &samples_val, uint32_t raw_sample_ct,
                                      const SampleInfo *sample_info, const string &func_name);

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

//! Parse "chr:start-end" region string and resolve to variant index range.
//! Assumes variants are sorted by (CHROM, POS) as required by PLINK format.
//! Returns a range with start_idx == end_idx if no variants match.
VariantRange ParseRegion(const string &region_str, const VariantMetadata &variants, const string &func_name);

} // namespace duckdb
