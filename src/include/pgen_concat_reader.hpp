// ConcatPgenReader: read a SAMPLE-SHARDED cohort (same variants, disjoint samples, one
// .pgen per shard) as if it were one .pgen whose samples are the shards' samples
// concatenated in shard order. Used by read_pfile's combine_samples := 'concatenate'
// (P13-008, docs/planning/P13-008-sample-concat-reader).
//
// Each Get* method has the same output contract as the pgenlib call it mirrors
// (PgrGet / PgrGetCounts / PgrGetD / PgrGetP), but over the COMBINED sample space
// [0, Σ shard sample_ct): the caller passes no sample_include/pssi — the reader owns
// one per shard, split from the combined subset given at Open. Every shard holds the
// same variant at the same vidx (positional trust; the variant COUNT is checked at
// bind), so one vidx addresses all shards.
//
// One reader per scan thread (a PgenReader's FILE* is not thread-safe).
#pragma once

#include "pgen_header_cache.hpp"
#include "duckdb.hpp"
#include "plink_common.hpp"

namespace duckdb {

//! One shard of a sample-sharded cohort.
struct PgenSampleShard {
	string pgen_path;
	uint32_t raw_sample_ct = 0; //! this .pgen's sample count (header)
	uint32_t sample_offset = 0; //! first combined-space sample index of this shard
};

//! Per-call tally of how each shard decoded under GetDifflistOrGenovec (P13-011).
//! A test reads these to PROVE both fallback branches ran in the concat reader
//! under test — the on-disk record encoding, not the data distribution, decides
//! which one fires, so a single-file stand-in could not prove it.
struct ConcatSparseStats {
	uint32_t difflist_shards = 0; //! shards decoded hom_ref-major (difflist consumed)
	uint32_t dense_shards = 0;    //! shards that returned a dense genovec (common_geno == UINT32_MAX)
	uint32_t altmajor_shards = 0; //! shards with an ALT/missing-major difflist (common_geno in {1,2,3}), re-decoded
};

class ConcatPgenReader {
public:
	ConcatPgenReader() = default;
	~ConcatPgenReader();
	ConcatPgenReader(const ConcatPgenReader &) = delete;
	ConcatPgenReader &operator=(const ConcatPgenReader &) = delete;

	//! Open every shard. `subset_sorted` holds sorted-unique COMBINED-space sample
	//! indices when `has_subset`, else all samples are read. `use_vfs` routes the
	//! .pgen opens through DuckDB's VFS (Path V) while opening.
	void Open(ClientContext &context, bool use_vfs, const vector<PgenSampleShard> &shards, uint32_t raw_variant_ct,
	          bool has_subset, const vector<uint32_t> &subset_sorted);

	//! Number of output (selected) samples across all shards.
	uint32_t OutputSampleCt() const {
		return output_sample_ct;
	}

	//! PgrGet over the combined space: `genovec` receives OutputSampleCt() 2-bit calls.
	plink2::PglErr Get(uint32_t vidx, uintptr_t *genovec);
	//! PgrGetCounts over the combined space: the per-shard counts are SUMMED.
	plink2::PglErr GetCounts(uint32_t vidx, STD_ARRAY_REF(uint32_t, 4) genocounts);
	//! PgrGetD over the combined space (dosage_main entries appended in sample order).
	plink2::PglErr GetD(uint32_t vidx, uintptr_t *genovec, uintptr_t *dosage_present, uint16_t *dosage_main,
	                    uint32_t *dosage_ct);
	//! PgrGetP over the combined space. A shard without phase contributes zero bits.
	plink2::PglErr GetP(uint32_t vidx, uintptr_t *genovec, uintptr_t *phasepresent, uintptr_t *phaseinfo,
	                    uint32_t *phasepresent_ct);

	//! Combined carrier list for the sparse (difflist) sample-count path (P13-011 §4.2).
	//! Returns, via `out_ids` / `out_geno` (both reader-owned, valid until the next call)
	//! and `out_len`, one (combined_sample_idx, genotype) pair for every non-hom_ref sample
	//! across all shards. Per shard it runs PgrGetDifflistOrGenovec with that shard's own
	//! SampleInclude()/pssi/out_ct and a per-shard threshold max(1, out_ct/8): a hom_ref-major
	//! difflist contributes only its carriers; a dense genovec or an ALT/missing-major difflist
	//! is walked densely (the latter re-decoded with the shard's own reader). hom_ref (g==0)
	//! never enters the list. `stats` (optional) tallies how each shard decoded.
	plink2::PglErr GetDifflistOrGenovec(uint32_t vidx, const uint32_t **out_ids, const uint8_t **out_geno,
	                                    uint32_t *out_len, ConcatSparseStats *stats);

private:
	struct Member {
		shared_ptr<PgenHeader> header; //! decoded .pgen header (cache tier 1); pgr points into it
		plink2::PgenReader pgr;
		AlignedBuffer pgr_alloc;
		bool open = false;

		uint32_t raw_sample_ct = 0;
		uint32_t out_ct = 0;     //! selected samples in this shard (0 = shard skipped)
		uint32_t out_offset = 0; //! first output position of this shard's samples
		bool has_subset = false; //! a strict subset of this shard is selected
		SampleSubset subset;     //! shard-LOCAL indices (valid when has_subset)
		plink2::PgrSampleSubsetIndex pssi;

		// Scratch (per shard, sized to out_ct; allocated on first use).
		AlignedBuffer genovec;
		AlignedBuffer bits_a; // dosage_present / phasepresent
		AlignedBuffer bits_b; // phaseinfo
		AlignedBuffer dosage_main;
		AlignedBuffer raregeno;     // sparse difflist: packed 2-bit genotypes (P13-011)
		AlignedBuffer difflist_ids; // sparse difflist: shard-local sample ids (uint32)

		const uintptr_t *SampleInclude() const {
			return has_subset ? subset.SampleInclude() : nullptr;
		}
		const uintptr_t *InterleavedVec() const {
			return has_subset ? subset.InterleavedVec() : nullptr;
		}
	};

	void EnsureScratch(Member &m, bool need_bits_a, bool need_bits_b, bool need_dosage_main);
	//! Allocate the sparse-path scratch (genovec + raregeno + difflist_ids) for a shard.
	void EnsureSparseScratch(Member &m);

	vector<unique_ptr<Member>> members; // only shards with out_ct > 0
	uint32_t output_sample_ct = 0;

	// Combined carrier output for GetDifflistOrGenovec, sized to OutputSampleCt(): a dense
	// or ALT-major shard can append ~out_ct entries, so the list can reach the full osc
	// (P13-011 §5.2 — never reuse the single-file smp_* buffers sized to osc/8).
	vector<uint32_t> combined_ids;
	vector<uint8_t> combined_geno;
};

//! OR the first `nbits` bits of `src` into `dst` starting at bit `dst_bit_offset`.
//! The destination bits at [dst_bit_offset, dst_bit_offset + nbits) must be zero.
//! Never touches a dst word past the last one the result occupies. Exposed for the
//! 2-bit genovec (nbits = 2 * sample_ct) and 1-bit arrays alike.
void OrBitsAtOffset(uintptr_t *dst, uint64_t dst_bit_offset, const uintptr_t *src, uint64_t nbits);

} // namespace duckdb
