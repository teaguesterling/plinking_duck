#include "pgen_concat_reader.hpp"
#include "pgen_vfs_opener.hpp"

#include <cstring>

namespace duckdb {

void OrBitsAtOffset(uintptr_t *dst, uint64_t dst_bit_offset, const uintptr_t *src, uint64_t nbits) {
	constexpr uint32_t W = plink2::kBitsPerWord;
	if (nbits == 0) {
		return;
	}
	const uint64_t word_idx = dst_bit_offset / W;
	const uint32_t shift = static_cast<uint32_t>(dst_bit_offset % W);
	const uint64_t src_word_ct = (nbits + W - 1) / W;
	const uint32_t tail_bits = static_cast<uint32_t>(nbits % W);
	for (uint64_t i = 0; i < src_word_ct; i++) {
		uintptr_t w = src[i];
		if (i + 1 == src_word_ct && tail_bits != 0) {
			// pgenlib does not promise zeroed trailing bits; never let them leak
			// into the next shard's slice.
			w &= (static_cast<uintptr_t>(1) << tail_bits) - 1;
		}
		dst[word_idx + i] |= w << shift;
		if (shift != 0) {
			uintptr_t hi = w >> (W - shift);
			// A non-zero carry lies inside the result, so its word exists; a zero
			// carry may point one word past the destination buffer — skip it.
			if (hi != 0) {
				dst[word_idx + i + 1] |= hi;
			}
		}
	}
}

ConcatPgenReader::~ConcatPgenReader() {
	for (auto &m : members) {
		if (m->open) {
			plink2::PglErr ce = plink2::kPglRetSuccess;
			plink2::CleanupPgr(&m->pgr, &ce);
			plink2::CleanupPgfi(&m->pgfi, &ce);
		}
	}
}

void ConcatPgenReader::Open(ClientContext &context, bool use_vfs, const vector<PgenSampleShard> &shards,
                            uint32_t raw_variant_ct, bool has_subset, const vector<uint32_t> &subset_sorted) {
	PgenVfsScope pgen_vfs_scope(context, use_vfs);
	members.clear();
	output_sample_ct = 0;
	idx_t subset_pos = 0; // subset_sorted is ascending; shards are in offset order
	for (auto &shard : shards) {
		auto m = make_uniq<Member>();
		m->raw_sample_ct = shard.raw_sample_ct;
		m->out_offset = output_sample_ct;
		vector<uint32_t> local;
		if (has_subset) {
			const uint32_t lo = shard.sample_offset;
			const uint32_t hi = shard.sample_offset + shard.raw_sample_ct;
			while (subset_pos < subset_sorted.size() && subset_sorted[subset_pos] < hi) {
				D_ASSERT(subset_sorted[subset_pos] >= lo);
				local.push_back(subset_sorted[subset_pos] - lo);
				subset_pos++;
			}
			m->out_ct = static_cast<uint32_t>(local.size());
		} else {
			m->out_ct = shard.raw_sample_ct;
		}
		if (m->out_ct == 0) {
			continue; // no selected sample lives here — never call pgenlib with sample_ct 0
		}
		m->has_subset = m->out_ct < shard.raw_sample_ct;
		if (m->has_subset) {
			m->subset = BuildSampleSubset(shard.raw_sample_ct, local);
		}

		plink2::PreinitPgfi(&m->pgfi);
		plink2::PreinitPgr(&m->pgr);
		char errstr_buf[plink2::kPglErrstrBufBlen];
		plink2::PgenHeaderCtrl header_ctrl;
		uintptr_t pgfi_alloc_cacheline_ct = 0;
		plink2::PglErr err =
		    plink2::PgfiInitPhase1(shard.pgen_path.c_str(), nullptr, raw_variant_ct, shard.raw_sample_ct, &header_ctrl,
		                           &m->pgfi, &pgfi_alloc_cacheline_ct, errstr_buf);
		if (err != plink2::kPglRetSuccess) {
			plink2::PglErr ce = plink2::kPglRetSuccess;
			plink2::CleanupPgfi(&m->pgfi, &ce);
			throw IOException("read_pfile: failed to open sample shard '%s': %s", shard.pgen_path, errstr_buf);
		}
		if (pgfi_alloc_cacheline_ct > 0) {
			m->pgfi_alloc.Allocate(pgfi_alloc_cacheline_ct * plink2::kCacheline);
		}
		uint32_t max_vrec_width = 0;
		uintptr_t pgr_alloc_cacheline_ct = 0;
		err = plink2::PgfiInitPhase2(header_ctrl, 0, 0, 0, 0, m->pgfi.raw_variant_ct, &max_vrec_width, &m->pgfi,
		                             m->pgfi_alloc.As<unsigned char>(), &pgr_alloc_cacheline_ct, errstr_buf);
		if (err != plink2::kPglRetSuccess) {
			plink2::PglErr ce = plink2::kPglRetSuccess;
			plink2::CleanupPgfi(&m->pgfi, &ce);
			throw IOException("read_pfile: failed to initialize sample shard '%s' (phase 2): %s", shard.pgen_path,
			                  errstr_buf);
		}
		if (pgr_alloc_cacheline_ct > 0) {
			m->pgr_alloc.Allocate(pgr_alloc_cacheline_ct * plink2::kCacheline);
		}
		err = plink2::PgrInit(shard.pgen_path.c_str(), max_vrec_width, &m->pgfi, &m->pgr,
		                      m->pgr_alloc.As<unsigned char>());
		if (err != plink2::kPglRetSuccess) {
			plink2::PglErr ce = plink2::kPglRetSuccess;
			plink2::CleanupPgr(&m->pgr, &ce);
			plink2::CleanupPgfi(&m->pgfi, &ce);
			throw IOException("read_pfile: PgrInit failed for sample shard '%s'", shard.pgen_path);
		}
		m->open = true;
		if (m->has_subset) {
			plink2::PgrSetSampleSubsetIndex(m->subset.CumulativePopcounts(), &m->pgr, &m->pssi);
		} else {
			plink2::PgrClearSampleSubsetIndex(&m->pgr, &m->pssi);
		}
		output_sample_ct += m->out_ct;
		members.push_back(std::move(m));
	}
}

void ConcatPgenReader::EnsureScratch(Member &m, bool need_bits_a, bool need_bits_b, bool need_dosage_main) {
	if (!m.genovec.ptr) {
		uintptr_t wc = plink2::NypCtToAlignedWordCt(m.out_ct);
		m.genovec.Allocate(wc * sizeof(uintptr_t));
		std::memset(m.genovec.ptr, 0, wc * sizeof(uintptr_t));
	}
	uintptr_t bit_wc = plink2::BitCtToAlignedWordCt(m.out_ct);
	if (need_bits_a && !m.bits_a.ptr) {
		m.bits_a.Allocate(bit_wc * sizeof(uintptr_t));
		std::memset(m.bits_a.ptr, 0, bit_wc * sizeof(uintptr_t));
	}
	if (need_bits_b && !m.bits_b.ptr) {
		m.bits_b.Allocate(bit_wc * sizeof(uintptr_t));
		std::memset(m.bits_b.ptr, 0, bit_wc * sizeof(uintptr_t));
	}
	if (need_dosage_main && !m.dosage_main.ptr) {
		m.dosage_main.Allocate(static_cast<idx_t>(m.out_ct) * sizeof(uint16_t));
	}
}

plink2::PglErr ConcatPgenReader::Get(uint32_t vidx, uintptr_t *genovec) {
	std::memset(genovec, 0, plink2::NypCtToWordCt(output_sample_ct) * sizeof(uintptr_t));
	for (auto &mp : members) {
		auto &m = *mp;
		EnsureScratch(m, false, false, false);
		auto *g = m.genovec.As<uintptr_t>();
		plink2::PglErr err = plink2::PgrGet(m.SampleInclude(), m.pssi, m.out_ct, vidx, &m.pgr, g);
		if (err != plink2::kPglRetSuccess) {
			return err;
		}
		OrBitsAtOffset(genovec, 2ULL * m.out_offset, g, 2ULL * m.out_ct);
	}
	return plink2::kPglRetSuccess;
}

plink2::PglErr ConcatPgenReader::GetCounts(uint32_t vidx, STD_ARRAY_REF(uint32_t, 4) genocounts) {
	genocounts[0] = genocounts[1] = genocounts[2] = genocounts[3] = 0;
	for (auto &mp : members) {
		auto &m = *mp;
		STD_ARRAY_DECL(uint32_t, 4, part);
		plink2::PglErr err =
		    plink2::PgrGetCounts(m.SampleInclude(), m.InterleavedVec(), m.pssi, m.out_ct, vidx, &m.pgr, part);
		if (err != plink2::kPglRetSuccess) {
			return err;
		}
		for (int i = 0; i < 4; i++) {
			genocounts[i] += part[i];
		}
	}
	return plink2::kPglRetSuccess;
}

plink2::PglErr ConcatPgenReader::GetD(uint32_t vidx, uintptr_t *genovec, uintptr_t *dosage_present,
                                      uint16_t *dosage_main, uint32_t *dosage_ct) {
	std::memset(genovec, 0, plink2::NypCtToWordCt(output_sample_ct) * sizeof(uintptr_t));
	std::memset(dosage_present, 0, plink2::BitCtToWordCt(output_sample_ct) * sizeof(uintptr_t));
	uint32_t total = 0;
	for (auto &mp : members) {
		auto &m = *mp;
		EnsureScratch(m, true, false, true);
		auto *g = m.genovec.As<uintptr_t>();
		uint32_t part_ct = 0;
		plink2::PglErr err = plink2::PgrGetD(m.SampleInclude(), m.pssi, m.out_ct, vidx, &m.pgr, g,
		                                     m.bits_a.As<uintptr_t>(), m.dosage_main.As<uint16_t>(), &part_ct);
		if (err != plink2::kPglRetSuccess) {
			return err;
		}
		OrBitsAtOffset(genovec, 2ULL * m.out_offset, g, 2ULL * m.out_ct);
		if (part_ct != 0) {
			// dosage_main holds one entry per SET dosage_present bit, in sample order, so
			// shard slices append in shard order.
			OrBitsAtOffset(dosage_present, m.out_offset, m.bits_a.As<uintptr_t>(), m.out_ct);
			std::memcpy(dosage_main + total, m.dosage_main.As<uint16_t>(), part_ct * sizeof(uint16_t));
			total += part_ct;
		}
	}
	*dosage_ct = total;
	return plink2::kPglRetSuccess;
}

plink2::PglErr ConcatPgenReader::GetP(uint32_t vidx, uintptr_t *genovec, uintptr_t *phasepresent, uintptr_t *phaseinfo,
                                      uint32_t *phasepresent_ct) {
	const uintptr_t bit_bytes = plink2::BitCtToWordCt(output_sample_ct) * sizeof(uintptr_t);
	std::memset(genovec, 0, plink2::NypCtToWordCt(output_sample_ct) * sizeof(uintptr_t));
	std::memset(phasepresent, 0, bit_bytes);
	std::memset(phaseinfo, 0, bit_bytes);
	uint32_t total = 0;
	for (auto &mp : members) {
		auto &m = *mp;
		EnsureScratch(m, true, true, false);
		auto *g = m.genovec.As<uintptr_t>();
		uint32_t part_ct = 0;
		plink2::PglErr err = plink2::PgrGetP(m.SampleInclude(), m.pssi, m.out_ct, vidx, &m.pgr, g,
		                                     m.bits_a.As<uintptr_t>(), m.bits_b.As<uintptr_t>(), &part_ct);
		if (err != plink2::kPglRetSuccess) {
			return err;
		}
		OrBitsAtOffset(genovec, 2ULL * m.out_offset, g, 2ULL * m.out_ct);
		// pgenlib leaves the phase buffers untouched for a variant this shard has no
		// phase for (part_ct == 0); its slice then stays zero (unphased).
		if (part_ct != 0) {
			OrBitsAtOffset(phasepresent, m.out_offset, m.bits_a.As<uintptr_t>(), m.out_ct);
			OrBitsAtOffset(phaseinfo, m.out_offset, m.bits_b.As<uintptr_t>(), m.out_ct);
			total += part_ct;
		}
	}
	*phasepresent_ct = total;
	return plink2::kPglRetSuccess;
}

} // namespace duckdb
