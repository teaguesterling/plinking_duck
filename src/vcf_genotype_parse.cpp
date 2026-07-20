// vcf_genotype_parse.cpp — VCF GT field parsing for biallelic variants
//
// ============================ VENDORED-CODE PIN =============================
// Upstream file:      third_party/plink-ng/2.0/plink2_import.cc
// plink-ng submodule: pinned commit 4ce97faa08bc370bedb30dcc82b4eeeef1c7c1f4
//                     (verify: git -C third_party/plink-ng rev-parse HEAD)
// Extracted pieces + upstream line ranges (at the pinned commit):
//   VcfImportBaseContext struct           L869-881
//   VcfCheckQuals                         L898-920
//   VcfParseErr enum                      L1109-1116
//   VcfConvertUnphasedBiallelicLine       L1480-1579
//   VcfConvertPhasedBiallelicLine         L1962-2091
// Extract type:       function BODIES are identical to upstream MODULO two documented
//   renamings applied here (a local enum in src/include/vcf_genotype_parse.hpp is used):
//     VcfHalfCall (type)  -> uint32_t
//     kVcfHalfCall* (enum constants) -> kHalfCall*
//   Reformatted by clang-format; wrapped in namespace duckdb. The other functions in
//   this file (FillBaseContext, TranslateResult, ParseUnphasedBiallelicGT,
//   ParsePhasedBiallelicGT) are PlinkingDuck glue, NOT upstream copies.
// Why extracted (not linkable): plink2_import.cc is an 18k-line CLI/import-pipeline TU
//   and cannot be linked wholesale; the copied functions depend only on header-inline
//   utilities from plink2_base.h / plink2_string.h.
// Drift canary:       scripts/check_vendored_drift.sh (run it in CI / before release;
//                     it compares these function bodies to the pinned upstream, applying
//                     the two renamings above so the comparison stays an exact match).
// Re-sync procedure:  when the plink-ng submodule pin is bumped and the canary fails,
//   re-copy the named functions from the ranges above (adjusting for line shifts),
//   reapply the VcfHalfCall/kVcfHalfCall renamings, then update this pin block's commit
//   hash and line numbers.
// Copyright (C) 2005-2026 Shaun Purcell, Christopher Chang.  Licensed under GPL v3+.
// ===========================================================================

#include "vcf_genotype_parse.hpp"

#include <plink2_base.h>
#include <plink2_string.h>

// All plink2 internal types/functions are in the plink2 namespace (or global in C mode).
// We use them unqualified within the anonymous namespace below.
using namespace plink2;

namespace duckdb {

namespace {

// --------------------------------------------------------------------------
// Copied from plink2_import.cc — VcfImportBaseContext struct (lines 869-881)
// --------------------------------------------------------------------------
typedef struct VcfImportContextBaseStruct {
	uint32_t sample_ct;
	uint32_t halfcall_mode; // VcfHalfCall values (uint32_t enum)
	uint32_t error_on_polyploid;
	uint32_t gt_exists;
	STD_ARRAY_DECL(uint32_t, 2, qual_field_skips);
	STD_ARRAY_DECL(int32_t, 2, qual_line_mins);
	STD_ARRAY_DECL(int32_t, 2, qual_line_maxs);
	uint32_t qual_field_ct;
} VcfImportBaseContext;

// --------------------------------------------------------------------------
// Copied from plink2_import.cc — VcfParseErr enum (lines 1109-1116)
// --------------------------------------------------------------------------
ENUM_U31_DEF_START()
kVcfParseOk, kVcfParseMissingTokens, kVcfParseInvalidGt, kVcfParseHalfCallError, kVcfParseInvalidDosage,
    kVcfParsePolyploidError ENUM_U31_DEF_END(VcfParseErr);

// --------------------------------------------------------------------------
// Copied from plink2_import.cc — VcfCheckQuals (lines 898-920)
// --------------------------------------------------------------------------
// returns 1 if a quality check failed
// assumes either 1 or 2 qual fields, otherwise change this to a loop
uint32_t VcfCheckQuals(STD_ARRAY_KREF(uint32_t, 2) qual_field_skips, STD_ARRAY_KREF(int32_t, 2) qual_line_mins,
                       STD_ARRAY_KREF(int32_t, 2) qual_line_maxs, const char *gtext_iter, const char *gtext_end,
                       uint32_t qual_field_ct) {
	const uint32_t skip0 = qual_field_skips[0];
	if (skip0) {
		gtext_iter = AdvToNthDelimChecked(gtext_iter, gtext_end, skip0, ':');
		if (!gtext_iter) {
			return 0;
		}
		++gtext_iter;
	}
	int32_t ii;
	if ((!ScanInt32(gtext_iter, &ii)) && ((ii < qual_line_mins[0]) || (ii > qual_line_maxs[0]))) {
		return 1;
	}
	if (qual_field_ct == 1) {
		return 0;
	}
	gtext_iter = AdvToNthDelimChecked(gtext_iter, gtext_end, qual_field_skips[1], ':');
	if (!gtext_iter) {
		return 0;
	}
	++gtext_iter;
	return (!ScanInt32(gtext_iter, &ii)) && ((ii < qual_line_mins[1]) || (ii > qual_line_maxs[1]));
}

// --------------------------------------------------------------------------
// Copied from plink2_import.cc — VcfConvertUnphasedBiallelicLine (lines 1480-1578)
// --------------------------------------------------------------------------
VcfParseErr VcfConvertUnphasedBiallelicLine(const VcfImportBaseContext *vibcp, const char *linebuf_iter,
                                            uintptr_t *genovec) {
	const uint32_t sample_ct = vibcp->sample_ct;
	const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
	const uint32_t halfcall_mode = vibcp->halfcall_mode;
	const uint32_t error_on_polyploid = vibcp->error_on_polyploid;
	STD_ARRAY_KREF(uint32_t, 2) qual_field_skips = vibcp->qual_field_skips;
	STD_ARRAY_KREF(int32_t, 2) qual_line_mins = vibcp->qual_line_mins;
	STD_ARRAY_KREF(int32_t, 2) qual_line_maxs = vibcp->qual_line_maxs;
	const uint32_t qual_field_ct = vibcp->qual_field_ct;

	uint32_t inner_loop_last = kBitsPerWordD2 - 1;
	for (uint32_t widx = 0;; ++widx) {
		if (widx >= sample_ctl2_m1) {
			if (widx > sample_ctl2_m1) {
				break;
			}
			inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
		}
		uintptr_t genovec_word = 0;
		for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
			const char *cur_gtext_end = FirstPrespace(linebuf_iter);
			if (unlikely((*cur_gtext_end != '\t') &&
			             ((sample_idx_lowbits != inner_loop_last) || (widx != sample_ctl2_m1)))) {
				return kVcfParseMissingTokens;
			}
			uintptr_t cur_geno;
			if (qual_field_ct && VcfCheckQuals(qual_field_skips, qual_line_mins, qual_line_maxs, linebuf_iter,
			                                   cur_gtext_end, qual_field_ct)) {
				cur_geno = 3;
			} else {
				const uint32_t is_haploid = (linebuf_iter[1] != '/') && (linebuf_iter[1] != '|');
				cur_geno = ctow(*linebuf_iter) - 48;
				if (cur_geno <= 1) {
					if (is_haploid) {
						cur_geno *= 2;
					} else {
						const char polyploid_char = linebuf_iter[3];
						if ((polyploid_char == '/') || (polyploid_char == '|')) {
							if (unlikely(error_on_polyploid)) {
								return kVcfParsePolyploidError;
							}
							cur_geno = 3;
						} else {
							const uintptr_t second_allele_idx = ctow(linebuf_iter[2]) - 48;
							if (second_allele_idx <= 1) {
								cur_geno += second_allele_idx;
							} else if (unlikely(second_allele_idx != (~k0LU) * 2)) {
								return kVcfParseInvalidGt;
							} else if (halfcall_mode == kHalfCallMissing) {
								cur_geno = 3;
							} else if (unlikely(halfcall_mode == kHalfCallError)) {
								return kVcfParseHalfCallError;
							} else {
								cur_geno <<= halfcall_mode;
							}
						}
					}
				} else {
					if (unlikely(cur_geno != (~k0LU) * 2)) {
						return kVcfParseInvalidGt;
					}
					cur_geno = 3;
					if (!is_haploid) {
						const char polyploid_char = linebuf_iter[3];
						if ((polyploid_char == '/') || (polyploid_char == '|')) {
							if (unlikely(error_on_polyploid)) {
								return kVcfParsePolyploidError;
							}
						} else {
							const char second_allele_char = linebuf_iter[2];
							if ((second_allele_char != '.') && (halfcall_mode != kHalfCallMissing)) {
								cur_geno = ctow(second_allele_char) - 48;
								if (unlikely(cur_geno > 1)) {
									return kVcfParseInvalidGt;
								}
								if (unlikely(halfcall_mode == kHalfCallError)) {
									return kVcfParseHalfCallError;
								}
								cur_geno <<= halfcall_mode;
							}
						}
					}
				}
			}
			genovec_word |= cur_geno << (2 * sample_idx_lowbits);
			linebuf_iter = &(cur_gtext_end[1]);
		}
		genovec[widx] = genovec_word;
	}
	return kVcfParseOk;
}

// --------------------------------------------------------------------------
// Copied from plink2_import.cc — VcfConvertPhasedBiallelicLine (lines 1962-2091)
// --------------------------------------------------------------------------
VcfParseErr VcfConvertPhasedBiallelicLine(const VcfImportBaseContext *vibcp, const char *linebuf_iter,
                                          uintptr_t *genovec, uintptr_t *phasepresent, uintptr_t *phaseinfo) {
	const uint32_t sample_ct = vibcp->sample_ct;
	const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
	const uint32_t halfcall_mode = vibcp->halfcall_mode;
	const uint32_t error_on_polyploid = vibcp->error_on_polyploid;
	STD_ARRAY_KREF(uint32_t, 2) qual_field_skips = vibcp->qual_field_skips;
	STD_ARRAY_KREF(int32_t, 2) qual_line_mins = vibcp->qual_line_mins;
	STD_ARRAY_KREF(int32_t, 2) qual_line_maxs = vibcp->qual_line_maxs;
	const uint32_t qual_field_ct = vibcp->qual_field_ct;
	Halfword *phasepresent_alias = R_CAST(Halfword *, phasepresent);
	Halfword *phaseinfo_alias = R_CAST(Halfword *, phaseinfo);

	uint32_t inner_loop_last = kBitsPerWordD2 - 1;
	for (uint32_t widx = 0;; ++widx) {
		if (widx >= sample_ctl2_m1) {
			if (widx > sample_ctl2_m1) {
				if (widx % 2) {
					phasepresent_alias[widx] = 0;
				}
				break;
			}
			inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
		}
		uintptr_t genovec_word = 0;
		uint32_t phasepresent_hw = 0;
		uint32_t phaseinfo_hw = 0;
		for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
			const char *cur_gtext_end = FirstPrespace(linebuf_iter);
			if (unlikely((*cur_gtext_end != '\t') &&
			             ((sample_idx_lowbits != inner_loop_last) || (widx != sample_ctl2_m1)))) {
				return kVcfParseMissingTokens;
			}
			uintptr_t cur_geno;
			if (qual_field_ct && VcfCheckQuals(qual_field_skips, qual_line_mins, qual_line_maxs, linebuf_iter,
			                                   cur_gtext_end, qual_field_ct)) {
				cur_geno = 3;
			} else {
				const uint32_t is_phased = (linebuf_iter[1] == '|');
				const uint32_t is_haploid = (!is_phased) && (linebuf_iter[1] != '/');
				cur_geno = ctow(*linebuf_iter) - 48;
				if (cur_geno <= 1) {
					if (is_haploid) {
						cur_geno *= 2;
					} else {
						const char polyploid_char = linebuf_iter[3];
						if ((polyploid_char == '/') || (polyploid_char == '|')) {
							if (unlikely(error_on_polyploid)) {
								return kVcfParsePolyploidError;
							}
							cur_geno = 3;
						} else {
							const uintptr_t second_allele_idx = ctow(linebuf_iter[2]) - 48;
							if (second_allele_idx <= 1) {
								cur_geno += second_allele_idx;
								if (is_phased && (cur_geno == 1)) {
									const uint32_t shifted_bit = 1U << sample_idx_lowbits;
									phasepresent_hw |= shifted_bit;
#ifdef USE_AVX2
									if (!second_allele_idx) {
										phaseinfo_hw |= shifted_bit;
									}
#else
									phaseinfo_hw |= shifted_bit & (second_allele_idx - 1);
#endif
								}
							} else if (unlikely(second_allele_idx != (~k0LU) * 2)) {
								return kVcfParseInvalidGt;
							} else if (halfcall_mode == kHalfCallMissing) {
								cur_geno = 3;
							} else if (unlikely(halfcall_mode == kHalfCallError)) {
								return kVcfParseHalfCallError;
							} else {
								if (is_phased && (halfcall_mode == kHalfCallReference) && (cur_geno == 1)) {
									const uint32_t shifted_bit = 1U << sample_idx_lowbits;
									phasepresent_hw |= shifted_bit;
									phaseinfo_hw |= shifted_bit;
								} else {
									cur_geno <<= halfcall_mode;
								}
							}
						}
					}
				} else if (unlikely(cur_geno != (~k0LU) * 2)) {
					return kVcfParseInvalidGt;
				} else {
					cur_geno = 3;
					if (!is_haploid) {
						const char second_allele_char = linebuf_iter[2];
						const char polyploid_char = linebuf_iter[3];
						if ((polyploid_char == '/') || (polyploid_char == '|')) {
							if (unlikely(error_on_polyploid)) {
								return kVcfParsePolyploidError;
							}
						} else if ((second_allele_char != '.') && (halfcall_mode != kHalfCallMissing)) {
							cur_geno = ctow(second_allele_char) - 48;
							if (unlikely(cur_geno > 1)) {
								return kVcfParseInvalidGt;
							}
							if (unlikely(halfcall_mode == kHalfCallError)) {
								return kVcfParseHalfCallError;
							}
							if (is_phased && (halfcall_mode == kHalfCallReference) && (cur_geno == 1)) {
								const uint32_t shifted_bit = 1U << sample_idx_lowbits;
								phasepresent_hw |= shifted_bit;
							} else {
								cur_geno <<= halfcall_mode;
							}
						}
					}
				}
			}
			genovec_word |= cur_geno << (2 * sample_idx_lowbits);
			linebuf_iter = &(cur_gtext_end[1]);
		}
		genovec[widx] = genovec_word;
		phasepresent_alias[widx] = phasepresent_hw;
		phaseinfo_alias[widx] = phaseinfo_hw;
	}
	return kVcfParseOk;
}

// --------------------------------------------------------------------------
// Helper: populate VcfImportBaseContext from our public VcfParseContext
// --------------------------------------------------------------------------
void FillBaseContext(const VcfParseContext &ctx, VcfImportBaseContext *vibcp) {
	vibcp->sample_ct = ctx.sample_ct;
	vibcp->halfcall_mode = static_cast<uint32_t>(ctx.halfcall_mode);
	vibcp->error_on_polyploid = 1; // always error on polyploid for our use
	vibcp->gt_exists = 1;
	vibcp->qual_field_skips[0] = ctx.qual_field_skips[0];
	vibcp->qual_field_skips[1] = ctx.qual_field_skips[1];
	vibcp->qual_line_mins[0] = ctx.qual_line_mins[0];
	vibcp->qual_line_mins[1] = ctx.qual_line_mins[1];
	vibcp->qual_line_maxs[0] = ctx.qual_line_maxs[0];
	vibcp->qual_line_maxs[1] = ctx.qual_line_maxs[1];
	vibcp->qual_field_ct = ctx.qual_field_ct;
}

// --------------------------------------------------------------------------
// Convert VcfParseErr to VcfGenoParseResult
// --------------------------------------------------------------------------
VcfGenoParseResult TranslateResult(VcfParseErr err) {
	switch (err) {
	case kVcfParseOk:
		return VcfGenoParseResult::OK;
	case kVcfParseMissingTokens:
		return VcfGenoParseResult::MISSING_TOKENS;
	case kVcfParseInvalidGt:
		return VcfGenoParseResult::INVALID_GT;
	case kVcfParseHalfCallError:
		return VcfGenoParseResult::HALFCALL_ERROR;
	case kVcfParsePolyploidError:
		return VcfGenoParseResult::POLYPLOID_ERROR;
	default:
		return VcfGenoParseResult::INVALID_GT;
	}
}

} // anonymous namespace

// ==========================================================================
// Public API
// ==========================================================================

VcfGenoParseResult ParseUnphasedBiallelicGT(const VcfParseContext &ctx, const char *linebuf_iter, uintptr_t *genovec) {
	VcfImportBaseContext vibc;
	FillBaseContext(ctx, &vibc);
	VcfParseErr err = VcfConvertUnphasedBiallelicLine(&vibc, linebuf_iter, genovec);
	return TranslateResult(err);
}

VcfGenoParseResult ParsePhasedBiallelicGT(const VcfParseContext &ctx, const char *linebuf_iter, uintptr_t *genovec,
                                          uintptr_t *phasepresent, uintptr_t *phaseinfo) {
	VcfImportBaseContext vibc;
	FillBaseContext(ctx, &vibc);
	VcfParseErr err = VcfConvertPhasedBiallelicLine(&vibc, linebuf_iter, genovec, phasepresent, phaseinfo);
	return TranslateResult(err);
}

} // namespace duckdb
