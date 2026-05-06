// vcf_genotype_parse.hpp — VCF GT field parsing for biallelic variants
//
// Extracted from plink-ng/2.0/plink2_import.cc (GPLv3)
// Copyright (C) 2005-2026 Shaun Purcell, Christopher Chang
// Adapted for DuckDB extension use.

#pragma once

#include <plink2_base.h>

#include <cstdint>

namespace duckdb {

// Half-call handling modes (values match plink2::VcfHalfCall in plink2_import.h)
enum VcfHalfCallMode : uint32_t {
	kHalfCallReference = 0,
	kHalfCallHaploid = 1,
	kHalfCallMissing = 2,
	kHalfCallError = 3
};

// Parse result codes
enum class VcfGenoParseResult : uint8_t { OK, MISSING_TOKENS, INVALID_GT, HALFCALL_ERROR, POLYPLOID_ERROR };

// Context for VCF GT parsing — maps to VcfImportBaseContext fields
struct VcfParseContext {
	uint32_t sample_ct;
	VcfHalfCallMode halfcall_mode;
	uint32_t qual_field_ct;       // 0 = no quality filtering
	uint32_t qual_field_skips[2]; // colon-delimited skips to reach each quality field
	int32_t qual_line_mins[2];
	int32_t qual_line_maxs[2];
};

// Parse unphased biallelic GT fields from a VCF data line.
// linebuf_iter should point to the first sample's GT field (after FORMAT tab).
// genovec must be allocated with NypCtToAlignedWordCt(sample_ct) words.
// Output: 2-bit packed genotype vector (0=hom_ref, 1=het, 2=hom_alt, 3=missing).
VcfGenoParseResult ParseUnphasedBiallelicGT(const VcfParseContext &ctx, const char *linebuf_iter, uintptr_t *genovec);

// Parse phased biallelic GT fields from a VCF data line.
// linebuf_iter should point to the first sample's GT field (after FORMAT tab).
// genovec must be allocated with NypCtToAlignedWordCt(sample_ct) words.
// phasepresent must be allocated with BitCtToAlignedWordCt(sample_ct) words (zeroed).
// phaseinfo must be allocated with BitCtToAlignedWordCt(sample_ct) words (zeroed).
// Output:
//   genovec: 2-bit packed genotype (same as unphased)
//   phasepresent: bit set for each het sample with explicit phase
//   phaseinfo: bit set for each het where first allele is ALT (1|0)
VcfGenoParseResult ParsePhasedBiallelicGT(const VcfParseContext &ctx, const char *linebuf_iter, uintptr_t *genovec,
                                          uintptr_t *phasepresent, uintptr_t *phaseinfo);

} // namespace duckdb
