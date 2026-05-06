# read_plink_vcf() Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `read_plink_vcf()` table function that reads VCF files using plink-ng's optimized GT parsing helpers and outputs genotypes in the same format as `read_pfile()`.

**Architecture:** DuckDB FileSystem handles file I/O (plain/gzip/bgzf). Header parsing extracts sample names. Per-line, fixed VCF columns are split, then plink-ng's `VcfConvertUnphasedBiallelicLine` / `VcfConvertPhasedBiallelicLine` parse GT fields into genovec. Shared output utilities (extracted from pfile_reader) format genovec into ARRAY/LIST/columns output vectors.

**Tech Stack:** C++17, DuckDB extension API, plink-ng header-inline utilities (plink2_base.h, plink2_string.h), existing plink_common infrastructure.

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `src/vcf_genotype_parse.cpp` | Create | Extracted plink-ng VCF GT parsing functions with wrapper API |
| `src/include/vcf_genotype_parse.hpp` | Create | Public interface for VCF genotype parsing |
| `src/vcf_reader.cpp` | Create | Table function: bind/init/scan for VCF files |
| `src/include/vcf_reader.hpp` | Create | Registration function declaration |
| `src/include/plink_common.hpp` | Modify | Add `FillGenotypeVector` shared output utility declaration |
| `src/plink_common.cpp` | Modify | Implement `FillGenotypeVector` |
| `src/pfile_reader.cpp` | Modify | Replace inline output logic with `FillGenotypeVector` calls |
| `src/plinking_duck_extension.cpp` | Modify | Register `read_plink_vcf` |
| `CMakeLists.txt` | Modify | Add new source files |
| `test/data/vcf_example.vcf` | Create | Test VCF with GQ/DP fields |
| `test/data/vcf_phased.vcf` | Create | Test VCF with phased genotypes |
| `test/sql/read_plink_vcf.test` | Create | Positive tests |
| `test/sql/read_plink_vcf_negative.test` | Create | Error case tests |

---

### Task 1: Create VCF Test Fixtures

**Files:**
- Create: `test/data/vcf_example.vcf`
- Create: `test/data/vcf_phased.vcf`

These fixtures have the same genotype values as `test/data/pfile_example` (4 variants × 4 samples) to enable cross-validation. The existing `test/data/example.vcf` already matches, but we need variants with GQ/DP fields and phased genotypes.

- [ ] **Step 1: Create `test/data/vcf_example.vcf` with quality fields**

```
##fileformat=VCFv4.3
##contig=<ID=1,length=100000>
##contig=<ID=2,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3	SAMPLE4
1	10000	rs1	A	G	.	.	.	GT:GQ:DP	0/0:30:10	0/1:25:8	1/1:35:15	./.:0:0
1	20000	rs2	C	T	.	.	.	GT:GQ:DP	0/1:40:12	0/1:5:3	0/0:30:10	1/1:28:9
1	30000	rs3	G	A	.	.	.	GT:GQ:DP	1/1:45:20	./.:0:0	0/1:20:7	0/0:35:11
2	15000	rs4	T	C	.	.	.	GT:GQ:DP	0/0:50:25	0/0:30:10	0/1:15:5	1/1:40:18
```

Genotype values (matching pfile_example):
- rs1: [0, 1, 2, NULL]
- rs2: [1, 1, 0, 2]
- rs3: [2, NULL, 1, 0]
- rs4: [0, 0, 1, 2]

GQ values designed so min_gq=20 filters: rs2/SAMPLE2 (GQ=5), rs4/SAMPLE3 (GQ=15)
DP values designed so min_dp=8 filters: rs2/SAMPLE2 (DP=3), rs3/SAMPLE3 (DP=7), rs4/SAMPLE3 (DP=5)

- [ ] **Step 2: Create `test/data/vcf_phased.vcf` with phased genotypes**

```
##fileformat=VCFv4.3
##contig=<ID=1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3	SAMPLE4
1	10000	rs1	A	G	.	.	.	GT	0|0	0|1	1|0	1|1
1	20000	rs2	C	T	.	.	.	GT	0|1	1|0	0|0	./.
1	30000	rs3	G	A	.	.	.	GT	1|1	0|1	1|0	0|0
```

Phase encoding (matches existing phased_example):
- `0|1` → alleles [0, 1] (REF|ALT)
- `1|0` → alleles [1, 0] (ALT|REF)
- `0|0` → alleles [0, 0]
- `1|1` → alleles [1, 1]
- `./.` → NULL

- [ ] **Step 3: Commit test fixtures**

```bash
git add test/data/vcf_example.vcf test/data/vcf_phased.vcf
git commit -m "test: add VCF fixtures for read_plink_vcf"
```

---

### Task 2: Extract plink-ng VCF GT Parsing Helpers

**Files:**
- Create: `src/vcf_genotype_parse.cpp`
- Create: `src/include/vcf_genotype_parse.hpp`

Extract `VcfConvertUnphasedBiallelicLine`, `VcfConvertPhasedBiallelicLine`, and `VcfCheckQuals` from `third_party/plink-ng/2.0/plink2_import.cc` (lines 898-920, 1480-1578, 1962-2091). These functions depend only on header-inline utilities already compiled via pgenlib.

- [ ] **Step 1: Create `src/include/vcf_genotype_parse.hpp`**

```cpp
#pragma once

#include <pgenlib_misc.h>
#include <plink2_base.h>
#include <plink2_string.h>

#include <cstdint>

namespace duckdb {

// Re-export plink2 VcfHalfCall enum for use in our API
using plink2::VcfHalfCall;
using plink2::kVcfHalfCallReference;
using plink2::kVcfHalfCallHaploid;
using plink2::kVcfHalfCallMissing;
using plink2::kVcfHalfCallError;

enum class VcfGenoParseResult : uint8_t {
    OK,
    MISSING_TOKENS,
    INVALID_GT,
    HALFCALL_ERROR,
    POLYPLOID_ERROR
};

struct VcfParseContext {
    uint32_t sample_ct;
    VcfHalfCall halfcall_mode;
    // Quality filtering: field positions within FORMAT (0-based index after GT)
    // Set qual_field_ct = 0 to disable quality filtering
    uint32_t qual_field_ct;
    // qual_field_skips[i] = number of ':' delimiters to skip from current position to reach qual field i
    uint32_t qual_field_skips[2];
    int32_t qual_line_mins[2];   // minimum acceptable value (-1 = no min)
    int32_t qual_line_maxs[2];   // maximum acceptable value (INT32_MAX = no max)
};

// Parse unphased biallelic GT fields from sample data portion of VCF line.
// linebuf_iter points to first character of first sample's data (after the last FORMAT\t).
// genovec must be allocated via NypCtToAlignedWordCt(sample_ct).
VcfGenoParseResult ParseUnphasedBiallelicGT(
    const VcfParseContext &ctx,
    const char *linebuf_iter,
    uintptr_t *genovec);

// Parse phased biallelic GT fields from sample data portion of VCF line.
// phasepresent and phaseinfo must be allocated to BitCtToAlignedWordCt(sample_ct).
VcfGenoParseResult ParsePhasedBiallelicGT(
    const VcfParseContext &ctx,
    const char *linebuf_iter,
    uintptr_t *genovec,
    uintptr_t *phasepresent,
    uintptr_t *phaseinfo);

} // namespace duckdb
```

- [ ] **Step 2: Create `src/vcf_genotype_parse.cpp`**

This file adapts the plink-ng functions to our wrapper interface. The core logic is copied from `plink2_import.cc` (GPLv3, attributed). The functions use only header-inline helpers (`FirstPrespace`, `ctow`, `AdvToNthDelimChecked`, `ScanInt32` from plink2_base.h/plink2_string.h).

```cpp
// VCF genotype parsing helpers extracted from plink-ng (plink2_import.cc).
// Copyright (C) 2005-2026 Shaun Purcell, Christopher Chang.
// Licensed under GNU General Public License v3.
//
// These functions parse VCF GT fields for biallelic variants into plink2's
// 2-bit genovec format. Only the biallelic parsers are included; multiallelic
// support is planned for a future version.

#include "vcf_genotype_parse.hpp"

namespace duckdb {

// Internal: maps VcfParseContext to plink2's VcfImportBaseContext for the call
static plink2::VcfImportBaseContext MakeVibc(const VcfParseContext &ctx) {
    plink2::VcfImportBaseContext vibc;
    vibc.sample_ct = ctx.sample_ct;
    vibc.halfcall_mode = ctx.halfcall_mode;
    vibc.error_on_polyploid = 0;  // polyploid → missing (not error)
    vibc.gt_exists = 1;
    vibc.qual_field_ct = ctx.qual_field_ct;
    vibc.qual_field_skips[0] = ctx.qual_field_skips[0];
    vibc.qual_field_skips[1] = ctx.qual_field_skips[1];
    vibc.qual_line_mins[0] = ctx.qual_line_mins[0];
    vibc.qual_line_mins[1] = ctx.qual_line_mins[1];
    vibc.qual_line_maxs[0] = ctx.qual_line_maxs[0];
    vibc.qual_line_maxs[1] = ctx.qual_line_maxs[1];
    return vibc;
}

// --- Copied from plink2_import.cc (VcfCheckQuals, lines 898-920) ---
// [paste the ~25 lines of VcfCheckQuals here]

// --- Copied from plink2_import.cc (VcfConvertUnphasedBiallelicLine, lines 1480-1578) ---
// [paste the ~100 lines here]

// --- Copied from plink2_import.cc (VcfConvertPhasedBiallelicLine, lines 1962-2091) ---
// [paste the ~130 lines here]

// --- Public wrapper functions ---

static VcfGenoParseResult TranslateResult(plink2::VcfParseErr err) {
    switch (err) {
    case plink2::kVcfParseOk: return VcfGenoParseResult::OK;
    case plink2::kVcfParseMissingTokens: return VcfGenoParseResult::MISSING_TOKENS;
    case plink2::kVcfParseInvalidGt: return VcfGenoParseResult::INVALID_GT;
    case plink2::kVcfParseHalfCallError: return VcfGenoParseResult::HALFCALL_ERROR;
    case plink2::kVcfParsePolyploidError: return VcfGenoParseResult::POLYPLOID_ERROR;
    default: return VcfGenoParseResult::INVALID_GT;
    }
}

VcfGenoParseResult ParseUnphasedBiallelicGT(
    const VcfParseContext &ctx,
    const char *linebuf_iter,
    uintptr_t *genovec) {
    plink2::VcfImportBaseContext vibc = MakeVibc(ctx);
    auto err = VcfConvertUnphasedBiallelicLine(&vibc, linebuf_iter, genovec);
    return TranslateResult(err);
}

VcfGenoParseResult ParsePhasedBiallelicGT(
    const VcfParseContext &ctx,
    const char *linebuf_iter,
    uintptr_t *genovec,
    uintptr_t *phasepresent,
    uintptr_t *phaseinfo) {
    plink2::VcfImportBaseContext vibc = MakeVibc(ctx);
    auto err = VcfConvertPhasedBiallelicLine(&vibc, linebuf_iter, genovec, phasepresent, phaseinfo);
    return TranslateResult(err);
}

} // namespace duckdb
```

The actual step is: copy the three functions verbatim from `third_party/plink-ng/2.0/plink2_import.cc` into the file body, wrapped in an anonymous namespace. The `VcfParseErr` enum also needs to be copied (lines 1110-1116) since it's local to plink2_import.cc.

- [ ] **Step 3: Verify it compiles**

Add to CMakeLists.txt temporarily and build:

```bash
# In CMakeLists.txt EXTENSION_SOURCES, add:
#     src/vcf_genotype_parse.cpp
make -j4
```

Fix any missing includes or type issues. The likely issues:
- `STD_ARRAY_KREF` macro — defined in `plink2_base.h`, should be available
- `VcfParseErr` enum — needs to be copied locally since it's defined inside plink2_import.cc
- `Halfword` type — defined in `plink2_base.h`

- [ ] **Step 4: Commit**

```bash
git add src/vcf_genotype_parse.cpp src/include/vcf_genotype_parse.hpp CMakeLists.txt
git commit -m "feat: extract plink-ng VCF GT parsing helpers for read_plink_vcf"
```

---

### Task 3: Extract Shared Genotype Output Utility

**Files:**
- Modify: `src/include/plink_common.hpp`
- Modify: `src/plink_common.cpp`
- Modify: `src/pfile_reader.cpp` (lines ~1900-2043)

Extract the ARRAY/LIST genotype output logic from `pfile_reader.cpp` into a shared function in `plink_common`. This function writes genotype data into a DuckDB Vector for one row, handling all modes: unphased ARRAY, unphased LIST, phased ARRAY, phased LIST.

The COLUMNS mode is NOT extracted here — it's handled separately in the column-dispatch switch and stays inline.

- [ ] **Step 1: Add declaration to `src/include/plink_common.hpp`**

After the `UnpackPhasedGenotypes` declaration (line 457), add:

```cpp
// ---------------------------------------------------------------------------
// Shared genotype output vector filling
// ---------------------------------------------------------------------------

//! Fill a genotype output vector for one row in ARRAY or LIST mode.
//! Handles unphased (genotype_bytes) and phased (phased_pairs) output.
//! For unphased: writes int8 genotype values, -9 → NULL.
//! For phased: writes ARRAY(TINYINT, 2) pairs, -9 → NULL pair.
//! Does NOT handle COLUMNS mode (caller handles that per-column).
void FillGenotypeVector(
    Vector &vec,
    idx_t row_idx,
    GenotypeMode mode,
    uint32_t output_sample_ct,
    const int8_t *genotype_bytes,
    const int8_t *phased_pairs,
    bool include_phased);
```

- [ ] **Step 2: Implement in `src/plink_common.cpp`**

```cpp
void FillGenotypeVector(
    Vector &vec,
    idx_t row_idx,
    GenotypeMode mode,
    uint32_t output_sample_ct,
    const int8_t *genotype_bytes,
    const int8_t *phased_pairs,
    bool include_phased) {

    if (include_phased) {
        if (mode == GenotypeMode::ARRAY) {
            auto &pair_vec = ArrayVector::GetEntry(vec);
            auto &allele_vec = ArrayVector::GetEntry(pair_vec);
            auto *allele_data = FlatVector::GetData<int8_t>(allele_vec);
            auto &pair_validity = FlatVector::Validity(pair_vec);

            idx_t pair_base = row_idx * static_cast<idx_t>(output_sample_ct);
            for (idx_t s = 0; s < output_sample_ct; s++) {
                idx_t pair_idx = pair_base + s;
                idx_t allele_base = pair_idx * 2;
                int8_t a1 = phased_pairs[s * 2];
                int8_t a2 = phased_pairs[s * 2 + 1];
                if (a1 == -9) {
                    pair_validity.SetInvalid(pair_idx);
                    allele_data[allele_base] = 0;
                    allele_data[allele_base + 1] = 0;
                } else {
                    allele_data[allele_base] = a1;
                    allele_data[allele_base + 1] = a2;
                }
            }
        } else {
            // LIST mode
            auto list_offset = ListVector::GetListSize(vec);
            ListVector::Reserve(vec, list_offset + output_sample_ct);
            auto &pair_vec = ListVector::GetEntry(vec);
            auto &allele_vec = ArrayVector::GetEntry(pair_vec);
            auto *allele_data = FlatVector::GetData<int8_t>(allele_vec);
            auto &pair_validity = FlatVector::Validity(pair_vec);

            for (idx_t s = 0; s < output_sample_ct; s++) {
                idx_t pair_idx = list_offset + s;
                idx_t allele_base = pair_idx * 2;
                int8_t a1 = phased_pairs[s * 2];
                int8_t a2 = phased_pairs[s * 2 + 1];
                if (a1 == -9) {
                    pair_validity.SetInvalid(pair_idx);
                    allele_data[allele_base] = 0;
                    allele_data[allele_base + 1] = 0;
                } else {
                    allele_data[allele_base] = a1;
                    allele_data[allele_base + 1] = a2;
                }
            }

            auto *list_data = FlatVector::GetData<list_entry_t>(vec);
            list_data[row_idx].offset = list_offset;
            list_data[row_idx].length = output_sample_ct;
            ListVector::SetListSize(vec, list_offset + output_sample_ct);
        }
    } else {
        // Unphased
        if (mode == GenotypeMode::ARRAY) {
            auto array_size = static_cast<idx_t>(output_sample_ct);
            auto &child = ArrayVector::GetEntry(vec);
            auto *child_data = FlatVector::GetData<int8_t>(child);
            auto &child_validity = FlatVector::Validity(child);

            idx_t base = row_idx * array_size;
            for (idx_t s = 0; s < array_size; s++) {
                int8_t geno = genotype_bytes[s];
                if (geno == -9) {
                    child_validity.SetInvalid(base + s);
                    child_data[base + s] = 0;
                } else {
                    child_data[base + s] = geno;
                }
            }
        } else {
            // LIST mode
            auto list_offset = ListVector::GetListSize(vec);
            ListVector::Reserve(vec, list_offset + output_sample_ct);
            auto &child = ListVector::GetEntry(vec);
            auto *child_data = FlatVector::GetData<int8_t>(child);
            auto &child_validity = FlatVector::Validity(child);
            for (idx_t s = 0; s < output_sample_ct; s++) {
                int8_t geno = genotype_bytes[s];
                if (geno == -9) {
                    child_validity.SetInvalid(list_offset + s);
                    child_data[list_offset + s] = 0;
                } else {
                    child_data[list_offset + s] = geno;
                }
            }
            auto *list_data = FlatVector::GetData<list_entry_t>(vec);
            list_data[row_idx].offset = list_offset;
            list_data[row_idx].length = output_sample_ct;
            ListVector::SetListSize(vec, list_offset + output_sample_ct);
        }
    }
}
```

- [ ] **Step 3: Replace inline logic in `pfile_reader.cpp`**

In `pfile_reader.cpp`, the variant-orient scan function's `case PfileBindData::GENOTYPES_COL:` section (around lines 1909-2043) has the ARRAY/LIST output code. Replace the phased + unphased ARRAY/LIST blocks with calls to `FillGenotypeVector`.

Keep the following inline (they don't apply to vcf_reader):
- COLUMNS mode dispatch (in the `default:` case)
- STRUCT mode
- COUNTS/STATS mode (uses PgrGetCounts, not genotype bytes)
- Dosage ARRAY/LIST mode
- GenotypeRangeFilter per-element checks (v1 vcf_reader doesn't have this)

The replacement for the unphased+phased ARRAY/LIST block (~lines 1909-2043):

```cpp
if (bind_data.include_phased) {
    FillGenotypeVector(vec, rows_emitted, bind_data.genotype_mode,
                       output_sample_ct, nullptr, lstate.phased_pairs.data(), true);
} else {
    FillGenotypeVector(vec, rows_emitted, bind_data.genotype_mode,
                       output_sample_ct, lstate.genotype_bytes.data(), nullptr, false);
}
```

**Important:** The existing code has `genotype_filter` (GenotypeRangeFilter) checks interleaved with output. Since vcf_reader v1 doesn't need genotype_range, and removing it from pfile_reader would break existing behavior, we have two choices:
1. Keep the filter logic inline in pfile_reader and don't extract it (simpler, vcf_reader just doesn't call with filter)
2. Add filter parameter to FillGenotypeVector

**Decision:** Keep filter inline in pfile_reader for now. The shared function handles the common no-filter case. pfile_reader keeps its existing code for filtered paths. vcf_reader calls the shared function directly.

This means the refactor is more conservative: we extract the *no-filter* path into the shared function, and pfile_reader uses it only when `!genotype_filter.active`. When genotype_filter IS active, pfile_reader keeps its existing inline code.

- [ ] **Step 4: Build and run existing tests**

```bash
make -j4 && make test
```

All existing pfile tests must pass unchanged.

- [ ] **Step 5: Commit**

```bash
git add src/include/plink_common.hpp src/plink_common.cpp src/pfile_reader.cpp
git commit -m "refactor: extract shared FillGenotypeVector from pfile_reader"
```

---

### Task 4: Implement VCF Reader — Bind Phase

**Files:**
- Create: `src/include/vcf_reader.hpp`
- Create: `src/vcf_reader.cpp` (partial — bind function only)

- [ ] **Step 1: Create `src/include/vcf_reader.hpp`**

```cpp
#pragma once

#include "duckdb.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

void RegisterPlinkVcfReader(ExtensionLoader &loader);

} // namespace duckdb
```

- [ ] **Step 2: Create `src/vcf_reader.cpp` with bind data structures and bind function**

The bind function:
1. Opens the VCF file via DuckDB FileSystem
2. Reads header lines (## metadata) — scans for FORMAT=GT/GQ/DP definitions
3. Parses #CHROM header line → extracts sample names
4. Resolves genotype output mode (auto/array/list/columns)
5. Builds return schema: chrom(VARCHAR), pos(INTEGER), id(VARCHAR), ref(VARCHAR), alt(VARCHAR), genotypes(type depends on mode)

Key data structures:

```cpp
struct VcfBindData : public TableFunctionData {
    string file_path;
    vector<string> sample_names;
    uint32_t sample_ct;

    // Output configuration
    GenotypeMode genotype_mode;
    bool include_phased;

    // Region filter (parsed from 'region' parameter)
    string filter_chrom;     // empty = no filter
    int32_t filter_start;    // 1-based, inclusive
    int32_t filter_end;      // 1-based, inclusive

    // Quality filter config
    int32_t min_gq;          // -1 = disabled
    int32_t min_dp;          // -1 = disabled
    int32_t max_dp;          // -1 = disabled

    // Half-call mode
    VcfHalfCall halfcall_mode;

    // Column layout (for columns mode)
    vector<string> genotype_column_names;
    idx_t columns_mode_first_geno_col;
    idx_t columns_mode_geno_col_count;

    // Column IDs for projection pushdown
    static constexpr idx_t CHROM_COL = 0;
    static constexpr idx_t POS_COL = 1;
    static constexpr idx_t ID_COL = 2;
    static constexpr idx_t REF_COL = 3;
    static constexpr idx_t ALT_COL = 4;
    static constexpr idx_t GENOTYPES_COL = 5;
};
```

The bind function parses the VCF header to extract sample names and validates the file format. It does NOT read any data lines.

- [ ] **Step 3: Build to verify compilation**

```bash
make -j4
```

- [ ] **Step 4: Commit**

```bash
git add src/vcf_reader.cpp src/include/vcf_reader.hpp
git commit -m "feat: add VCF reader bind phase (header parsing, schema setup)"
```

---

### Task 5: Implement VCF Reader — Init and Scan

**Files:**
- Modify: `src/vcf_reader.cpp` (add init + scan functions)

- [ ] **Step 1: Add global/local state structs and init functions**

```cpp
struct VcfGlobalState : public GlobalTableFunctionState {
    mutex lock;
    unique_ptr<FileHandle> file_handle;
    bool done;
    uint64_t multiallelic_skipped;
};

struct VcfLocalState : public LocalTableFunctionState {
    // Genotype parsing buffers
    AlignedBuffer genovec_buf;
    AlignedBuffer phasepresent_buf;
    AlignedBuffer phaseinfo_buf;
    vector<int8_t> genotype_bytes;
    vector<int8_t> phased_pairs;

    // Parse context
    VcfParseContext parse_ctx;

    // Line buffer
    string line;
};
```

InitGlobal opens the file, seeks past the header (re-reads header lines to position at first data line). InitLocal allocates the genovec buffer using `NypCtToAlignedWordCt(sample_ct)`.

- [ ] **Step 2: Implement the scan function**

The scan function reads lines from the file (serialized via mutex on the global state's file handle), parses each line:

1. Split on tabs to find the 9 fixed columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
2. Check region filter (compare chrom + pos)
3. Check if multiallelic (ALT contains ',') → skip
4. Parse FORMAT field to locate GT position and quality field positions
5. Position pointer at first sample data
6. Call `ParseUnphasedBiallelicGT` or `ParsePhasedBiallelicGT`
7. Convert genovec → genotype_bytes via `GenoarrToBytesMinus9`
8. If phased: call `UnpackPhasedGenotypes`
9. Fill output vectors:
   - Fixed columns: write CHROM/POS/ID/REF/ALT directly
   - Genotypes: call `FillGenotypeVector` (or fill columns inline for COLUMNS mode)

Projection pushdown: check `column_ids` — if GENOTYPES_COL is not referenced, skip steps 4-8 entirely.

- [ ] **Step 3: Register the function**

In the same file, add the registration function:

```cpp
void RegisterPlinkVcfReader(ExtensionLoader &loader) {
    TableFunction func("read_plink_vcf", {LogicalType::VARCHAR}, VcfScanFunction,
                       VcfBind, VcfInitGlobal, VcfInitLocal);
    func.named_parameters["genotypes"] = LogicalType::VARCHAR;
    func.named_parameters["phased"] = LogicalType::BOOLEAN;
    func.named_parameters["region"] = LogicalType::VARCHAR;
    func.named_parameters["min_gq"] = LogicalType::INTEGER;
    func.named_parameters["min_dp"] = LogicalType::INTEGER;
    func.named_parameters["max_dp"] = LogicalType::INTEGER;
    func.named_parameters["halfcall"] = LogicalType::VARCHAR;
    func.projection_pushdown = true;

    ExtensionUtil::RegisterFunction(loader, func);
}
```

- [ ] **Step 4: Wire into extension entry point**

In `src/plinking_duck_extension.cpp`, add:
```cpp
#include "vcf_reader.hpp"
```
And in `PlinkingDuckExtension::Load()`:
```cpp
RegisterPlinkVcfReader(loader);
```

- [ ] **Step 5: Update CMakeLists.txt**

Add `src/vcf_reader.cpp` and `src/vcf_genotype_parse.cpp` to `EXTENSION_SOURCES`:

```cmake
set(EXTENSION_SOURCES
    src/plinking_duck_extension.cpp
    src/pvar_reader.cpp
    src/psam_reader.cpp
    src/pgen_reader.cpp
    src/pfile_reader.cpp
    src/plink_common.cpp
    src/plink_freq.cpp
    src/plink_hardy.cpp
    src/plink_missing.cpp
    src/plink_ld.cpp
    src/plink_score.cpp
    src/plink_glm.cpp
    src/plink2_glm_logistic_math.cpp
    src/vcf_genotype_parse.cpp
    src/vcf_reader.cpp
)
```

- [ ] **Step 6: Build**

```bash
make -j4
```

- [ ] **Step 7: Commit**

```bash
git add src/vcf_reader.cpp src/plinking_duck_extension.cpp CMakeLists.txt
git commit -m "feat: implement VCF reader init/scan with plink-ng GT parsing"
```

---

### Task 6: Write Positive Tests

**Files:**
- Create: `test/sql/read_plink_vcf.test`

- [ ] **Step 1: Write basic positive tests**

```sql
# name: test/sql/read_plink_vcf.test
# description: Positive tests for read_plink_vcf table function
# group: [sql]

require plinking_duck

# --- Basic read (uses existing example.vcf which has same genotypes as pfile_example) ---

# Row count
query I
SELECT COUNT(*) FROM read_plink_vcf('test/data/example.vcf');
----
4

# Variant metadata columns
query TITTT
SELECT CHROM, POS, ID, REF, ALT FROM read_plink_vcf('test/data/example.vcf') ORDER BY POS;
----
1	10000	rs1	A	G
1	20000	rs2	C	T
1	30000	rs3	G	A
2	15000	rs4	T	C

# Genotype array length matches sample count
query I
SELECT len(genotypes) FROM read_plink_vcf('test/data/example.vcf') LIMIT 1;
----
4

# Known genotype values match pfile_example exactly
query TT
SELECT ID, genotypes FROM read_plink_vcf('test/data/example.vcf') WHERE ID = 'rs1';
----
rs1	[0, 1, 2, NULL]

query TT
SELECT ID, genotypes FROM read_plink_vcf('test/data/example.vcf') WHERE ID = 'rs2';
----
rs2	[1, 1, 0, 2]

query TT
SELECT ID, genotypes FROM read_plink_vcf('test/data/example.vcf') WHERE ID = 'rs3';
----
rs3	[2, NULL, 1, 0]

query TT
SELECT ID, genotypes FROM read_plink_vcf('test/data/example.vcf') WHERE ID = 'rs4';
----
rs4	[0, 0, 1, 2]

# --- LIST mode ---

query TT
SELECT ID, genotypes FROM read_plink_vcf('test/data/example.vcf', genotypes := 'list') WHERE ID = 'rs1';
----
rs1	[0, 1, 2, NULL]

# --- Projection pushdown (no genotypes column) ---

query TI
SELECT CHROM, POS FROM read_plink_vcf('test/data/example.vcf') ORDER BY POS;
----
1	10000
1	20000
1	30000
2	15000

# --- Columns mode ---

query IIII
SELECT SAMPLE1, SAMPLE2, SAMPLE3, SAMPLE4 FROM read_plink_vcf('test/data/example.vcf', genotypes := 'columns') WHERE ID = 'rs1';
----
0	1	2	NULL

# --- Region filter ---

query I
SELECT COUNT(*) FROM read_plink_vcf('test/data/example.vcf', region := '1');
----
3

query I
SELECT COUNT(*) FROM read_plink_vcf('test/data/example.vcf', region := '1:15000-25000');
----
1

# --- Phased mode ---

query TT
SELECT ID, genotypes FROM read_plink_vcf('test/data/vcf_phased.vcf', phased := true) WHERE ID = 'rs1';
----
rs1	[[0, 0], [0, 1], [1, 0], [1, 1]]

query TT
SELECT ID, genotypes FROM read_plink_vcf('test/data/vcf_phased.vcf', phased := true) WHERE ID = 'rs2';
----
rs2	[[0, 1], [1, 0], [0, 0], NULL]

# --- Quality filtering ---

# min_gq=20 should set rs2/SAMPLE2 (GQ=5) and rs4/SAMPLE3 (GQ=15) to missing
query TT
SELECT ID, genotypes FROM read_plink_vcf('test/data/vcf_example.vcf', min_gq := 20) WHERE ID = 'rs2';
----
rs2	[1, NULL, 0, 2]

query TT
SELECT ID, genotypes FROM read_plink_vcf('test/data/vcf_example.vcf', min_gq := 20) WHERE ID = 'rs4';
----
rs4	[0, 0, NULL, 2]
```

- [ ] **Step 2: Run tests**

```bash
make test
```

All tests should pass. Debug any failures.

- [ ] **Step 3: Commit**

```bash
git add test/sql/read_plink_vcf.test
git commit -m "test: add positive tests for read_plink_vcf"
```

---

### Task 7: Write Negative Tests

**Files:**
- Create: `test/sql/read_plink_vcf_negative.test`

- [ ] **Step 1: Write error case tests**

```sql
# name: test/sql/read_plink_vcf_negative.test
# description: Negative tests for read_plink_vcf table function
# group: [sql]

require plinking_duck

# --- File not found ---
statement error
SELECT * FROM read_plink_vcf('test/data/nonexistent.vcf');
----
IO Error

# --- Invalid genotypes parameter ---
statement error
SELECT * FROM read_plink_vcf('test/data/example.vcf', genotypes := 'invalid');
----
Invalid Input Error

# --- Invalid halfcall parameter ---
statement error
SELECT * FROM read_plink_vcf('test/data/example.vcf', halfcall := 'invalid');
----
Invalid Input Error

# --- Invalid region format ---
statement error
SELECT * FROM read_plink_vcf('test/data/example.vcf', region := 'chr1:abc-def');
----
Invalid Input Error

# --- Phased + columns is valid (unlike phased + dosages) ---
query I
SELECT COUNT(*) FROM read_plink_vcf('test/data/vcf_phased.vcf', phased := true, genotypes := 'columns');
----
3
```

- [ ] **Step 2: Run tests**

```bash
make test
```

- [ ] **Step 3: Commit**

```bash
git add test/sql/read_plink_vcf_negative.test
git commit -m "test: add negative tests for read_plink_vcf"
```

---

### Task 8: Multiallelic Skip + Cross-Validation

**Files:**
- Create: `test/data/vcf_multiallelic.vcf`
- Modify: `test/sql/read_plink_vcf.test` (add multiallelic + cross-validation tests)

- [ ] **Step 1: Create test VCF with mixed biallelic/multiallelic variants**

```
##fileformat=VCFv4.3
##contig=<ID=1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3	SAMPLE4
1	10000	rs1	A	G	.	.	.	GT	0/0	0/1	1/1	./.
1	15000	rs_multi	A	G,T	.	.	.	GT	0/1	0/2	1/2	0/0
1	20000	rs2	C	T	.	.	.	GT	0/1	0/1	0/0	1/1
```

- [ ] **Step 2: Add tests verifying multiallelic rows are skipped**

```sql
# Multiallelic variants are skipped (only biallelic rows emitted)
query I
SELECT COUNT(*) FROM read_plink_vcf('test/data/vcf_multiallelic.vcf');
----
2

query T
SELECT ID FROM read_plink_vcf('test/data/vcf_multiallelic.vcf') ORDER BY POS;
----
rs1
rs2
```

- [ ] **Step 3: Add cross-validation test (VCF genotypes = pfile genotypes)**

```sql
# Cross-validate: read_plink_vcf on example.vcf should produce identical
# genotypes to read_pfile on pfile_example
query I
SELECT COUNT(*) FROM (
    SELECT v.ID, v.genotypes as vcf_geno, p.genotypes as pfile_geno
    FROM read_plink_vcf('test/data/example.vcf') v
    JOIN read_pfile('test/data/pfile_example') p ON v.ID = p.ID
    WHERE v.genotypes != p.genotypes OR (v.genotypes IS NULL) != (p.genotypes IS NULL)
);
----
0
```

- [ ] **Step 4: Run tests**

```bash
make test
```

- [ ] **Step 5: Commit**

```bash
git add test/data/vcf_multiallelic.vcf test/sql/read_plink_vcf.test
git commit -m "test: add multiallelic skip and cross-validation tests"
```

---

### Task 9: Final Integration + Cleanup

**Files:**
- Verify all tests pass
- Clean up any TODO comments or debug output

- [ ] **Step 1: Full build + test**

```bash
make clean && make -j4 && make test
```

- [ ] **Step 2: Verify no warnings**

Check build output for any compiler warnings in the new files. Fix any warnings.

- [ ] **Step 3: Test with compressed VCF**

Create a gzipped version of the test VCF and verify it works:

```bash
gzip -k test/data/example.vcf  # creates example.vcf.gz
```

Then in duckdb:
```sql
SELECT COUNT(*) FROM read_plink_vcf('test/data/example.vcf.gz');
-- Should return 4
```

If DuckDB's FileSystem doesn't auto-detect .vcf.gz, we may need to use `FileOpenFlags::FILE_FLAGS_READ | FileOpenFlags::FILE_FLAGS_ALLOW_COMPRESSION` or similar. Check how pvar_reader opens files.

- [ ] **Step 4: Add compressed VCF test to test file**

```sql
# Compressed VCF (.vcf.gz)
query I
SELECT COUNT(*) FROM read_plink_vcf('test/data/example.vcf.gz');
----
4
```

- [ ] **Step 5: Final commit**

```bash
git add -A
git commit -m "feat: complete read_plink_vcf implementation (closes #11)"
```
