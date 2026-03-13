# P4-003 Phase 2: genotype_range Filter Pushdown — Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `genotype_range` parameter that filters individual genotype values with a pre-decompression optimization using PgrGetCounts to skip entire variants.

**Architecture:** Reuse the `RangeFilter`/`ParseRangeFilter` infrastructure from Phase 1. Add `CheckGenotypeRange` that maps genotype values {0,1,2} to genocounts to determine `any_pass`/`all_pass` without decompression. Integrate the filter into all scan paths: variant orient (parallel), genotype orient (single-threaded), and sample orient (bind-time pre-pass).

**Tech Stack:** C++, DuckDB extension API, pgenlib (PgrGetCounts, PgrGet, PgrGetP)

**Spec:** `docs/superpowers/specs/2026-03-12-genotype-range-filter-design.md`

---

## File Structure

| File | Responsibility |
|------|---------------|
| `src/include/plink_common.hpp` | Add `GenotypeRangeResult` struct, `CheckGenotypeRange` declaration, `GenotypeRangeFilter` struct |
| `src/plink_common.cpp` | Implement `CheckGenotypeRange` |
| `src/pfile_reader.cpp` | Add `genotype_filter` to bind data, parse `genotype_range`, extend `need_pgen_reader`, integrate into all 3 scan paths + sample orient pre-pass |
| `src/pgen_reader.cpp` | Mirror changes for variant orient only |
| `test/sql/read_pfile_filter.test` | Add genotype_range test cases |
| `test/sql/read_pgen_filter.test` | Add genotype_range test cases |

---

## Chunk 1: Core Infrastructure (plink_common)

### Task 1: GenotypeRangeFilter struct and CheckGenotypeRange

**Files:**
- Modify: `src/include/plink_common.hpp:222-249` (after CountFilter section)
- Modify: `src/plink_common.cpp:615` (after VariantPassesCountFilter)
- Test: `test/sql/read_pfile_filter.test` (add genotype_range tests after existing)
- Test: `test/sql/read_pgen_filter.test` (add genotype_range tests after existing)

- [ ] **Step 1: Add GenotypeRangeFilter struct and CheckGenotypeRange declaration to header**

In `src/include/plink_common.hpp`, after the CountFilter section (after line 249), add:

```cpp
// ---------------------------------------------------------------------------
// Genotype range filtering (genotype_range)
// ---------------------------------------------------------------------------

struct GenotypeRangeResult {
    bool any_pass;  // at least one non-missing sample has value in [min, max]
    bool all_pass;  // every non-missing sample has value in [min, max]
};

struct GenotypeRangeFilter {
    RangeFilter range;  // reuses Phase 1 RangeFilter
    bool active = false;
};

GenotypeRangeResult CheckGenotypeRange(const RangeFilter &filter,
                                        const STD_ARRAY_REF(uint32_t, 4) genocounts);
```

- [ ] **Step 2: Implement CheckGenotypeRange in plink_common.cpp**

In `src/plink_common.cpp`, after `VariantPassesCountFilter` (after line 615), add:

```cpp
// ---------------------------------------------------------------------------
// Genotype range filtering (genotype_range)
// ---------------------------------------------------------------------------

GenotypeRangeResult CheckGenotypeRange(const RangeFilter &filter,
                                        const STD_ARRAY_REF(uint32_t, 4) genocounts) {
    GenotypeRangeResult result;
    result.any_pass = false;
    result.all_pass = true;

    // genocounts: [hom_ref(0), het(1), hom_alt(2), missing(3)]
    // Genotype values map directly to array indices 0, 1, 2
    for (int g = 0; g <= 2; g++) {
        bool in_range = filter.Passes(static_cast<double>(g));
        if (in_range && genocounts[g] > 0) {
            result.any_pass = true;
        }
        if (!in_range && genocounts[g] > 0) {
            result.all_pass = false;
        }
    }

    return result;
}
```

- [ ] **Step 3: Build to verify compilation**

Run: `blq run build`
Expected: Clean compilation

- [ ] **Step 4: Commit**

```bash
git add src/include/plink_common.hpp src/plink_common.cpp
git commit -m "feat: add GenotypeRangeFilter struct and CheckGenotypeRange function

Core infrastructure for genotype_range filter pushdown (P4-003 Phase 2).
CheckGenotypeRange maps genotype values {0,1,2} to genocounts to
determine any_pass/all_pass without decompression."
```

---

## Chunk 2: Parameter Parsing and Bind Integration

### Task 2: Add genotype_filter to PfileBindData and parse genotype_range

**Files:**
- Modify: `src/pfile_reader.cpp:229-231` (bind data, after count_filter)
- Modify: `src/pfile_reader.cpp:684-703` (parse section, after count filter parsing)
- Modify: `src/pfile_reader.cpp:397-398` (parameter loop comment)
- Modify: `src/pfile_reader.cpp:1190` (need_pgen_reader)
- Modify: `src/pfile_reader.cpp:2175-2176` (registration)

- [ ] **Step 1: Add genotype_filter field to PfileBindData**

In `src/pfile_reader.cpp`, after the `count_filter_subset` field (line 231), add:

```cpp
    // Genotype range filtering (genotype_range)
    GenotypeRangeFilter genotype_filter;
```

- [ ] **Step 2: Parse genotype_range parameter in PfileBind**

In `src/pfile_reader.cpp`, after the count filter parsing block (after line 703, before `// --- Build output schema ---`), add:

```cpp
    // --- Parse genotype_range filter ---
    {
        auto gr_it = input.named_parameters.find("genotype_range");
        if (gr_it != input.named_parameters.end()) {
            if (bind_data->include_dosages) {
                throw InvalidInputException("read_pfile: genotype_range is incompatible with dosages := true");
            }
            bind_data->genotype_filter.range = ParseRangeFilter(
                gr_it->second, "genotype_range", 0.0, 2.0, "read_pfile");
            bind_data->genotype_filter.active = bind_data->genotype_filter.range.active;
        }
    }
```

- [ ] **Step 3: Extend need_pgen_reader in PfileInitGlobal**

In `src/pfile_reader.cpp`, change line 1190 from:
```cpp
    state->need_pgen_reader = state->need_genotypes || bind_data.count_filter.HasFilter();
```
to:
```cpp
    state->need_pgen_reader = state->need_genotypes
                           || bind_data.count_filter.HasFilter()
                           || bind_data.genotype_filter.active;
```

- [ ] **Step 4: Register genotype_range named parameter**

In `src/pfile_reader.cpp`, after the `ac_range` registration (line 2176), add:

```cpp
    read_pfile.named_parameters["genotype_range"] = LogicalType::ANY;
```

- [ ] **Step 5: Build to verify compilation**

Run: `blq run build`
Expected: Clean compilation

- [ ] **Step 6: Commit**

```bash
git add src/pfile_reader.cpp
git commit -m "feat: parse genotype_range parameter in read_pfile bind

Adds GenotypeRangeFilter to PfileBindData, parses the parameter with
validation (domain [0,2], incompatible with dosages), extends
need_pgen_reader condition."
```

### Task 3: Add genotype_filter to PgenBindData and parse genotype_range

**Files:**
- Modify: `src/pgen_reader.cpp:46-48` (bind data)
- Modify: `src/pgen_reader.cpp:254-273` (parse section)
- Modify: `src/pgen_reader.cpp:357` (need_pgen_reader)
- Modify: `src/pgen_reader.cpp:841-842` (registration)

- [ ] **Step 1: Add genotype_filter field to PgenBindData**

In `src/pgen_reader.cpp`, after `count_filter_subset` (line 48), add:

```cpp
    // Genotype range filtering (genotype_range)
    GenotypeRangeFilter genotype_filter;
```

- [ ] **Step 2: Parse genotype_range parameter in PgenBind**

In `src/pgen_reader.cpp`, after the count filter parsing block (after line 273, before `// --- Resolve genotype output mode ---`), add:

```cpp
    // --- Parse genotype_range filter ---
    {
        auto gr_it = input.named_parameters.find("genotype_range");
        if (gr_it != input.named_parameters.end()) {
            if (bind_data->include_dosages) {
                throw InvalidInputException("read_pgen: genotype_range is incompatible with dosages := true");
            }
            bind_data->genotype_filter.range = ParseRangeFilter(
                gr_it->second, "genotype_range", 0.0, 2.0, "read_pgen");
            bind_data->genotype_filter.active = bind_data->genotype_filter.range.active;
        }
    }
```

- [ ] **Step 3: Extend need_pgen_reader in PgenInitGlobal**

In `src/pgen_reader.cpp`, change line 357 from:
```cpp
    state->need_pgen_reader = state->need_genotypes || bind_data.count_filter.HasFilter();
```
to:
```cpp
    state->need_pgen_reader = state->need_genotypes
                           || bind_data.count_filter.HasFilter()
                           || bind_data.genotype_filter.active;
```

- [ ] **Step 4: Register genotype_range named parameter**

In `src/pgen_reader.cpp`, after the `ac_range` registration (line 842), add:

```cpp
    read_pgen.named_parameters["genotype_range"] = LogicalType::ANY;
```

- [ ] **Step 5: Build to verify compilation**

Run: `blq run build`
Expected: Clean compilation

- [ ] **Step 6: Commit**

```bash
git add src/pgen_reader.cpp
git commit -m "feat: parse genotype_range parameter in read_pgen bind"
```

---

## Chunk 3: Variant Orient Scan Integration

### Task 4: Add genotype_range filtering to PfileDefaultScan (variant orient)

**Files:**
- Modify: `src/pfile_reader.cpp:1355-1371` (PfileDefaultScan count filter block)
- Modify: `src/pfile_reader.cpp:1581-1618` (unphased genotype emit loop)
- Modify: `src/pfile_reader.cpp:1494-1541` (phased genotype emit loop)
- Test: `test/sql/read_pfile_filter.test`

The integration has two parts:
1. **Pre-decompression**: After PgrGetCounts (reuse the existing count filter call), check `any_pass`/`all_pass`
2. **Per-element filtering**: In the genotype emit loops, NULL out values outside the range when `!all_pass`

- [ ] **Step 1: Extend the count filter block to include genotype_range pre-decompression check**

In `src/pfile_reader.cpp`, replace the count filter block in PfileDefaultScan (lines 1355-1371) with a combined block that handles both count filter and genotype range:

```cpp
            // Count filter + genotype range pre-decompression check
            bool geno_range_all_pass = true;
            if ((bind_data.count_filter.HasFilter() || bind_data.genotype_filter.active) && lstate.initialized) {
                STD_ARRAY_DECL(uint32_t, 4, genocounts);
                const uintptr_t *cf_si = (bind_data.has_sample_subset && bind_data.count_filter_subset)
                                             ? bind_data.count_filter_subset->SampleInclude() : nullptr;
                const uintptr_t *cf_iv = (bind_data.has_sample_subset && bind_data.count_filter_subset)
                                             ? bind_data.count_filter_subset->InterleavedVec() : nullptr;
                uint32_t cf_sc = bind_data.has_sample_subset ? bind_data.subset_sample_ct : bind_data.raw_sample_ct;
                plink2::PglErr cf_err = plink2::PgrGetCounts(cf_si, cf_iv, lstate.pssi, cf_sc, vidx,
                                                              &lstate.pgr, genocounts);
                if (cf_err != plink2::kPglRetSuccess) {
                    throw IOException("read_pfile: PgrGetCounts failed for variant %u", vidx);
                }
                if (bind_data.count_filter.HasFilter() &&
                    !VariantPassesCountFilter(bind_data.count_filter, genocounts, cf_sc)) {
                    continue;
                }
                if (bind_data.genotype_filter.active) {
                    auto gr = CheckGenotypeRange(bind_data.genotype_filter.range, genocounts);
                    if (!gr.any_pass) {
                        continue;
                    }
                    geno_range_all_pass = gr.all_pass;
                }
            }
```

- [ ] **Step 2: Add per-element genotype_range filtering in unphased ARRAY emit**

In the unphased ARRAY emit loop (currently starting around line 1582), after:
```cpp
                        int8_t geno = lstate.genotype_bytes[s];
                        if (geno == -9) {
```
Add genotype_range check. The modified loop becomes:
```cpp
                            idx_t base = rows_emitted * array_size;
                            for (idx_t s = 0; s < array_size; s++) {
                                int8_t geno = lstate.genotype_bytes[s];
                                if (geno == -9 ||
                                    (bind_data.genotype_filter.active && !geno_range_all_pass &&
                                     !bind_data.genotype_filter.range.Passes(static_cast<double>(geno)))) {
                                    child_validity.SetInvalid(base + s);
                                    child_data[base + s] = 0;
                                } else {
                                    child_data[base + s] = geno;
                                }
                            }
```

Apply the same pattern to:
- Unphased LIST emit loop (around line 1600)
- Phased ARRAY emit loop (around line 1495) — check `allele1 + allele2` against range, NULL the pair
- Phased LIST emit loop (around line 1516) — same as phased ARRAY
- COLUMNS mode emit (around line 1621) — check scalar geno against range
- Dosage paths are excluded (incompatible with genotype_range, checked at bind)

For phased mode, the range check is on the diploid sum:
```cpp
                                int8_t a1 = lstate.phased_pairs[s * 2];
                                int8_t a2 = lstate.phased_pairs[s * 2 + 1];
                                if (a1 == -9 ||
                                    (bind_data.genotype_filter.active && !geno_range_all_pass &&
                                     !bind_data.genotype_filter.range.Passes(static_cast<double>(a1 + a2)))) {
                                    pair_validity.SetInvalid(pair_idx);
                                    allele_data[allele_base] = 0;
                                    allele_data[allele_base + 1] = 0;
                                } else {
                                    allele_data[allele_base] = a1;
                                    allele_data[allele_base + 1] = a2;
                                }
```

- [ ] **Step 3: Write test: carriers only variant orient**

Append to `test/sql/read_pfile_filter.test`:

```sql
# ========== genotype_range tests ==========

# --- Test data reference (pgen_example, 4 variants x 4 samples) ---
# rs1: genotypes [0, 1, 2, NULL]
# rs2: genotypes [1, 1, 0, 2]
# rs3: genotypes [2, NULL, 1, 0]
# rs4: genotypes [0, 0, 1, 2]

# --- Carriers only (min: 1): all 4 variants survive, hom-ref and missing NULLed ---

query TT
SELECT ID, genotypes FROM read_pfile('test/data/pfile_example', genotype_range := {min: 1}) ORDER BY ID;
----
rs1	[NULL, 1, 2, NULL]
rs2	[1, 1, NULL, 2]
rs3	[2, NULL, 1, NULL]
rs4	[NULL, NULL, 1, 2]

# --- Hets only (min: 1, max: 1) ---

query TT
SELECT ID, genotypes FROM read_pfile('test/data/pfile_example', genotype_range := {min: 1, max: 1}) ORDER BY ID;
----
rs1	[NULL, 1, NULL, NULL]
rs2	[1, 1, NULL, NULL]
rs3	[NULL, NULL, 1, NULL]
rs4	[NULL, NULL, 1, NULL]

# --- Hom-alt only (min: 2, max: 2) ---

query TT
SELECT ID, genotypes FROM read_pfile('test/data/pfile_example', genotype_range := {min: 2, max: 2}) ORDER BY ID;
----
rs1	[NULL, NULL, 2, NULL]
rs2	[NULL, NULL, NULL, 2]
rs3	[2, NULL, NULL, NULL]
rs4	[NULL, NULL, NULL, 2]

# --- Hom-ref only (min: 0, max: 0): missing excluded ---

query TT
SELECT ID, genotypes FROM read_pfile('test/data/pfile_example', genotype_range := {min: 0, max: 0}) ORDER BY ID;
----
rs1	[0, NULL, NULL, NULL]
rs2	[NULL, NULL, 0, NULL]
rs3	[NULL, NULL, NULL, 0]
rs4	[0, 0, NULL, NULL]
```

- [ ] **Step 4: Build and run tests**

Run: `blq run build && blq run test`
Expected: All tests pass including new genotype_range tests

- [ ] **Step 5: Commit**

```bash
git add src/pfile_reader.cpp test/sql/read_pfile_filter.test
git commit -m "feat: integrate genotype_range into PfileDefaultScan (variant orient)

Pre-decompression: PgrGetCounts → CheckGenotypeRange → skip if !any_pass.
Per-element: NULL out-of-range values when !all_pass. Works with
unphased, phased, ARRAY, LIST, and COLUMNS modes."
```

### Task 5: Add genotype_range filtering to PgenScan

**Files:**
- Modify: `src/pgen_reader.cpp:514-530` (PgenScan count filter block)
- Modify: `src/pgen_reader.cpp:743-781` (unphased ARRAY/LIST genotype emit loops)
- Modify: `src/pgen_reader.cpp:654-703` (phased ARRAY/LIST genotype emit loops)
- Modify: `src/pgen_reader.cpp:787-811` (COLUMNS mode emit in default case)
- Test: `test/sql/read_pgen_filter.test`

- [ ] **Step 1: Extend PgenScan count filter block with genotype_range check**

Same pattern as Task 4 Step 1. Replace lines 514-530 with combined block including `geno_range_all_pass` flag.

- [ ] **Step 2: Add per-element filtering in genotype emit loops**

Same pattern as Task 4 Step 2. Apply to all emit paths in PgenScan (ARRAY, LIST, COLUMNS, phased variants).

- [ ] **Step 3: Write tests**

Append to `test/sql/read_pgen_filter.test`:

```sql
# ========== genotype_range tests ==========

# --- Carriers only (min: 1) ---

query TT
SELECT ID, genotypes FROM read_pgen('test/data/pgen_example.pgen', genotype_range := {min: 1}) ORDER BY ID;
----
rs1	[NULL, 1, 2, NULL]
rs2	[1, 1, NULL, 2]
rs3	[2, NULL, 1, NULL]
rs4	[NULL, NULL, 1, 2]

# --- Hets only ---

query TT
SELECT ID, genotypes FROM read_pgen('test/data/pgen_example.pgen', genotype_range := {min: 1, max: 1}) ORDER BY ID;
----
rs1	[NULL, 1, NULL, NULL]
rs2	[1, 1, NULL, NULL]
rs3	[NULL, NULL, 1, NULL]
rs4	[NULL, NULL, 1, NULL]

# --- Hom-alt only ---

query TT
SELECT ID, genotypes FROM read_pgen('test/data/pgen_example.pgen', genotype_range := {min: 2, max: 2}) ORDER BY ID;
----
rs1	[NULL, NULL, 2, NULL]
rs2	[NULL, NULL, NULL, 2]
rs3	[2, NULL, NULL, NULL]
rs4	[NULL, NULL, NULL, 2]
```

- [ ] **Step 4: Build and run tests**

Run: `blq run build && blq run test`
Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add src/pgen_reader.cpp test/sql/read_pgen_filter.test
git commit -m "feat: integrate genotype_range into PgenScan (variant orient)"
```

---

## Chunk 4: Genotype Orient and Sample Orient

### Task 6: Add genotype_range filtering to PfileTidyScan (genotype orient)

**Files:**
- Modify: `src/pfile_reader.cpp:1730-1748` (PfileTidyScan count filter block)
- Modify: `src/pfile_reader.cpp:1793-1906` (sample emit loop)
- Test: `test/sql/read_pfile_filter.test`

In genotype orient, filtering means **skipping rows** (not NULLing elements) when genotype is out of range or missing.

- [ ] **Step 1: Extend PfileTidyScan count filter block with genotype_range check**

Replace lines 1730-1748 with combined block. Store `geno_range_all_pass` as a local variable accessible during sample iteration:

```cpp
        // Count filter + genotype range pre-decompression check
        bool geno_range_all_pass = true;
        if ((bind_data.count_filter.HasFilter() || bind_data.genotype_filter.active)
            && lstate.initialized && !lstate.geno_variant_loaded) {
            STD_ARRAY_DECL(uint32_t, 4, genocounts);
            const uintptr_t *cf_si = (bind_data.has_sample_subset && bind_data.count_filter_subset)
                                         ? bind_data.count_filter_subset->SampleInclude() : nullptr;
            const uintptr_t *cf_iv = (bind_data.has_sample_subset && bind_data.count_filter_subset)
                                         ? bind_data.count_filter_subset->InterleavedVec() : nullptr;
            uint32_t cf_sc = bind_data.has_sample_subset ? bind_data.subset_sample_ct : bind_data.raw_sample_ct;
            plink2::PglErr cf_err = plink2::PgrGetCounts(cf_si, cf_iv, lstate.pssi, cf_sc, vidx,
                                                          &lstate.pgr, genocounts);
            if (cf_err != plink2::kPglRetSuccess) {
                throw IOException("read_pfile: PgrGetCounts failed for variant %u", vidx);
            }
            if (bind_data.count_filter.HasFilter() &&
                !VariantPassesCountFilter(bind_data.count_filter, genocounts, cf_sc)) {
                lstate.geno_current_variant_pos++;
                lstate.geno_current_sample = 0;
                continue;
            }
            if (bind_data.genotype_filter.active) {
                auto gr = CheckGenotypeRange(bind_data.genotype_filter.range, genocounts);
                if (!gr.any_pass) {
                    lstate.geno_current_variant_pos++;
                    lstate.geno_current_sample = 0;
                    continue;
                }
                geno_range_all_pass = gr.all_pass;
            }
        }
```

- [ ] **Step 2: Add per-sample genotype_range skip in the sample emit loop**

In the sample emit loop (line 1793), insert at the top of the `while (lstate.geno_current_sample < ...)` loop, before line 1799. This skips rows for missing genotypes (always when genotype_range active) and out-of-range values (when `!all_pass`):

Note: In genotype orient mode, genotypes are always projected (`gstate.need_genotypes` is always true), so the `gstate.need_genotypes && lstate.geno_variant_loaded` guard is always satisfied after loading.

```cpp
            // When genotype_range is active, skip rows for missing/out-of-range genotypes
            if (bind_data.genotype_filter.active && gstate.need_genotypes && lstate.geno_variant_loaded) {
                if (bind_data.include_phased) {
                    int8_t a1 = lstate.phased_pairs[lstate.geno_current_sample * 2];
                    if (a1 == -9) {
                        lstate.geno_current_sample++;
                        continue;
                    }
                    if (!geno_range_all_pass) {
                        int8_t a2 = lstate.phased_pairs[lstate.geno_current_sample * 2 + 1];
                        if (!bind_data.genotype_filter.range.Passes(static_cast<double>(a1 + a2))) {
                            lstate.geno_current_sample++;
                            continue;
                        }
                    }
                } else if (!bind_data.include_dosages) {
                    int8_t geno = lstate.genotype_bytes[lstate.geno_current_sample];
                    if (geno == -9) {
                        lstate.geno_current_sample++;
                        continue;
                    }
                    if (!geno_range_all_pass &&
                        !bind_data.genotype_filter.range.Passes(static_cast<double>(geno))) {
                        lstate.geno_current_sample++;
                        continue;
                    }
                }
            }
```

- [ ] **Step 3: Write test: carriers only genotype orient**

Append to `test/sql/read_pfile_filter.test`:

```sql
# --- Carriers only, genotype orient ---

query TIT
SELECT ID, genotype, IID FROM read_pfile('test/data/pfile_example', orient := 'genotype', genotype_range := {min: 1}) ORDER BY ID, IID;
----
rs1	1	HG00097
rs1	2	HG00099
rs2	1	HG00096
rs2	1	HG00097
rs2	2	HG00100
rs3	2	HG00096
rs3	1	HG00099
rs4	1	HG00099
rs4	2	HG00100

# --- Hom-alt only, genotype orient (fewer rows) ---

query TIT
SELECT ID, genotype, IID FROM read_pfile('test/data/pfile_example', orient := 'genotype', genotype_range := {min: 2, max: 2}) ORDER BY ID, IID;
----
rs1	2	HG00099
rs2	2	HG00100
rs3	2	HG00096
rs4	2	HG00100
```

- [ ] **Step 4: Build and run tests**

Run: `blq run build && blq run test`
Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add src/pfile_reader.cpp test/sql/read_pfile_filter.test
git commit -m "feat: integrate genotype_range into PfileTidyScan (genotype orient)

Genotype orient skips entire rows for out-of-range and missing genotypes
when genotype_range is active."
```

### Task 7: Add genotype_range filtering to sample orient pre-pass

**Files:**
- Modify: `src/pfile_reader.cpp:744-842` (sample orient count filter pre-pass)
- Modify: `src/pfile_reader.cpp:1017-1074` (sample orient pre-read loop)
- Test: `test/sql/read_pfile_filter.test`

The sample orient pre-pass already has a PgrGetCounts loop for count filters. We extend it to also check genotype_range. For partial-match variants, we store -9 sentinel for out-of-range values in the genotype matrix during pre-read (existing scan code already NULLs -9).

- [ ] **Step 1: Extend the condition for entering the pre-pass block**

Change line 744 from:
```cpp
        if (bind_data->count_filter.HasFilter()) {
```
to:
```cpp
        if (bind_data->count_filter.HasFilter() || bind_data->genotype_filter.active) {
```

- [ ] **Step 2: Extend the pre-pass PgrGetCounts loop with genotype_range check**

In the pre-pass loop (lines 813-830), after the count filter check and before `filtered_indices.push_back`, add the genotype range check. We also need to track per-variant `all_pass` flags for use during pre-read.

Add a vector to store all_pass flags before the loop:
```cpp
            vector<bool> genotype_range_all_pass; // per-variant all_pass flags
```

Extend the loop body:
```cpp
                if (bind_data->count_filter.HasFilter() &&
                    !VariantPassesCountFilter(bind_data->count_filter, genocounts, cf_sc)) {
                    continue;
                }

                bool all_pass = true;
                if (bind_data->genotype_filter.active) {
                    auto gr = CheckGenotypeRange(bind_data->genotype_filter.range, genocounts);
                    if (!gr.any_pass) {
                        continue;
                    }
                    all_pass = gr.all_pass;
                }

                filtered_indices.push_back(vidx);
                genotype_range_all_pass.push_back(all_pass);
```

- [ ] **Step 3: Apply per-element filtering during pre-read**

In the pre-read loop (starting around line 1023), after reading genotypes but before storing in genotype_matrix, apply the genotype_range filter for partial-match variants.

After the `GenoarrToBytesMinus9` call (around line 1071-1072), add:
```cpp
                // Apply genotype_range per-element filtering
                if (bind_data->genotype_filter.active && !genotype_range_all_pass[ev]) {
                    for (uint32_t s = 0; s < output_sample_ct; s++) {
                        int8_t geno = tmp_bytes[s];
                        if (geno != -9 &&
                            !bind_data->genotype_filter.range.Passes(static_cast<double>(geno))) {
                            tmp_bytes[s] = -9; // will be NULLed at emit time
                        }
                    }
                }
```

For phased pre-read (around line 1059-1061), apply similarly on diploid sum:
```cpp
                if (bind_data->genotype_filter.active && !genotype_range_all_pass[ev]) {
                    for (uint32_t s = 0; s < output_sample_ct; s++) {
                        int8_t a1 = preread_phased_pairs[s * 2];
                        int8_t a2 = preread_phased_pairs[s * 2 + 1];
                        if (a1 != -9 &&
                            !bind_data->genotype_filter.range.Passes(static_cast<double>(a1 + a2))) {
                            preread_phased_pairs[s * 2] = -9;
                            preread_phased_pairs[s * 2 + 1] = -9;
                        }
                    }
                }
```

**Scoping note:** The `genotype_range_all_pass` vector must be declared **before** the pre-pass `if` block (at line 742, immediately after the `} else if (bind_data->orient_mode == OrientMode::SAMPLE) {` branch entry) so it's accessible in both the pre-pass loop and the pre-read loop:

```cpp
        // Declare at sample orient scope level (line 742, before the pre-pass if-block)
        vector<bool> genotype_range_all_pass;
```

It's populated in the pre-pass loop (Step 2) and consumed in the pre-read loop (Step 3). When no pre-pass is needed, the vector stays empty and the pre-read loop won't reference it (guarded by `bind_data->genotype_filter.active`).

- [ ] **Step 4: Write test: sample orient with genotype_range**

Append to `test/sql/read_pfile_filter.test`:

```sql
# --- Sample orient with genotype_range ---

query I
SELECT len(genotypes) FROM read_pfile('test/data/pfile_example', orient := 'sample', genotype_range := {min: 1}) LIMIT 1;
----
4

# Verify NULLs in correct positions for sample orient
query TT
SELECT IID, genotypes FROM read_pfile('test/data/pfile_example', orient := 'sample', genotype_range := {min: 1}) ORDER BY IID;
----
HG00096	[NULL, 1, 2, NULL]
HG00097	[1, 1, NULL, NULL]
HG00099	[2, NULL, 1, 1]
HG00100	[NULL, 2, NULL, 2]
```

- [ ] **Step 5: Build and run tests**

Run: `blq run build && blq run test`
Expected: All tests pass

- [ ] **Step 6: Commit**

```bash
git add src/pfile_reader.cpp test/sql/read_pfile_filter.test
git commit -m "feat: integrate genotype_range into sample orient pre-pass

Any-pass check filters variant from matrix entirely.
Partial-match variants get -9 sentinel for out-of-range elements.
Existing -9 → NULL handling at emit time produces correct output."
```

---

## Chunk 5: Combined Filters, Validation, and Pre-decompression Tests

### Task 8: Combined genotype_range + af_range/ac_range tests

**Files:**
- Test: `test/sql/read_pfile_filter.test`
- Test: `test/sql/read_pgen_filter.test`

- [ ] **Step 1: Write combined filter tests for read_pfile**

Append to `test/sql/read_pfile_filter.test`:

```sql
# --- Combined af_range + genotype_range ---

query TT
SELECT ID, genotypes FROM read_pfile('test/data/pfile_example', af_range := {max: 0.4}, genotype_range := {min: 1}) ORDER BY ID;
----
rs4	[NULL, NULL, 1, 2]

# --- Columns mode with genotype_range ---

query TIIIII
SELECT ID, HG00096, HG00097, HG00099, HG00100 FROM read_pfile('test/data/pfile_example', genotypes := 'columns', genotype_range := {min: 1}) ORDER BY ID;
----
rs1	NULL	1	2	NULL
rs2	1	1	NULL	2
rs3	2	NULL	1	NULL
rs4	NULL	NULL	1	2

# --- Phased mode with genotype_range (spec test case 9) ---
# Diploid sum used for range check; pair NULLed as unit
# rs1: [0,1,2,NULL] → phased: [[0,0],[0,1],[1,1],NULL]
# genotype_range {min: 1}: sum>=1 keeps het(sum=1) and hom_alt(sum=2), excludes hom_ref(sum=0) and missing

query TT
SELECT ID, genotypes FROM read_pfile('test/data/pfile_example', phased := true, genotype_range := {min: 1}) ORDER BY ID;
----
rs1	[NULL, [0, 1], [1, 1], NULL]
rs2	[[0, 1], [0, 1], NULL, [1, 1]]
rs3	[[1, 1], NULL, [0, 1], NULL]
rs4	[NULL, NULL, [0, 1], [1, 1]]
```

- [ ] **Step 2: Write combined filter tests for read_pgen**

Append to `test/sql/read_pgen_filter.test`:

```sql
# --- Combined af_range + genotype_range ---

query TT
SELECT ID, genotypes FROM read_pgen('test/data/pgen_example.pgen', af_range := {max: 0.4}, genotype_range := {min: 1}) ORDER BY ID;
----
rs4	[NULL, NULL, 1, 2]
```

- [ ] **Step 3: Run tests**

Run: `blq run test`
Expected: All tests pass

- [ ] **Step 4: Commit**

```bash
git add test/sql/read_pfile_filter.test test/sql/read_pgen_filter.test
git commit -m "test: add combined genotype_range + af_range/ac_range tests"
```

### Task 9: Pre-decompression optimization tests

**Files:**
- Test: `test/sql/read_pfile_filter.test`
- Test: `test/sql/read_pgen_filter.test`

- [ ] **Step 1: Write pre-decompression optimization tests**

Append to `test/sql/read_pfile_filter.test`:

```sql
# --- Pre-decompression: all_missing fixture, any genotype_range → variant skipped ---

query I
SELECT COUNT(*) FROM read_pfile('test/data/all_missing', genotype_range := {min: 0, max: 2});
----
0

# --- Pre-decompression: variant skipping verification ---
# rs2 has genotypes [1, 1, 0, 2] — no missing, all values present
# genotype_range {min: 2, max: 2} means only hom_alt passes
# rs1 has hom_alt (sample 3), rs2 has hom_alt (sample 4), rs3 has hom_alt (sample 1), rs4 has hom_alt (sample 4)
# All 4 variants have at least one hom_alt, so all survive

# Use genotype orient to verify row-level skipping
query I
SELECT COUNT(*) FROM read_pfile('test/data/pfile_example', orient := 'genotype', genotype_range := {min: 2, max: 2});
----
4
```

Append to `test/sql/read_pgen_filter.test`:

```sql
# --- Pre-decompression: all_missing fixture ---

query I
SELECT COUNT(*) FROM read_pgen('test/data/all_missing.pgen', genotype_range := {min: 0, max: 2});
----
0
```

- [ ] **Step 2: Run tests**

Run: `blq run test`
Expected: All tests pass

- [ ] **Step 3: Commit**

```bash
git add test/sql/read_pfile_filter.test test/sql/read_pgen_filter.test
git commit -m "test: add pre-decompression optimization tests for genotype_range"
```

### Task 10: Validation error tests

**Files:**
- Test: `test/sql/read_pfile_filter.test`
- Test: `test/sql/read_pgen_filter.test`

- [ ] **Step 1: Write validation error tests**

Append to `test/sql/read_pfile_filter.test`:

```sql
# --- Validation errors for genotype_range ---

# min > max
statement error
SELECT * FROM read_pfile('test/data/pfile_example', genotype_range := {min: 2, max: 0});
----
min (2) > max (0)

# Value out of range [0, 2]
statement error
SELECT * FROM read_pfile('test/data/pfile_example', genotype_range := {max: 3});
----
out of range

# Negative value
statement error
SELECT * FROM read_pfile('test/data/pfile_example', genotype_range := {min: -1});
----
out of range

# Incompatible with dosages
statement error
SELECT * FROM read_pfile('test/data/pfile_example', genotype_range := {min: 1}, dosages := true);
----
genotype_range is incompatible with dosages
```

Append to `test/sql/read_pgen_filter.test`:

```sql
# --- Validation errors for genotype_range ---

statement error
SELECT * FROM read_pgen('test/data/pgen_example.pgen', genotype_range := {min: 2, max: 0});
----
min (2) > max (0)

statement error
SELECT * FROM read_pgen('test/data/pgen_example.pgen', genotype_range := {max: 3});
----
out of range

statement error
SELECT * FROM read_pgen('test/data/pgen_example.pgen', genotype_range := {min: 1}, dosages := true);
----
genotype_range is incompatible with dosages
```

- [ ] **Step 2: Run all tests**

Run: `blq run test`
Expected: All tests pass

- [ ] **Step 3: Commit**

```bash
git add test/sql/read_pfile_filter.test test/sql/read_pgen_filter.test
git commit -m "test: add genotype_range validation error tests"
```

---

## Chunk 6: CLAUDE.md Update and Final Verification

### Task 11: Update CLAUDE.md conventions

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Add genotype_range convention to CLAUDE.md**

In the Conventions section, after the existing af_range/ac_range bullets, add:

```
- `genotype_range := {min: 1, max: 2}` filters individual genotype values; uses PgrGetCounts for pre-decompression optimization (skip variant if no values in range); per-element NULLing for partial matches
- `genotype_range` incompatible with `dosages := true`; domain [0, 2]
- `GenotypeRangeResult` from `CheckGenotypeRange`: `any_pass` = skip variant entirely when false; `all_pass` = skip per-element check when true (existing -9 handling suffices)
```

- [ ] **Step 2: Run full test suite**

Run: `blq run build && blq run test`
Expected: All tests pass (existing + new)

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: add genotype_range conventions to CLAUDE.md"
```
