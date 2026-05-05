# read_plink_vcf() — Fast Genotype-Focused VCF Reader

**Issue:** #11 — Faster BCF/VCF reader using plink-ng's optimized genotype parsing  
**Date:** 2026-05-05  
**Status:** Design

## Motivation

DuckHTS (`read_bcf()`) provides full-featured VCF/BCF reading via htslib — all INFO fields, FORMAT columns, tidy mode, tabix region queries. But when users only need genotypes in plink-style encoding (0/1/2/missing as TINYINT), plink-ng's VCF parser is significantly faster because it:

1. Skips INFO/FORMAT field parsing beyond GT (and optionally GQ/DP for quality filtering)
2. Packs genotypes directly into 2-bit genovec format (SIMD-friendly)
3. Handles the common biallelic case with tight inner loops

`read_plink_vcf()` exposes this fast path: read a VCF file and output genotypes in the same format as `read_pfile()`, reusing plink-ng's battle-tested GT parsing helpers.

## Scope

### In scope (v1)
- Biallelic variants
- Plain VCF (.vcf) and compressed (.vcf.gz, .vcf.bgz) input via DuckDB FileSystem
- Variant orient mode (one row per variant)
- All genotype output modes: ARRAY, LIST, columns, phased
- Region filtering (linear scan, no index required)
- GQ/DP quality filtering (plink-ng's VcfCheckQuals)
- Projection pushdown (skip genotype parsing when genotypes column not referenced)

### Out of scope (future work)
- Multiallelic variants (skip with warning or emit NULL genotypes)
- BCF binary format (requires separate header parsing; separate follow-up)
- Genotype orient and sample orient modes
- Dosages (DS/GP FORMAT fields)
- Tabix/CSI index-based region queries
- Sample subsetting at parse time
- Count filters (af_range, ac_range, genotype_range) — requires full genotype decode anyway

## Function Signature

```sql
-- Basic usage
SELECT * FROM read_plink_vcf('input.vcf.gz');

-- With phased output
SELECT * FROM read_plink_vcf('variants.vcf', phased := true);

-- Region filter
SELECT * FROM read_plink_vcf('genome.vcf.gz', region := 'chr1:1000000-2000000');

-- Columns mode
SELECT * FROM read_plink_vcf('input.vcf.gz', genotypes := 'columns');

-- Quality filtering
SELECT * FROM read_plink_vcf('input.vcf.gz', min_gq := 20, min_dp := 10);
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| path | VARCHAR | required | Path to .vcf, .vcf.gz, or .vcf.bgz file |
| genotypes | VARCHAR | 'auto' | Output mode: 'array', 'list', 'columns', 'auto' |
| phased | BOOLEAN | false | Output haplotype pairs instead of collapsed genotypes |
| region | VARCHAR | NULL | Region filter as 'chrom:start-end' or 'chrom' |
| min_gq | INTEGER | -1 | Minimum genotype quality (GQ); -1 = no filter |
| min_dp | INTEGER | -1 | Minimum read depth (DP); -1 = no filter |
| max_dp | INTEGER | -1 | Maximum read depth (DP); -1 = no filter |
| halfcall | VARCHAR | 'missing' | Half-call handling: 'missing', 'reference', 'haploid', 'error' |

### Output Schema (variant orient)

| Column | Type | Description |
|--------|------|-------------|
| chrom | VARCHAR | Chromosome (VCF CHROM field) |
| pos | INTEGER | 1-based position (VCF POS field) |
| id | VARCHAR | Variant ID (VCF ID field; NULL if '.') |
| ref | VARCHAR | Reference allele |
| alt | VARCHAR | Alternate allele (first ALT only for biallelic) |
| genotypes | ARRAY(TINYINT, N) | Genotypes: 0=hom_ref, 1=het, 2=hom_alt, NULL=missing |

When `phased := true`: genotypes column becomes `ARRAY(ARRAY(TINYINT, 2), N)` with `[allele1, allele2]` pairs.

When `genotypes := 'columns'`: genotypes column replaced by N individual TINYINT columns named by sample IID.

When `genotypes := 'list'`: genotypes column is `LIST(TINYINT)`.

## Architecture

### Data Flow

```
VCF file (.vcf / .vcf.gz)
    │
    ▼ DuckDB FileSystem (handles gzip/bgzf transparently)
    │
    ▼ Header parsing (extract sample names, sample_ct)
    │
    ▼ Per-line: split fixed columns (CHROM/POS/ID/REF/ALT/QUAL/FILTER/INFO)
    │
    ▼ Per-line: locate GT in FORMAT, position linebuf_iter at first sample
    │
    ▼ plink-ng VcfConvert*BiallelicLine() → genovec (2-bit packed)
    │
    ▼ GenoarrToBytesMinus9 / UnpackPhasedGenotypes → intermediate arrays
    │
    ▼ Shared genotype output formatting → DuckDB output vectors
```

### New Files

#### `src/vcf_reader.cpp` / `src/include/vcf_reader.hpp`

Table function implementation following DuckDB lifecycle: Bind → InitGlobal → InitLocal → Scan.

**Bind phase:**
- Open VCF file via DuckDB FileSystem
- Parse header lines (## metadata — scan for FORMAT definitions)
- Parse #CHROM line to extract sample names and sample_ct
- Validate: confirm GT exists in at least one variant (scan first data line FORMAT)
- Resolve genotype mode, set up output schema
- Store variant metadata column types

**InitGlobal:**
- Record file offset after header (start of data lines)
- Set up atomic variant counter for parallel scan (future)

**InitLocal:**
- Allocate genovec buffer: `NypCtToAlignedWordCt(sample_ct)` words
- Allocate phasepresent/phaseinfo buffers if phased mode
- Allocate genotype_bytes / phased_pairs intermediate arrays
- Set up VcfImportBaseContext struct

**Scan:**
- Read lines from file (batch of STANDARD_VECTOR_SIZE)
- For each line:
  1. Parse fixed columns (tab-split first 9 fields)
  2. Check region filter (chrom + pos range)
  3. Skip multiallelic (ALT contains comma) → emit NULL genotypes or skip row
  4. Locate GT field position in FORMAT column
  5. Set up qual_field_skips if GQ/DP filtering enabled
  6. Call `VcfConvertUnphasedBiallelicLine` or `VcfConvertPhasedBiallelicLine`
  7. Convert genovec → output via shared formatting utilities

#### `src/vcf_genotype_parse.cpp` / `src/include/vcf_genotype_parse.hpp`

Thin compilation unit that compiles plink-ng's VCF parsing helpers in our build context.

**Included/wrapped functions from `plink2_import.cc`:**
- `VcfConvertUnphasedBiallelicLine`
- `VcfConvertPhasedBiallelicLine`
- `VcfCheckQuals`
- `VcfScanShortallelicLine` (for detecting if phase info present)

**Approach:** Either:
1. Copy these ~300 lines into `vcf_genotype_parse.cpp` with attribution (they depend only on header-inlines from plink2_base.h and plink2_string.h which we already compile), or
2. Include `plink2_import.cc` with preprocessor guards to only compile the needed functions

Option 1 is preferred: these functions are self-contained, stable, and the alternative would pull in `plink2_import.cc`'s 18K lines with all their dependencies (pgenlib_write, threading, bigstack).

**Exposed interface:**

```cpp
namespace duckdb {

// Wraps plink2 VcfImportBaseContext with a simpler construction API
struct VcfParseContext {
    uint32_t sample_ct;
    VcfHalfCall halfcall_mode;   // from plink2
    int32_t min_gq;              // -1 = disabled
    int32_t min_dp;              // -1 = disabled
    int32_t max_dp;              // -1 = disabled

    // Internal: populated per-line from FORMAT field
    plink2::VcfImportBaseContext vibc;

    void Initialize();
    void SetupQualFields(uint32_t gt_field_idx, uint32_t gq_field_idx, uint32_t dp_field_idx);
};

enum class VcfGenoParseResult : uint8_t {
    OK,
    MISSING_TOKENS,
    INVALID_GT,
    HALFCALL_ERROR,
    POLYPLOID_ERROR
};

// Parse unphased biallelic GT fields from FORMAT portion of VCF line
VcfGenoParseResult ParseUnphasedBiallelicGT(
    const VcfParseContext &ctx,
    const char *format_start,     // points to first sample's data after FORMAT\t
    uintptr_t *genovec);          // output: 2-bit packed genotypes

// Parse phased biallelic GT fields
VcfGenoParseResult ParsePhasedBiallelicGT(
    const VcfParseContext &ctx,
    const char *format_start,
    uintptr_t *genovec,
    uintptr_t *phasepresent,
    uintptr_t *phaseinfo);

} // namespace duckdb
```

### Modified Files

#### `src/plink_common.cpp` / `src/include/plink_common.hpp`

Extract genotype output formatting into shared utilities callable by both `pfile_reader` and `vcf_reader`:

```cpp
// Write genotype_bytes into a DuckDB output vector for one row
void FillGenotypeOutputVector(
    Vector &vec,
    idx_t row_idx,
    GenotypeMode mode,
    uint32_t output_sample_ct,
    const int8_t *genotype_bytes,       // for unphased
    const int8_t *phased_pairs,         // for phased (may be nullptr)
    const double *dosage_doubles,       // for dosages (may be nullptr)
    bool include_phased,
    bool include_dosages,
    const GenotypeFilter *filter);      // optional genotype_range filter
```

This consolidates the ~150 lines of ARRAY/LIST/phased/dosage output logic currently duplicated between the variant-orient and sample-orient scan paths in `pfile_reader.cpp`.

#### `src/pfile_reader.cpp`

Refactor the genotype output section (lines ~1900-2040) to call `FillGenotypeOutputVector` instead of inline logic. Behavior is unchanged.

#### `src/plinking_duck_extension.cpp`

Register `read_plink_vcf` table function.

#### `CMakeLists.txt`

Add `vcf_reader.cpp` and `vcf_genotype_parse.cpp` to the source list. No new third-party dependencies — the VCF parse helpers only need headers already in the include path.

## Multiallelic Handling (v1)

When a VCF line has multiple ALT alleles (comma in ALT field):
- **Default behavior:** Skip the variant entirely (do not emit a row)
- Log a warning with count of skipped multiallelic variants at end of scan
- Future: support via plink-ng's `VcfConvertUnphasedMultiallelicLine` with appropriate output encoding

## Compressed File Support

DuckDB's FileSystem handles gzip-compressed files transparently when opened. For `.vcf.gz` and `.vcf.bgz` files:
- Use `FileSystem::OpenFile` with standard flags — DuckDB auto-detects compression
- Sequential scan only (no seeking) — compatible with gzip streams
- bgzf is a valid gzip stream, so bgzf files work without special handling

## Error Handling

| Condition | Behavior |
|-----------|----------|
| File not found / unreadable | IOException at bind |
| No #CHROM header line | IOException at bind |
| No samples in header | IOException at bind (unless future no_samples_ok) |
| Multiallelic variant | Skip row (log warning) |
| Malformed GT field | IOException with line number |
| Polyploid GT (e.g., 0/0/1) | Set missing (default) or error per halfcall mode |
| Half-call (e.g., 0/.) | Per halfcall parameter: missing/reference/haploid/error |
| GQ/DP below threshold | Set genotype to missing for that sample |

## Performance Characteristics

Compared to `read_bcf()` (DuckHTS/htslib):
- **Faster** for genotype-only queries: skips INFO parsing, FORMAT field iteration beyond GT/GQ/DP
- **Lower memory**: no htslib record allocation, direct 2-bit packing
- **Comparable** for variant metadata: both must parse the tab-delimited fixed columns

Compared to `read_pfile()` (native .pgen):
- **Slower**: VCF is text, .pgen is binary with random access. VCF requires full sequential scan.
- **No parallelism** in v1: single-threaded sequential line reading (future: split file into chunks for bgzf)

Expected use case: "I have a VCF from an imputation server / collaborator and want genotypes in plink format without converting to .pgen first."

## Testing

- `test/sql/vcf_reader_positive.test` — basic reads, phased, columns, list modes
- `test/sql/vcf_reader_negative.test` — missing files, malformed GT, no header
- `test/sql/vcf_reader_quality.test` — GQ/DP filtering
- `test/sql/vcf_reader_region.test` — region filtering
- Test data: generate small VCF fixtures with known genotypes (4 variants × 4 samples, matching pgen_example data for cross-validation)

## Future Extensions

1. **Multiallelic support** — Use `VcfConvertUnphasedMultiallelicLine` / `VcfConvertPhasedMultiallelicLine`
2. **BCF binary format** — Requires BCF header parsing (plink-ng has `BcfToPgen` helpers)
3. **Dosages** — Parse DS/GP FORMAT fields via `VcfConvertPhasedBiallelicDosageLine`
4. **Genotype/sample orient modes** — Reuse shared orient machinery from pfile_reader
5. **Parallel scan** — For bgzf files, split into block-aligned chunks for multi-threaded parsing
6. **Count filters** — af_range/ac_range/genotype_range (would need genotype counting per variant)
7. **Sample subsetting** — Parse only requested samples' GT fields
8. **Index-based region queries** — Use .tbi/.csi for random access in bgzf files
