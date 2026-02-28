# PlinkingDuck

A DuckDB extension for reading [PLINK 2](https://www.cog-genomics.org/plink/2.0/) genomics file formats and running common genetic analyses directly in SQL.

**[Documentation](https://plinking-duck.readthedocs.io)** · **[Getting Started](https://plinking-duck.readthedocs.io/getting-started/)** · **[Tutorial](https://plinking-duck.readthedocs.io/tutorial/)**

PlinkingDuck brings PLINK genotype, variant, and sample data into DuckDB, letting you query genomics datasets with standard SQL instead of format-specific command-line tools. Read files, filter variants, compute allele frequencies, test Hardy-Weinberg equilibrium, measure missingness, calculate linkage disequilibrium, and run polygenic scoring — all from a SQL prompt.

> **Built on the shoulders of [DuckHTS](https://github.com/RGenomicsETL/duckhts).**
> DuckHTS pioneered the idea of querying genomics file formats directly in DuckDB using htslib. It supports VCF, BCF, BAM, CRAM, FASTA, FASTQ, GTF, GFF, and tabix-indexed files — the sequencing side of the house. PlinkingDuck picks up where DuckHTS leaves off, covering the PLINK genotype file formats used in population genetics and GWAS. If you work with both sequencing and genotype data, you'll want both extensions.

## What's Included

PlinkingDuck provides **four file readers** and **five analysis functions**:

| Function | Purpose |
|----------|---------|
| [`read_pvar`](#read_pvarpath) | Read `.pvar` / `.bim` variant metadata |
| [`read_psam`](#read_psampath) | Read `.psam` / `.fam` sample metadata |
| [`read_pgen`](#read_pgenpath--pvar-psam-samples) | Read `.pgen` binary genotypes |
| [`read_pfile`](#read_pfileprefix--pgen-pvar-psam-tidy-samples-variants-region) | Unified reader with tidy mode, filtering |
| [`plink_freq`](#plink_freqpath--pvar-psam-samples-region-counts) | Per-variant allele frequencies |
| [`plink_hardy`](#plink_hardypath--pvar-psam-samples-region-midp) | Hardy-Weinberg equilibrium test |
| [`plink_missing`](#plink_missingpath--pvar-psam-samples-region-mode) | Per-variant or per-sample missingness |
| [`plink_ld`](#plink_ldpath--variant1-variant2-window_kb-r2_threshold-inter_chr) | Pairwise linkage disequilibrium |
| [`plink_score`](#plink_scorepath--weights-no_mean_imputation) | Polygenic risk scoring |

All functions support **projection pushdown** (skip expensive genotype decoding when columns aren't referenced) and **parallel scanning** (multi-threaded variant processing with atomic batch claiming).

---

## File Readers

### `read_pvar(path)`

Read PLINK 2 `.pvar` (variant information) or legacy `.bim` files. Format is auto-detected.

```sql
-- Read a .pvar file
SELECT * FROM read_pvar('path/to/file.pvar');

-- Read a legacy .bim file (auto-detected, columns normalized to .pvar order)
SELECT * FROM read_pvar('path/to/file.bim');

-- Filter variants by chromosome
SELECT CHROM, POS, ID, REF, ALT
FROM read_pvar('data/variants.pvar')
WHERE CHROM = '22';
```

**Output columns (.pvar):** CHROM, POS, ID, REF, ALT, and any optional columns
present in the header (QUAL, FILTER, INFO, CM).

**Output columns (.bim):** CHROM, POS, ID, REF, ALT, CM (normalized from
the .bim column order of CHROM, ID, CM, POS, ALT, REF).

Dot values (`.`) in any column are returned as `NULL`.

### `read_psam(path)`

Read PLINK 2 `.psam` (sample information) or PLINK 1 `.fam` files.

```sql
-- Read a .psam file
SELECT * FROM read_psam('samples.psam');

-- Count samples per family
SELECT FID, COUNT(*) AS n FROM read_psam('cohort.psam') GROUP BY FID;

-- Join sample info with phenotype data
SELECT s.IID, s.SEX, p.BMI
FROM read_psam('cohort.psam') s
JOIN 'phenotypes.csv' p ON s.IID = p.IID;

-- Read a legacy .fam file
SELECT * FROM read_psam('cohort.fam');
```

**Supported formats:**

| Format | Header | Columns |
|--------|--------|---------|
| `.psam` (`#FID` header) | `#FID  IID  SEX  [phenotypes...]` | FID + IID + dynamic columns |
| `.psam` (`#IID` header) | `#IID  SEX  [phenotypes...]` | IID + dynamic columns (no FID) |
| `.fam` (no header) | Fixed 6 columns | FID, IID, PAT, MAT, SEX, PHENO1 |

**Type mapping:**

| Column | DuckDB Type | Notes |
|--------|-------------|-------|
| SEX | INTEGER | 1=male, 2=female; `0`/`NA`/`.` → NULL |
| PAT, MAT | VARCHAR | `0`/`.`/`NA` → NULL (unknown parent) |
| All others | VARCHAR | `.`/`NA` → NULL |

### `read_pgen(path [, pvar, psam, samples])`

Read PLINK 2 `.pgen` binary genotype files. Auto-discovers companion `.pvar` and `.psam` files.

```sql
-- Basic usage (auto-discovers .pvar and .psam)
SELECT * FROM read_pgen('path/to/file.pgen');

-- Explicit companion file paths
SELECT * FROM read_pgen('file.pgen', pvar := 'other.pvar', psam := 'other.psam');

-- Only variant metadata (skips genotype decoding — fast)
SELECT chrom, pos, id FROM read_pgen('file.pgen');

-- Subset samples by 0-based index
SELECT id, genotypes FROM read_pgen('file.pgen', samples := [0, 2]);
```

**Output columns:** CHROM, POS, ID, REF, ALT (from `.pvar`), genotypes (from `.pgen`).

**Projection pushdown:** If the `genotypes` column is not referenced in the query,
genotype decoding is skipped entirely, making metadata-only queries very fast.

**Optional .psam:** The `.psam` file is optional for `read_pgen`. Without it,
genotype lists are sized from the `.pgen` header and the `samples` parameter
only accepts integer indices.

### `read_pfile(prefix [, pgen, pvar, psam, tidy, samples, variants, region])`

Read a complete PLINK 2 fileset (`.pgen` + `.pvar` + `.psam`) from a common prefix.
Supports tidy output mode, sample subsetting, variant filtering, and region filtering.

```sql
-- Read all variants (one row per variant, genotypes as list)
SELECT * FROM read_pfile('data/example');

-- Tidy mode (one row per variant x sample)
SELECT chrom, pos, iid, genotype
FROM read_pfile('data/example', tidy := true);

-- Sample subset (by name or 0-based index)
SELECT * FROM read_pfile('data/example',
    tidy := true, samples := ['SAMPLE1', 'SAMPLE3']);

-- Region filter
SELECT * FROM read_pfile('data/example', region := '22:1-50000000');

-- Variant filter
SELECT * FROM read_pfile('data/example', variants := ['rs1', 'rs2']);

-- Combined filters
SELECT iid, genotype
FROM read_pfile('data/example', tidy := true,
    region := '1:10000-50000', samples := ['SAMPLE1']);
```

**Default mode output:** CHROM, POS, ID, REF, ALT, genotypes `LIST(TINYINT)` — same as `read_pgen`.

**Tidy mode output:** CHROM, POS, ID, REF, ALT, plus all `.psam` columns (FID, IID, SEX, etc.),
plus scalar `genotype` (TINYINT). One row per variant x sample combination.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `pgen` | VARCHAR | Explicit `.pgen` path (overrides prefix) |
| `pvar` | VARCHAR | Explicit `.pvar`/`.bim` path (overrides prefix) |
| `psam` | VARCHAR | Explicit `.psam`/`.fam` path (overrides prefix) |
| `tidy` | BOOLEAN | Tidy output: one row per (variant, sample). Default: false |
| `samples` | LIST(VARCHAR) or LIST(INTEGER) | Filter to specific samples by IID or 0-based index |
| `variants` | LIST(VARCHAR) or LIST(INTEGER) | Filter to specific variants by ID or 0-based index |
| `region` | VARCHAR | Filter to genomic region (`chr:start-end` or `chr`) |

---

## Analysis Functions

These functions use pgenlib's fast counting paths (no genotype decompression), making them efficient even on large datasets. All support sample subsetting, region filtering, and parallel execution.

### `plink_freq(path [, pvar, psam, samples, region, counts])`

Compute per-variant allele frequencies.

```sql
-- Basic allele frequencies
SELECT * FROM plink_freq('data/example.pgen');

-- With genotype counts
SELECT * FROM plink_freq('data/example.pgen', counts := true);

-- MAF filter
SELECT * FROM plink_freq('data/example.pgen')
WHERE ALT_FREQ > 0.01 AND ALT_FREQ < 0.99;
```

**Output columns:** CHROM, POS, ID, REF, ALT, ALT_FREQ, OBS_CT.
With `counts := true`: also HOM_REF_CT, HET_CT, HOM_ALT_CT, MISSING_CT.

### `plink_hardy(path [, pvar, psam, samples, region, midp])`

Compute per-variant Hardy-Weinberg equilibrium exact test p-values.

```sql
-- Basic usage
SELECT * FROM plink_hardy('data/example.pgen');

-- With mid-p correction
SELECT * FROM plink_hardy('data/example.pgen', midp := true);

-- QC workflow: keep variants passing HWE filter
SELECT * FROM plink_hardy('data/example.pgen')
WHERE P_HWE > 1e-6;
```

**Output columns:** CHROM, POS, ID, REF, ALT, A1, HOM_REF_CT, HET_CT,
HOM_ALT_CT, O_HET, E_HET, P_HWE.

### `plink_missing(path [, pvar, psam, samples, region, mode])`

Compute per-variant or per-sample missingness rates.

```sql
-- Per-variant missingness
SELECT * FROM plink_missing('data/example.pgen');

-- Per-sample missingness
SELECT * FROM plink_missing('data/example.pgen', mode := 'sample');

-- QC: flag samples with >5% missingness
SELECT IID, F_MISS
FROM plink_missing('data/example.pgen', mode := 'sample')
WHERE F_MISS > 0.05;
```

**Variant mode output (default):** CHROM, POS, ID, REF, ALT, MISSING_CT, OBS_CT, F_MISS.

**Sample mode output (`mode := 'sample'`):** FID, IID, MISSING_CT, OBS_CT, F_MISS.

### `plink_ld(path [, variant1, variant2, window_kb, r2_threshold, inter_chr])`

Compute pairwise linkage disequilibrium (r², D').

```sql
-- LD between two specific variants
SELECT * FROM plink_ld('data/example.pgen',
    variant1 := 'rs123', variant2 := 'rs456');

-- Windowed LD scan (1 Mb window, R² > 0.2)
SELECT * FROM plink_ld('data/example.pgen',
    window_kb := 1000, r2_threshold := 0.2);

-- Cross-chromosome LD (rare, but sometimes needed)
SELECT * FROM plink_ld('data/example.pgen', inter_chr := true);
```

**Output columns:** CHROM_A, POS_A, ID_A, CHROM_B, POS_B, ID_B, R2, D_PRIME, OBS_CT.

**Modes:**
- **Pairwise:** Specify `variant1` and `variant2` for a single pair.
- **Windowed:** Scan all pairs within `window_kb` (default: 1000 kb), optionally filtered by `r2_threshold`.

### `plink_score(path [, weights, no_mean_imputation])`

Compute per-sample polygenic risk scores from variant weights.

```sql
-- Positional weights (one weight per variant, in file order)
SELECT * FROM plink_score('data/example.pgen',
    weights := [1.0, 0.5, -0.5, 2.0]);

-- ID-keyed weights with allele specification
SELECT * FROM plink_score('data/example.pgen',
    weights := [
      {'id': 'rs1', 'allele': 'G', 'weight': 1.0},
      {'id': 'rs2', 'allele': 'T', 'weight': 0.5}
    ]);

-- Without mean imputation for missing genotypes
SELECT * FROM plink_score('data/example.pgen',
    weights := [1.0, 0.5], no_mean_imputation := true);
```

**Output columns:** FID, IID, ALLELE_CT, DENOM, NAMED_ALLELE_DOSAGE_SUM, SCORE_SUM, SCORE_AVG.

**Missing genotype handling:** By default, missing genotypes are mean-imputed using
the average dosage across non-missing samples. Use `no_mean_imputation := true` to
skip missing samples instead.

---

## Common Parameters

These parameters are shared across analysis functions and `read_pfile`:

| Parameter | Type | Description |
|-----------|------|-------------|
| `pvar` | VARCHAR | Explicit `.pvar`/`.bim` path (overrides auto-discovery) |
| `psam` | VARCHAR | Explicit `.psam`/`.fam` path (overrides auto-discovery) |
| `samples` | LIST(VARCHAR) or LIST(INTEGER) | Filter to specific samples by IID or 0-based index |
| `region` | VARCHAR | Filter to genomic region: `'chr'`, `'chr:start-end'`, or `'chr:start-'` |

## Genotype Encoding

All genotype data uses `TINYINT` values:

| Value | Meaning |
|-------|---------|
| 0 | Homozygous reference (0/0) |
| 1 | Heterozygous (0/1) |
| 2 | Homozygous alternate (1/1) |
| NULL | Missing (./.) |

In `read_pgen` and `read_pfile` default mode, genotypes are returned as `LIST(TINYINT)`
with one element per sample. In `read_pfile` tidy mode, genotype is a scalar `TINYINT` column.

## File Format Reference

PLINK 2 uses a three-file fileset to represent genotyped samples:

| File | Contents | Description |
|------|----------|-------------|
| `.pgen` | Binary genotypes | Compressed genotype matrix (variants x samples) |
| `.pvar` | Variant metadata | Tab-delimited: CHROM, POS, ID, REF, ALT (+ optional cols) |
| `.psam` | Sample metadata | Tab-delimited: FID, IID, SEX (+ optional phenotypes) |

All three files share the same prefix (e.g., `data/cohort.pgen`, `data/cohort.pvar`,
`data/cohort.psam`). Legacy PLINK 1 equivalents (`.bed`/`.bim`/`.fam`) are partially
supported: `.bim` and `.fam` files can be read through `read_pvar` and `read_psam`
respectively.

## Building from Source

```sh
make
```

This produces:
```
./build/release/duckdb                                              # DuckDB shell with extension loaded
./build/release/test/unittest                                       # Test runner
./build/release/extension/plinking_duck/plinking_duck.duckdb_extension  # Loadable extension binary
```

## Running Tests

```sh
make test
```

SQL tests live in `./test/sql/` and test data in `./test/data/`. Tests use DuckDB's
[sqllogictest](https://duckdb.org/docs/dev/sqllogictest/intro.html) framework.

## Documentation

Full documentation is available at **[plinking-duck.readthedocs.io](https://plinking-duck.readthedocs.io)**, including:

- [Getting Started](https://plinking-duck.readthedocs.io/getting-started/) — installation and first queries
- [Tutorial](https://plinking-duck.readthedocs.io/tutorial/) — hands-on genomic analysis walkthrough
- [Function Reference](https://plinking-duck.readthedocs.io/functions/) — detailed docs for every function
- [Quality Control Guide](https://plinking-duck.readthedocs.io/guides/quality-control/) — QC workflows with SQL
- [Performance Guide](https://plinking-duck.readthedocs.io/guides/performance/) — optimization tips
