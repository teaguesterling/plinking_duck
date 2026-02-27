# PlinkingDuck

A DuckDB extension for reading [PLINK 2](https://www.cog-genomics.org/plink/2.0/) genomics file formats directly in SQL.

PlinkingDuck brings PLINK genotype, variant, and sample data into DuckDB, letting you query genomics datasets with standard SQL instead of format-specific tools.

## Functions

### `read_pvar(path)`

Read PLINK2 `.pvar` (variant information) or legacy `.bim` files into SQL-queryable tables.
Format is auto-detected: files with a `#CHROM` header line are parsed as `.pvar`,
files without a header are parsed as legacy `.bim` (6-column, tab-delimited).

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
the .bim file order of CHROM, ID, CM, POS, ALT, REF).

Dot values (`.`) in any column are returned as `NULL`. Projection pushdown is
supported — only columns referenced in the query are parsed.

### `read_psam(path)`

Read PLINK 2 `.psam` (sample information) or PLINK 1 `.fam` files as DuckDB tables.

```sql
-- Read a .psam file
SELECT * FROM read_psam('samples.psam');

-- Count samples per family
SELECT FID, COUNT(*) AS n FROM read_psam('cohort.psam') GROUP BY FID;

-- Filter by sex
SELECT IID FROM read_psam('cohort.psam') WHERE SEX = 2;

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

**Projection pushdown:** Genotype decoding is skipped when genotype columns
are not referenced in the query, making metadata-only queries fast.

### `plink_hardy(path [, pvar, psam, samples, region, midp])`

Compute per-variant Hardy-Weinberg equilibrium exact test p-values.
Uses pgenlib's fast counting path (no genotype decompression).

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

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pvar` | VARCHAR | auto-discovered | Explicit `.pvar`/`.bim` path |
| `psam` | VARCHAR | auto-discovered | Explicit `.psam`/`.fam` path |
| `samples` | LIST(VARCHAR) or LIST(INTEGER) | all | Filter to specific sample IDs |
| `region` | VARCHAR | all | Filter to genomic region (`chr:start-end`) |
| `midp` | BOOLEAN | false | Use mid-p correction for exact test |

## Genotype Encoding

All genotype data uses `TINYINT` values with the following encoding:

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
| `.pgen` | Binary genotypes | Compressed genotype matrix (variants × samples) |
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
