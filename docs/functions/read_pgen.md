# read_pgen

Read PLINK 2 `.pgen` binary genotype files.

## Synopsis

```sql
read_pgen(path VARCHAR [, pvar := ..., psam := ..., samples := ...,
          genotypes := ..., phased := ..., dosages := ...,
          af_range := ..., ac_range := ..., genotype_range := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `path` | `VARCHAR` | *(required)* | Path to the `.pgen` file |
| `pvar` | `VARCHAR` | Auto-discovered | Path to `.pvar` or `.bim` companion file |
| `psam` | `VARCHAR` | Auto-discovered | Path to `.psam` or `.fam` companion file |
| `samples` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All samples | Subset to specific samples |
| `genotypes` | `VARCHAR` | `'array'` | Genotype output format: `'array'`, `'list'`, or `'columns'` |
| `phased` | `BOOLEAN` | `false` | Output phased haplotype pairs |
| `dosages` | `BOOLEAN` | `false` | Output dosage values |
| `af_range` | `STRUCT(min DOUBLE, max DOUBLE)` | All | Filter variants by allele frequency |
| `ac_range` | `STRUCT(min INTEGER, max INTEGER)` | All | Filter variants by allele count |
| `genotype_range` | `STRUCT(min TINYINT, max TINYINT)` | All | Filter individual genotype values |

See [Common Parameters](../common-parameters.md) for details on `pvar`, `psam`, `samples`, and filter parameters.

## Output Columns

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | `VARCHAR` | Chromosome (from `.pvar`/`.bim`) |
| `POS` | `INTEGER` | Base-pair position (from `.pvar`/`.bim`) |
| `ID` | `VARCHAR` | Variant identifier (from `.pvar`/`.bim`) |
| `REF` | `VARCHAR` | Reference allele (from `.pvar`/`.bim`) |
| `ALT` | `VARCHAR` | Alternate allele (from `.pvar`/`.bim`) |
| `genotypes` | `ARRAY(TINYINT, N)` | Genotype calls, one element per sample |

Each element in the `genotypes` list is a [genotype value](../genotype-encoding.md): 0 (hom ref), 1 (het), 2 (hom alt), or NULL (missing).

## Description

`read_pgen` reads binary genotype data from `.pgen` files and joins it with variant metadata from a companion `.pvar` or `.bim` file. Each output row represents one variant with a list of genotype calls across all (or selected) samples.

### Companion File Discovery

The function requires a `.pvar` or `.bim` file for variant metadata. If not specified via the `pvar` parameter, it auto-discovers by replacing the `.pgen` extension.

The `.psam` file is **optional**. Without it:

- Genotype list length is determined from the `.pgen` header
- The `samples` parameter only accepts `LIST(INTEGER)` (no IID matching)

### Projection Pushdown

If the `genotypes` column is not referenced in the query, genotype decoding is skipped entirely. This makes metadata-only queries very fast:

```sql
-- Fast: no genotype decoding
SELECT CHROM, POS, ID FROM read_pgen('data.pgen');

-- Requires genotype decoding
SELECT ID, genotypes FROM read_pgen('data.pgen');
```

### Parallel Scan

Variant processing is parallelized across multiple threads. Each thread gets its own pgenlib reader instance. The degree of parallelism scales with the number of variants.

### Filter Pushdown

`af_range` and `ac_range` filter variants by allele frequency or count using fast genotype counting (no decompression needed). `genotype_range` filters individual genotype values, setting non-matching values to NULL.

See [Common Parameters](../common-parameters.md) for details on filter parameters.

## Examples

```sql
-- Basic usage (auto-discovers .pvar and .psam)
SELECT * FROM read_pgen('data/example.pgen');
```

```sql
-- Explicit companion file paths
SELECT * FROM read_pgen('data/example.pgen',
    pvar := 'other/variants.pvar',
    psam := 'other/samples.psam');
```

```sql
-- Metadata-only query (fast, skips genotype decoding)
SELECT CHROM, POS, ID FROM read_pgen('data/example.pgen');
```

```sql
-- Subset to specific samples by name
SELECT ID, genotypes
FROM read_pgen('data/example.pgen', samples := ['SAMPLE1', 'SAMPLE3']);
```

```sql
-- Subset to specific samples by 0-based index
SELECT ID, genotypes
FROM read_pgen('data/example.pgen', samples := [0, 2]);
```

```sql
-- Phased haplotype output
SELECT ID, genotypes
FROM read_pgen('data/example.pgen', phased := true);
```

```sql
-- Filter to common variants
SELECT * FROM read_pgen('data/example.pgen',
    af_range := {min: 0.01, max: 0.5});
```

```sql
-- Columns mode: one column per sample
SELECT * FROM read_pgen('data/example.pgen',
    genotypes := 'columns');
```

## See Also

- [read_pfile](read_pfile.md) -- complete fileset reader with orient modes
- [read_pvar](read_pvar.md) -- read variant metadata only
- [read_psam](read_psam.md) -- read sample metadata only
- [Genotype Encoding](../genotype-encoding.md) -- encoding and phased output reference
