# read_pgen

Read PLINK 2 `.pgen` binary genotype files.

## Synopsis

```sql
read_pgen(path VARCHAR [, pvar := ..., psam := ..., samples := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `path` | `VARCHAR` | *(required)* | Path to the `.pgen` file |
| `pvar` | `VARCHAR` | Auto-discovered | Path to `.pvar` or `.bim` companion file |
| `psam` | `VARCHAR` | Auto-discovered | Path to `.psam` or `.fam` companion file |
| `samples` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All samples | Subset to specific samples |

See [Common Parameters](../common-parameters.md) for details on `pvar`, `psam`, and `samples`.

### Reserved Parameters (Not Yet Implemented)

| Name | Type | Description |
|------|------|-------------|
| `dosages` | `BOOLEAN` | Include dosage data (raises error if `true`) |
| `phased` | `BOOLEAN` | Include phasing data (raises error if `true`) |

## Output Columns

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | `VARCHAR` | Chromosome (from `.pvar`/`.bim`) |
| `POS` | `INTEGER` | Base-pair position (from `.pvar`/`.bim`) |
| `ID` | `VARCHAR` | Variant identifier (from `.pvar`/`.bim`) |
| `REF` | `VARCHAR` | Reference allele (from `.pvar`/`.bim`) |
| `ALT` | `VARCHAR` | Alternate allele (from `.pvar`/`.bim`) |
| `genotypes` | `LIST(TINYINT)` | Genotype calls, one element per sample |

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

## See Also

- [read_pfile](read_pfile.md) -- complete fileset reader with tidy mode
- [read_pvar](read_pvar.md) -- read variant metadata only
- [read_psam](read_psam.md) -- read sample metadata only
- [Genotype Encoding](../genotype-encoding.md) -- encoding reference
