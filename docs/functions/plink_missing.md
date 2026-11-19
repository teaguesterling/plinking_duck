# plink_missing

Compute per-variant or per-sample missingness rates.

## Synopsis

```sql
plink_missing([path VARCHAR] [, catalog := ..., pvar := ..., psam := ..., samples := ...,
              region := ..., selection := ..., mode := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `path` | `VARCHAR` | *(required unless `catalog :=`)* | Path to the `.pgen` file, or to a [catalog](#catalogs) |
| `catalog` | `VARCHAR` or `BOOLEAN` | Discovery | A `.pfile.json` [catalog](#catalogs) path; `false` turns discovery off; `true` requires a catalog |
| `pvar` | `VARCHAR` | Auto-discovered | Path to `.pvar` or `.bim` file |
| `psam` | `VARCHAR` | Auto-discovered | Path to `.psam` or `.fam` file |
| `samples` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All | Subset to specific samples |
| `region` | `VARCHAR` | All | Filter to genomic region (`chr:start-end`) |
| `selection` | `STRUCT` | None | Composable per-axis [selection](../common-parameters.md#selection) struct (`samples`/`region`/`variant_range`) |
| `mode` | `VARCHAR` | `'variant'` | `'variant'` for per-variant or `'sample'` for per-sample |

See [Common Parameters](../common-parameters.md) for details on `pvar`, `psam`, `samples`, `region`, and the composable `selection` struct.

## Output Columns

### Variant Mode (`mode := 'variant'`)

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | `VARCHAR` | Chromosome |
| `POS` | `INTEGER` | Base-pair position |
| `ID` | `VARCHAR` | Variant identifier |
| `REF` | `VARCHAR` | Reference allele |
| `ALT` | `VARCHAR` | Alternate allele |
| `MISSING_CT` | `INTEGER` | Number of missing genotypes |
| `OBS_CT` | `INTEGER` | Number of non-missing genotypes |
| `F_MISS` | `DOUBLE` | Fraction missing (`MISSING_CT / total_samples`) |

### Sample Mode (`mode := 'sample'`)

| Column | Type | Description |
|--------|------|-------------|
| `FID` | `VARCHAR` | Family ID (NULL if unavailable) |
| `IID` | `VARCHAR` | Individual ID |
| `MISSING_CT` | `INTEGER` | Number of missing genotypes across variants |
| `OBS_CT` | `INTEGER` | Number of non-missing genotypes across variants |
| `F_MISS` | `DOUBLE` | Fraction missing (`MISSING_CT / total_variants`) |

## Description

`plink_missing` computes missingness statistics using pgenlib's `PgrGetMissingness` fast path, which extracts missingness without full genotype decompression.

### Variant Mode

In variant mode (the default), each row represents one variant with its missingness count across all (or selected) samples. This mode supports multi-threaded parallel scanning.

### Sample Mode

In sample mode, the function scans all variants and accumulates per-sample missing counts, then outputs one row per sample. Sample mode:

- **Requires** a `.psam` or `.fam` file (for IID output)
- Is **single-threaded** (scans all variants, then emits sample rows)
- Supports `region` filtering (only counts missingness for variants in the region)

### Projection Pushdown

If only metadata columns are selected (variant identifiers in variant mode, FID/IID in sample mode), the missingness computation is skipped.

## Examples

```sql
-- Per-variant missingness
SELECT * FROM plink_missing('data/example.pgen');
```

```sql
-- Per-sample missingness
SELECT * FROM plink_missing('data/example.pgen', mode := 'sample');
```

```sql
-- Find variants with high missingness
SELECT ID, F_MISS
FROM plink_missing('data/example.pgen')
WHERE F_MISS > 0.05;
```

```sql
-- Find samples with high missingness
SELECT IID, F_MISS
FROM plink_missing('data/example.pgen', mode := 'sample')
WHERE F_MISS > 0.1;
```

```sql
-- Missingness in a specific region
SELECT ID, F_MISS
FROM plink_missing('data/example.pgen',
    region := '1:1-50000000');
```

## Catalogs

`plink_missing` reads a pfile catalog (`<prefix>.pfile.json`, v1 or v2) as one cohort. All catalog kinds work:
variant-sharded, sample-sharded and matrix. The rows are the rows of the same call on the un-sharded
fileset.

The inputs name a catalog as they do for [`read_pfile`](read_pfile.md):

- a positional `'<x>.pfile.json'`;
- a positional local prefix `'<p>'` (not ending in `.pgen`) when `<p>.pfile.json` exists;
- `catalog := '<x>.pfile.json'`, with no positional argument.

An explicit `.pgen` path is always a file. `catalog := false` turns discovery off, and
`catalog := true` requires a catalog at the prefix.

- Variant blocks concatenate in catalog order.
- Counts add across sample blocks. In sample mode, the sample rows concatenate across sample blocks, and each sample's counts add
  across variant blocks.
- SEX comes from each sample block's `.psam`. The writer makes IID-only `.psam` files, so add SEX
  to them before you read chrX, chrY or chrMT. Without SEX in every block, the cohort has no SEX.
- `samples` takes IIDs, or 0-based indices over the whole cohort. `region` applies to each
  variant block. `selection := {variant_range: [lo, hi)}` takes cohort-wide ordinals.
- `mode` works as for one fileset.
- `pvar` and `psam` are errors with a catalog.

```sql
SELECT * FROM plink_missing('cohort.pfile.json');
SELECT * FROM plink_missing(catalog := 'data/cohort.pfile.json');
```

## Sex chromosomes and ploidy

Missingness counts **samples, not alleles**, so chrX and chrMT need no ploidy
adjustment -- a hemizygous sample either has a call or it does not.

**chrY is the exception**, since females have no chrY at all. Matching plink2:

- *Variant mode:* the denominator is the **male count** and only males' calls are
  counted, so a fully genotyped chrY variant reports `F_MISS` 0, not 0.5.
- *Sample mode:* chrY variants leave **both** the numerator and the denominator
  for females and unknown-sex samples, so their `OBS_CT` is lower than males' by
  the number of chrY variants in range.

Unknown-sex samples are treated as not carrying chrY. Without a **SEX** column,
chrY rows are `NULL` in variant mode and sample mode spanning chrY is an error.

See [Sex-chromosome handling](https://github.com/teaguesterling/plinking_duck#sex-chromosome-handling).

## See Also

- [plink_freq](plink_freq.md) -- allele frequencies (includes MISSING_CT with `counts := true`)
- [plink_hardy](plink_hardy.md) -- HWE test
- [Quality Control Guide](../guides/quality-control.md) -- missingness thresholds in QC
