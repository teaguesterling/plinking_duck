# plink_missing

Compute per-variant or per-sample missingness rates.

## Synopsis

```sql
plink_missing(path VARCHAR [, pvar := ..., psam := ..., samples := ...,
              region := ..., mode := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `path` | `VARCHAR` | *(required)* | Path to the `.pgen` file |
| `pvar` | `VARCHAR` | Auto-discovered | Path to `.pvar` or `.bim` file |
| `psam` | `VARCHAR` | Auto-discovered | Path to `.psam` or `.fam` file |
| `samples` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All | Subset to specific samples |
| `region` | `VARCHAR` | All | Filter to genomic region (`chr:start-end`) |
| `mode` | `VARCHAR` | `'variant'` | `'variant'` for per-variant or `'sample'` for per-sample |

See [Common Parameters](../common-parameters.md) for details on `pvar`, `psam`, `samples`, and `region`.

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

## See Also

- [plink_freq](plink_freq.md) -- allele frequencies (includes MISSING_CT with `counts := true`)
- [plink_hardy](plink_hardy.md) -- HWE test
- [Quality Control Guide](../guides/quality-control.md) -- missingness thresholds in QC
