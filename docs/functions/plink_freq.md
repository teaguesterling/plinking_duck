# plink_freq

Compute per-variant allele frequencies using pgenlib's fast counting path.

## Synopsis

```sql
plink_freq(path VARCHAR [, pvar := ..., psam := ..., samples := ...,
           region := ..., counts := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `path` | `VARCHAR` | *(required)* | Path to the `.pgen` file |
| `pvar` | `VARCHAR` | Auto-discovered | Path to `.pvar` or `.bim` file |
| `psam` | `VARCHAR` | Auto-discovered | Path to `.psam` or `.fam` file |
| `samples` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All | Subset to specific samples |
| `region` | `VARCHAR` | All | Filter to genomic region (`chr:start-end`) |
| `counts` | `BOOLEAN` | `false` | Include genotype count columns |

See [Common Parameters](../common-parameters.md) for details on `pvar`, `psam`, `samples`, and `region`.

### Reserved Parameters (Not Yet Implemented)

| Name | Type | Description |
|------|------|-------------|
| `dosage` | `BOOLEAN` | Dosage-based frequency computation (raises error if `true`) |

## Output Columns

### Default Columns

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | `VARCHAR` | Chromosome |
| `POS` | `INTEGER` | Base-pair position |
| `ID` | `VARCHAR` | Variant identifier |
| `REF` | `VARCHAR` | Reference allele |
| `ALT` | `VARCHAR` | Alternate allele |
| `ALT_FREQ` | `DOUBLE` | Alternate allele frequency (NULL if all missing) |
| `OBS_CT` | `INTEGER` | Number of observed alleles (2 * non-missing samples) |

### Additional Columns with `counts := true`

| Column | Type | Description |
|--------|------|-------------|
| `HOM_REF_CT` | `INTEGER` | Homozygous reference count |
| `HET_CT` | `INTEGER` | Heterozygous count |
| `HOM_ALT_CT` | `INTEGER` | Homozygous alternate count |
| `MISSING_CT` | `INTEGER` | Missing genotype count |

## Description

`plink_freq` computes allele frequencies for each variant using pgenlib's `PgrGetCounts` fast path, which counts genotypes without full decompression. This is significantly faster than extracting individual genotypes.

### Frequency Calculation

The alternate allele frequency is computed as:

```
ALT_FREQ = (HET_CT + 2 * HOM_ALT_CT) / (2 * (HOM_REF_CT + HET_CT + HOM_ALT_CT))
```

When all genotypes for a variant are missing, `ALT_FREQ` is `NULL` and `OBS_CT` is `0`.

### Projection Pushdown

If only metadata columns (CHROM, POS, ID, REF, ALT) are selected, genotype counting is skipped entirely.

## Examples

```sql
-- Basic allele frequencies
SELECT ID, ALT_FREQ FROM plink_freq('data/example.pgen');
```

```sql
-- With genotype counts
SELECT ID, ALT_FREQ, HOM_REF_CT, HET_CT, HOM_ALT_CT
FROM plink_freq('data/example.pgen', counts := true);
```

```sql
-- Filter rare variants (MAF < 1%)
SELECT ID, ALT_FREQ
FROM plink_freq('data/example.pgen')
WHERE ALT_FREQ < 0.01 OR ALT_FREQ > 0.99;
```

```sql
-- Frequency in a sample subset
SELECT ID, ALT_FREQ
FROM plink_freq('data/example.pgen',
    samples := ['CASE1', 'CASE2', 'CASE3']);
```

```sql
-- Region-restricted frequencies
SELECT ID, ALT_FREQ
FROM plink_freq('data/example.pgen',
    region := '22:1-50000000');
```

## See Also

- [plink_hardy](plink_hardy.md) -- HWE test (also uses genotype counts)
- [plink_missing](plink_missing.md) -- missingness analysis
- [Quality Control Guide](../guides/quality-control.md) -- QC workflows using frequency data
