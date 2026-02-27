# plink_hardy

Compute per-variant Hardy-Weinberg equilibrium exact test p-values.

## Synopsis

```sql
plink_hardy(path VARCHAR [, pvar := ..., psam := ..., samples := ...,
            region := ..., midp := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `path` | `VARCHAR` | *(required)* | Path to the `.pgen` file |
| `pvar` | `VARCHAR` | Auto-discovered | Path to `.pvar` or `.bim` file |
| `psam` | `VARCHAR` | Auto-discovered | Path to `.psam` or `.fam` file |
| `samples` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All | Subset to specific samples |
| `region` | `VARCHAR` | All | Filter to genomic region (`chr:start-end`) |
| `midp` | `BOOLEAN` | `false` | Use mid-p correction for the exact test |

See [Common Parameters](../common-parameters.md) for details on `pvar`, `psam`, `samples`, and `region`.

## Output Columns

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | `VARCHAR` | Chromosome |
| `POS` | `INTEGER` | Base-pair position |
| `ID` | `VARCHAR` | Variant identifier |
| `REF` | `VARCHAR` | Reference allele |
| `ALT` | `VARCHAR` | Alternate allele |
| `A1` | `VARCHAR` | Tested allele (the alternate allele) |
| `HOM_REF_CT` | `INTEGER` | Homozygous reference count |
| `HET_CT` | `INTEGER` | Heterozygous count |
| `HOM_ALT_CT` | `INTEGER` | Homozygous alternate count |
| `O_HET` | `DOUBLE` | Observed heterozygosity (NULL if no observations) |
| `E_HET` | `DOUBLE` | Expected heterozygosity under HWE (NULL if no observations) |
| `P_HWE` | `DOUBLE` | HWE exact test p-value (NULL if no observations) |

## Description

`plink_hardy` performs a Hardy-Weinberg equilibrium exact test for each variant, implementing the algorithm described by Wigginton et al. (2005). The test enumerates all valid heterozygote counts for fixed allele counts and computes p-values by summing probabilities of configurations as likely or less likely than observed.

### Heterozygosity

- **O_HET** (observed): `HET_CT / (HOM_REF_CT + HET_CT + HOM_ALT_CT)`
- **E_HET** (expected): `2 * p * q` where `p = (2*HOM_REF_CT + HET_CT) / (2*total)` and `q = 1 - p`

### Mid-p Correction

When `midp := true`, the p-value is adjusted by subtracting half the probability of the observed configuration. This reduces conservativeness of the exact test and is recommended for large sample sizes.

### Projection Pushdown

If only metadata columns (CHROM through ALT) are selected, genotype counting and the HWE test are skipped.

## Examples

```sql
-- Basic HWE test
SELECT * FROM plink_hardy('data/example.pgen');
```

```sql
-- With mid-p correction
SELECT ID, P_HWE
FROM plink_hardy('data/example.pgen', midp := true);
```

```sql
-- QC filter: find variants failing HWE
SELECT ID, P_HWE, O_HET, E_HET
FROM plink_hardy('data/example.pgen')
WHERE P_HWE < 1e-6;
```

```sql
-- HWE test restricted to a region
SELECT ID, P_HWE
FROM plink_hardy('data/example.pgen',
    region := '1:1-50000000');
```

```sql
-- HWE test in a sample subset (e.g., controls only)
SELECT ID, P_HWE
FROM plink_hardy('data/example.pgen',
    samples := ['CTRL1', 'CTRL2', 'CTRL3']);
```

## See Also

- [plink_freq](plink_freq.md) -- allele frequencies
- [plink_missing](plink_missing.md) -- missingness analysis
- [Quality Control Guide](../guides/quality-control.md) -- HWE filtering in QC pipelines
