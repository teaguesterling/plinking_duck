# plink_score

Compute polygenic risk scores (PRS) from genotype data and variant weights.

## Synopsis

```sql
plink_score(path VARCHAR, weights := ...,
            [, pvar := ..., psam := ..., samples := ..., region := ...,
            center := ..., no_mean_imputation := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `path` | `VARCHAR` | *(required)* | Path to the `.pgen` file |
| `weights` | `LIST(DOUBLE)` or `LIST(STRUCT)` | *(required)* | Variant scoring weights |
| `pvar` | `VARCHAR` | Auto-discovered | Path to `.pvar` or `.bim` file |
| `psam` | `VARCHAR` | Auto-discovered | Path to `.psam` or `.fam` file (required) |
| `samples` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All | Subset to specific samples |
| `region` | `VARCHAR` | All | Filter to genomic region (`chr:start-end`) |
| `center` | `BOOLEAN` | `false` | Variance-standardized scoring |
| `no_mean_imputation` | `BOOLEAN` | `false` | Skip mean imputation for missing genotypes |

See [Common Parameters](../common-parameters.md) for details on `pvar`, `psam`, `samples`, and `region`.

### Weights Parameter

The `weights` parameter accepts two formats:

#### Positional Mode: `LIST(DOUBLE)`

A list of numeric weights, one per variant in order. The list length must match the variant count (or the region-filtered variant count). Zero-weight variants are skipped.

```sql
plink_score('data.pgen', weights := [0.5, -0.3, 0.0, 1.2])
```

#### ID-keyed Mode: `LIST(STRUCT(id VARCHAR, allele VARCHAR, weight DOUBLE))`

A list of structs specifying variant ID, scored allele, and weight. Variants not found in the `.pvar` file are silently skipped. The `allele` field determines orientation:

- If `allele` matches ALT: dosage is used as-is
- If `allele` matches REF: dosage is flipped (`2 - alt_dosage`)
- If `allele` matches neither: the variant is skipped

```sql
plink_score('data.pgen', weights := [
    {'id': 'rs1', 'allele': 'A', 'weight': 0.5},
    {'id': 'rs2', 'allele': 'T', 'weight': -0.3}
])
```

## Output Columns

| Column | Type | Description |
|--------|------|-------------|
| `FID` | `VARCHAR` | Family ID (NULL if no FID column in `.psam`) |
| `IID` | `VARCHAR` | Individual ID |
| `ALLELE_CT` | `INTEGER` | Total alleles observed (2 per non-missing scored variant) |
| `DENOM` | `INTEGER` | Denominator for score averaging (same as ALLELE_CT) |
| `NAMED_ALLELE_DOSAGE_SUM` | `DOUBLE` | Sum of scored allele dosages across variants |
| `SCORE_SUM` | `DOUBLE` | Sum of weight * dosage across all scored variants |
| `SCORE_AVG` | `DOUBLE` | `SCORE_SUM / ALLELE_CT` (0 if ALLELE_CT is 0) |

## Description

`plink_score` computes a polygenic risk score for each sample by reading genotype dosages (via pgenlib's `PgrGetD`) and applying variant-specific weights.

### Scoring Formula

For each sample, the score is:

```
SCORE_SUM = sum(weight_i * dosage_i) for all scored variants i
SCORE_AVG = SCORE_SUM / ALLELE_CT
```

Where `dosage_i` is the allele dosage (0, 1, or 2 for hardcalled genotypes).

### Missing Data Handling

Three modes control how missing genotypes are handled:

| Mode | Behavior |
|------|----------|
| Default (mean imputation) | Missing dosages are replaced with the variant's mean dosage |
| `no_mean_imputation := true` | Missing samples are skipped (their ALLELE_CT reflects only non-missing variants) |
| `center := true` | Variance-standardized scoring: `(dosage - mean) / sd`. Missing samples are skipped. Monomorphic variants are excluded. |

`center` and `no_mean_imputation` cannot both be `true`.

### Requirements

- The `.psam` file is **required** (sample IDs are needed for output)
- The `weights` parameter is **required**

## Examples

```sql
-- Positional weights (one per variant)
SELECT IID, SCORE_SUM, SCORE_AVG
FROM plink_score('data/example.pgen',
    weights := [0.5, -0.3, 1.2, 0.0]);
```

```sql
-- ID-keyed weights with allele specification
SELECT IID, SCORE_SUM
FROM plink_score('data/example.pgen', weights := [
    {'id': 'rs1', 'allele': 'A', 'weight': 0.5},
    {'id': 'rs2', 'allele': 'T', 'weight': -0.3}
]);
```

```sql
-- Variance-standardized scoring
SELECT IID, SCORE_AVG
FROM plink_score('data/example.pgen',
    weights := [1.0, 1.0, 1.0, 1.0],
    center := true);
```

```sql
-- Without mean imputation
SELECT IID, SCORE_SUM, ALLELE_CT
FROM plink_score('data/example.pgen',
    weights := [0.5, -0.3, 1.2, 0.8],
    no_mean_imputation := true);
```

```sql
-- Score a subset of samples
SELECT IID, SCORE_AVG
FROM plink_score('data/example.pgen',
    weights := [0.5, 0.3, 0.2, 0.1],
    samples := ['SAMPLE1', 'SAMPLE2']);
```

## See Also

- [plink_freq](plink_freq.md) -- allele frequencies (useful for QC before scoring)
- [Genotype Encoding](../genotype-encoding.md) -- encoding reference
