# plink_glm

Run per-variant association testing (linear and logistic regression) by wrapping plink2's GLM implementation.

!!! note "Coming Soon"

    `plink_glm` is under active development. This page documents the planned interface and behavior.

## Synopsis

```sql
plink_glm(path VARCHAR, phenotype := ...,
          [, covariates := ..., pvar := ..., psam := ...,
          samples := ..., region := ..., model := ...,
          firth := ..., max_corr := ..., vif_thresh := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `path` | `VARCHAR` | *(required)* | Path to the `.pgen` file |
| `phenotype` | `LIST(DOUBLE)` | *(required)* | Quantitative or binary phenotype, positionally aligned with `.psam` |
| `covariates` | `STRUCT(LIST(DOUBLE), ...)` | None | Named covariate vectors (struct-of-lists) |
| `pvar` | `VARCHAR` | Auto-discovered | Path to `.pvar` or `.bim` file |
| `psam` | `VARCHAR` | Auto-discovered | Path to `.psam` or `.fam` file |
| `samples` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All | Subset to specific samples |
| `region` | `VARCHAR` | All | Filter to genomic region (`chr:start-end`) |
| `model` | `VARCHAR` | `'auto'` | Force `'linear'` or `'logistic'` (auto-detect from phenotype) |
| `firth` | `BOOLEAN` | `true` | Use Firth correction for logistic regression |
| `max_corr` | `DOUBLE` | `0.999` | Max covariate correlation before exclusion |
| `vif_thresh` | `DOUBLE` | `50.0` | Max variance inflation factor |

See [Common Parameters](../common-parameters.md) for details on `pvar`, `psam`, `samples`, and `region`.

### Phenotype Type Detection

The regression model is auto-detected from the phenotype values:

- If all non-NULL values are 0 or 1: **binary** phenotype, logistic regression
- If more than 2 distinct values: **quantitative** phenotype, linear regression
- Use `model := 'linear'` or `model := 'logistic'` to override auto-detection

### Covariates Parameter

The `covariates` parameter accepts a struct-of-lists, where each key is the covariate name and each value is a `LIST(DOUBLE)` positionally aligned with the `.psam` file:

```sql
covariates := {'sex': list(sex), 'pc1': list(pc1), 'pc2': list(pc2)}
```

## Output Columns

### Linear Regression

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | `VARCHAR` | Chromosome |
| `POS` | `INTEGER` | Base-pair position |
| `ID` | `VARCHAR` | Variant ID |
| `REF` | `VARCHAR` | Reference allele |
| `ALT` | `VARCHAR` | Alternate allele |
| `A1` | `VARCHAR` | Tested allele (ALT) |
| `A1_FREQ` | `DOUBLE` | Tested allele frequency |
| `TEST` | `VARCHAR` | Test name (`'ADD'`) |
| `OBS_CT` | `INTEGER` | Non-missing sample count for this variant |
| `BETA` | `DOUBLE` | Effect size (regression coefficient) |
| `SE` | `DOUBLE` | Standard error of beta |
| `T_STAT` | `DOUBLE` | T-statistic (`BETA / SE`) |
| `P` | `DOUBLE` | P-value |
| `ERRCODE` | `VARCHAR` | Error code (NULL if no error) |

### Logistic Regression

Same columns as linear regression, plus:

| Column | Type | Description |
|--------|------|-------------|
| `OR` | `DOUBLE` | Odds ratio (`exp(BETA)`) |
| `FIRTH_YN` | `VARCHAR` | `'Y'` if Firth correction was applied |

## Description

`plink_glm` performs per-variant association testing using plink2's battle-tested regression implementation. It wraps the actual plink2 GLM compute code rather than reimplementing the statistics, ensuring correctness and taking advantage of plink2's SIMD-optimized genotype handling.

The function produces one row per variant with the additive test result. Output column names match the standard plink2 `--glm` output format.

### Why Wrap plink2

plink2's `--glm` implementation is 18K lines of thoroughly tested, peer-reviewed statistical code with:

- SIMD-optimized genotype unpacking (2-bit to dosage, no DOUBLE conversion)
- Precomputed regression matrices for the no-missing-data fast path
- Firth correction for logistic regression on rare variants
- Proper handling of dosage data, covariates, and sex chromosomes

### Missing Data

Missing genotypes reduce the per-variant `OBS_CT`. The regression is computed only over non-missing samples for each variant.

## Examples

```sql
-- Simple linear regression: quantitative phenotype, no covariates
SELECT ID, BETA, SE, P
FROM plink_glm('data/cohort.pgen',
    phenotype := list(age));
```

```sql
-- Logistic regression: binary phenotype (auto-detected)
SELECT ID, OR, P
FROM plink_glm('data/cohort.pgen',
    phenotype := list(case_status));
```

```sql
-- With covariates
SELECT ID, BETA, SE, P
FROM plink_glm('data/cohort.pgen',
    phenotype := list(age),
    covariates := {'sex': list(sex), 'pc1': list(pc1), 'pc2': list(pc2)});
```

```sql
-- Phenotype from a parquet file
SET VARIABLE pheno = (SELECT list(age) FROM 'phenotypes/age.parquet');
SELECT ID, CHROM, POS, BETA, SE, P
FROM plink_glm('data/cohort.pgen',
    phenotype := getvariable('pheno'))
WHERE P < 5e-8;
```

```sql
-- Composable with other PlinkingDuck functions
SELECT g.ID, g.P, f.ALT_FREQ
FROM plink_glm('data/cohort.pgen',
    phenotype := getvariable('pheno')) g
JOIN plink_freq('data/cohort.pgen') f ON g.ID = f.ID
WHERE g.P < 5e-8 AND f.ALT_FREQ > 0.01;
```

```sql
-- Region-filtered GWAS
SELECT ID, BETA, P
FROM plink_glm('data/cohort.pgen',
    phenotype := list(trait),
    region := '22:1-50000000');
```

## See Also

- [plink_freq](plink_freq.md) -- allele frequencies (useful for QC before association testing)
- [plink_hardy](plink_hardy.md) -- HWE test (QC filter for association studies)
- [plink_score](plink_score.md) -- polygenic risk scoring from GWAS weights
