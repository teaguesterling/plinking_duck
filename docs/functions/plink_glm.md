# plink_glm

Run per-variant genome-wide association (GWAS) regression using plink2's proven regression engine.

## Synopsis

```sql
plink_glm(prefix VARCHAR, phenotype := ...,
          [, covariates := ..., model := ..., firth := ...,
          pvar := ..., psam := ..., samples := ..., region := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `prefix` | `VARCHAR` | *(required)* | Pfile prefix (e.g., `'data/cohort'` resolves to `.pgen`, `.pvar`, `.psam`) |
| `phenotype` | `LIST(DOUBLE)` | *(required)* | Phenotype values, one per sample in pgen order |
| `covariates` | `STRUCT(name LIST(DOUBLE), ...)` | None | Named covariate vectors |
| `model` | `VARCHAR` | `'auto'` | Regression model: `'auto'`, `'linear'`, or `'logistic'` |
| `firth` | `BOOLEAN` | `true` | Enable Firth correction fallback for logistic regression |
| `pvar` | `VARCHAR` | Auto-discovered | Path to `.pvar` or `.bim` file |
| `psam` | `VARCHAR` | Auto-discovered | Path to `.psam` or `.fam` file |
| `samples` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All | Subset to specific samples |
| `region` | `VARCHAR` | All | Filter to genomic region (`chr:start-end`) |

See [Common Parameters](../common-parameters.md) for details on `pvar`, `psam`, `samples`, and `region`.

### Prefix Resolution

The first argument is a **pfile prefix**, matching `read_pfile`'s convention:

- `'data/cohort'` resolves to `data/cohort.pgen`, `data/cohort.pvar`, `data/cohort.psam`
- If the literal path exists as a file (e.g., already has `.pgen` extension), it is used directly

### Phenotype

A list of numeric values, one per sample in pgen file order. Use `NULL` for missing samples (excluded from regression). The list length must match the sample count in the `.pgen` file (or the subset size if `samples` is specified).

### Model Selection

| `model` | Behavior |
|---------|----------|
| `'auto'` | Binary phenotype (only 0/1 or 1/2 values) → logistic; otherwise → linear |
| `'linear'` | Force OLS linear regression regardless of phenotype |
| `'logistic'` | Force logistic regression regardless of phenotype |

Binary 1/2 phenotypes are automatically remapped to 0/1 before logistic regression.

### Covariates

A struct of named covariate vectors. Each value must be a `LIST(DOUBLE)` with the same length as `phenotype`. NULL values in covariates are not permitted.

```sql
covariates := {'age': [25.0, 30.0, 35.0], 'bmi': [22.1, 24.5, 23.0]}
```

### Firth Correction

When `firth` is `true` (the default) and logistic regression fails to converge for a variant, plink2's Firth penalized logistic regression is used as a fallback. This handles quasi-complete separation, which is common with rare variants or small samples.

## Output Columns

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | `VARCHAR` | Chromosome |
| `POS` | `INTEGER` | Base-pair position |
| `ID` | `VARCHAR` | Variant identifier |
| `REF` | `VARCHAR` | Reference allele |
| `ALT` | `VARCHAR` | Alternate allele |
| `A1` | `VARCHAR` | Tested allele (always ALT) |
| `A1_FREQ` | `DOUBLE` | Tested allele frequency among non-missing samples |
| `TEST` | `VARCHAR` | Test type (always `ADD` for additive) |
| `OBS_CT` | `INTEGER` | Number of non-missing samples used |
| `BETA` | `DOUBLE` | Effect size estimate (log-odds for logistic) |
| `SE` | `DOUBLE` | Standard error of BETA |
| `T_STAT` | `DOUBLE` | Test statistic (t for linear, z for logistic) |
| `P` | `DOUBLE` | Association p-value |
| `ERRCODE` | `VARCHAR` | Error code if regression failed, NULL otherwise |
| `OR` | `DOUBLE` | Odds ratio `exp(BETA)` (logistic only, NULL for linear) |
| `FIRTH_YN` | `VARCHAR` | `'Y'` if Firth was used, `'N'` if standard logistic, NULL for linear |

### Error Codes

| Code | Meaning |
|------|---------|
| `TOO_FEW_SAMPLES` | Fewer than 3 non-missing samples (or fewer than predictors + 1) |
| `CONSTANT_GENOTYPE` | All non-missing samples have the same genotype |
| `SINGULAR_MATRIX` | Design matrix is singular (collinear predictors) |
| `NO_CONVERGENCE` | Logistic regression did not converge (and Firth disabled or also failed) |
| `SEPARATION` | Complete separation detected in logistic regression |

When an error occurs, `BETA`, `SE`, `T_STAT`, `P`, `OR` are all NULL.

## Description

`plink_glm` performs per-variant association testing — the core operation in a genome-wide association study (GWAS). For each variant, it fits a regression model with the genotype as the predictor of interest and reports the association statistics.

The function wraps plink2's actual regression engine:

- **Linear regression**: `LinearRegressionInv` with SVD-based matrix inversion
- **Logistic regression**: `LogisticRegressionF` (IRLS with SIMD-accelerated float arithmetic)
- **Firth correction**: `FirthRegressionF` (penalized likelihood, fallback for non-convergence)

This produces results matching `plink2 --glm` output.

### Missing Data

Samples with `NULL` phenotype values or missing genotypes (genotype = 3 in the pgen) are excluded on a per-variant basis. The `OBS_CT` column reports how many samples were actually used.

### Projection Pushdown

If only metadata columns (`CHROM`, `POS`, `ID`, `REF`, `ALT`) are selected, genotype loading and regression are skipped entirely.

## Examples

### Basic Linear GWAS

```sql
SELECT ID, BETA, SE, P
FROM plink_glm('data/cohort',
    phenotype := [1.5, 2.3, 3.7, 0.8])
ORDER BY P;
```

### Logistic GWAS (Binary Phenotype)

```sql
-- Auto-detected as logistic (0/1 phenotype)
SELECT ID, BETA, "OR", P, FIRTH_YN
FROM plink_glm('data/cohort',
    phenotype := [0, 1, 0, 1, 1, 0, 1, 0])
WHERE P < 0.05
ORDER BY P;
```

### With Covariates

```sql
SELECT ID, BETA, SE, P
FROM plink_glm('data/cohort',
    phenotype := [1.2, 3.4, 2.1, 5.6, 4.3, 0.9, 3.8, 2.7],
    covariates := {
        'age': [25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0],
        'bmi': [22.1, 24.5, 23.0, 28.3, 26.1, 21.5, 25.8, 23.2]
    })
WHERE P < 5e-8;
```

### Phenotype and Covariates from Parquet Files

For real datasets, phenotypes and covariates live in external files. Use DuckDB's `SET VARIABLE` to load them, then pass with `getvariable()`.

**Important**: The phenotype and covariate lists must be in the same order as samples in the `.pgen` file.

#### Approach 1: Pre-Sorted Data (Positional)

If your parquet file is already sorted to match pgen sample order (e.g., it was exported from the same pipeline), you can load directly without a join:

```sql
-- Phenotype file rows are in pgen sample order
SET VARIABLE pheno = (
    SELECT list(pheno_quant) FROM 'phenotypes.parquet'
);

SET VARIABLE covars = (
    SELECT {
        'age': list(age),
        'bmi': list(bmi),
        'pc1': list(pc1)
    } FROM 'covariates.parquet'
);

SELECT ID, BETA, P
FROM plink_glm('data/cohort',
    phenotype := getvariable('pheno'),
    covariates := getvariable('covars'))
WHERE P < 5e-8
ORDER BY P;
```

#### Approach 2: Join on IID via psam (Safe)

If your phenotype file uses sample IDs (IID) that may not match pgen order, join against the `.psam` file to ensure correct alignment:

```sql
-- Load psam to get pgen sample order
SET VARIABLE pheno = (
    SELECT list(p.pheno_quant ORDER BY s.rowid)
    FROM read_psam('data/cohort.psam') s
    JOIN 'phenotypes.parquet' p ON s.IID = p.IID
);

SET VARIABLE covars = (
    SELECT {
        'age': list(p.age ORDER BY s.rowid),
        'bmi': list(p.bmi ORDER BY s.rowid),
        'pc1': list(p.pc1 ORDER BY s.rowid)
    }
    FROM read_psam('data/cohort.psam') s
    JOIN 'covariates.parquet' p ON s.IID = p.IID
);

SELECT ID, CHROM, POS, BETA, SE, P
FROM plink_glm('data/cohort',
    phenotype := getvariable('pheno'),
    covariates := getvariable('covars'))
WHERE P < 5e-8
ORDER BY P;
```

The `ORDER BY s.rowid` ensures values are emitted in the physical `.psam` row order, which matches the `.pgen` sample order.

### Region-Restricted Analysis

```sql
-- Only test variants on chromosome 22
SELECT ID, BETA, P
FROM plink_glm('data/cohort',
    phenotype := getvariable('pheno'),
    region := '22:1-50000000')
ORDER BY P;
```

### Composing with Other Functions

```sql
-- Join GWAS results with allele frequencies
SELECT g.ID, g.BETA, g.P, f.ALT_FREQ
FROM plink_glm('data/cohort',
    phenotype := [1.5, 2.3, 3.7, 0.8]) g
JOIN plink_freq('data/cohort.pgen') f ON g.ID = f.ID
WHERE g.P < 0.05
ORDER BY g.P;
```

### Manhattan Plot Data

```sql
-- Prepare data for a Manhattan plot
SELECT CHROM, POS, -LOG10(P) AS neg_log_p
FROM plink_glm('data/cohort',
    phenotype := getvariable('pheno'),
    covariates := getvariable('covars'))
WHERE ERRCODE IS NULL
ORDER BY CHROM, POS;
```

## See Also

- [plink_freq](plink_freq.md) -- allele frequencies (useful for QC before GWAS)
- [plink_hardy](plink_hardy.md) -- HWE test (QC filter)
- [plink_score](plink_score.md) -- polygenic risk scoring (downstream of GWAS)
- [Quality Control Guide](../guides/quality-control.md) -- pre-GWAS QC workflows
