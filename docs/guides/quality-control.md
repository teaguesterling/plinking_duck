# Quality Control Workflows

This guide shows how to build common genomics QC pipelines using PlinkingDuck's SQL functions. These workflows mirror standard PLINK QC steps but expressed as SQL queries.

## Overview

A typical genotype QC pipeline filters on:

1. **Variant missingness** -- remove variants with too many missing calls
2. **Sample missingness** -- remove samples with too many missing calls
3. **Allele frequency** -- remove very rare or monomorphic variants
4. **Hardy-Weinberg equilibrium** -- remove variants that deviate from HWE

## Step 1: Variant Missingness

Remove variants where more than 5% of samples have missing genotypes.

```sql
-- Identify variants with high missingness
SELECT ID, MISSING_CT, OBS_CT, F_MISS
FROM plink_missing('cohort.pgen')
WHERE F_MISS > 0.05;
```

```sql
-- Count how many variants would be removed
SELECT
    COUNT(*) AS total_variants,
    COUNT(*) FILTER (WHERE F_MISS > 0.05) AS removed,
    COUNT(*) FILTER (WHERE F_MISS <= 0.05) AS kept
FROM plink_missing('cohort.pgen');
```

## Step 2: Sample Missingness

Remove samples where more than 10% of variants have missing genotypes.

```sql
-- Identify samples with high missingness
SELECT IID, FID, MISSING_CT, F_MISS
FROM plink_missing('cohort.pgen', mode := 'sample')
WHERE F_MISS > 0.10;
```

```sql
-- Get the list of samples passing QC
SELECT IID
FROM plink_missing('cohort.pgen', mode := 'sample')
WHERE F_MISS <= 0.10;
```

## Step 3: Allele Frequency Filtering

Remove monomorphic variants and very rare variants (MAF < 1%).

```sql
-- Find variants with MAF < 1%
-- MAF = min(ALT_FREQ, 1 - ALT_FREQ)
SELECT ID, ALT_FREQ,
    LEAST(ALT_FREQ, 1.0 - ALT_FREQ) AS MAF
FROM plink_freq('cohort.pgen')
WHERE ALT_FREQ IS NOT NULL
    AND LEAST(ALT_FREQ, 1.0 - ALT_FREQ) < 0.01;
```

```sql
-- Frequency distribution summary
SELECT
    COUNT(*) FILTER (WHERE ALT_FREQ IS NULL) AS monomorphic_or_all_missing,
    COUNT(*) FILTER (WHERE LEAST(ALT_FREQ, 1 - ALT_FREQ) < 0.01) AS rare,
    COUNT(*) FILTER (WHERE LEAST(ALT_FREQ, 1 - ALT_FREQ) >= 0.01) AS common
FROM plink_freq('cohort.pgen');
```

## Step 4: Hardy-Weinberg Equilibrium

Remove variants that deviate significantly from HWE (p < 1e-6).

```sql
-- Variants failing HWE
SELECT ID, P_HWE, O_HET, E_HET, HOM_REF_CT, HET_CT, HOM_ALT_CT
FROM plink_hardy('cohort.pgen')
WHERE P_HWE < 1e-6;
```

!!! tip
    For case-control studies, run HWE filtering on controls only to avoid removing variants associated with disease status.

```sql
-- HWE test on controls only (assuming you have a list of control IIDs)
SELECT ID, P_HWE
FROM plink_hardy('cohort.pgen',
    samples := ['CTRL1', 'CTRL2', 'CTRL3'])
WHERE P_HWE < 1e-6;
```

## Combined QC Report

Generate a single variant-level QC summary by joining frequency, missingness, and HWE results:

```sql
SELECT
    f.CHROM, f.POS, f.ID,
    f.ALT_FREQ,
    LEAST(f.ALT_FREQ, 1.0 - f.ALT_FREQ) AS MAF,
    m.F_MISS,
    h.P_HWE,
    -- QC pass/fail flags
    CASE WHEN m.F_MISS <= 0.05 THEN true ELSE false END AS pass_missingness,
    CASE WHEN LEAST(f.ALT_FREQ, 1 - f.ALT_FREQ) >= 0.01 THEN true ELSE false END AS pass_maf,
    CASE WHEN h.P_HWE >= 1e-6 THEN true ELSE false END AS pass_hwe
FROM plink_freq('cohort.pgen') f
JOIN plink_missing('cohort.pgen') m ON f.ID = m.ID
JOIN plink_hardy('cohort.pgen') h ON f.ID = h.ID;
```

## Applying QC Filters

After identifying variants and samples to keep, use them with the `samples` and `variants` parameters on downstream analyses:

```sql
-- Use QC-passing samples for frequency analysis
WITH good_samples AS (
    SELECT IID
    FROM plink_missing('cohort.pgen', mode := 'sample')
    WHERE F_MISS <= 0.10
)
SELECT f.ID, f.ALT_FREQ
FROM plink_freq('cohort.pgen',
    samples := (SELECT LIST(IID) FROM good_samples));
```

## Suggested Thresholds

These are common starting points; adjust based on your dataset:

| Filter | Typical Threshold | Function |
|--------|-------------------|----------|
| Variant missingness | F_MISS < 0.05 | `plink_missing` (variant mode) |
| Sample missingness | F_MISS < 0.10 | `plink_missing` (sample mode) |
| Minor allele frequency | MAF > 0.01 | `plink_freq` |
| HWE p-value | P_HWE > 1e-6 | `plink_hardy` |
