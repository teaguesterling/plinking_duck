# Tutorial: Genomic Analysis with SQL

This tutorial walks through a complete genomic analysis workflow using PlinkingDuck — from loading data through quality control to downstream analysis. By the end, you'll know how to explore genomic datasets, assess data quality, measure variant correlation, and compute polygenic risk scores, all in SQL.

No prior genomics knowledge is required. We'll explain concepts as they arise.

!!! info "Prerequisites"
    Before starting, make sure PlinkingDuck is installed and loaded. See the [Getting Started](getting-started.md) guide for setup instructions.

---

## 1. Meet the Dataset

Genomic data in PLINK 2 format is stored as a **fileset** — a group of files sharing a common prefix:

- **`.pvar`** — variant metadata (which genetic markers were measured)
- **`.psam`** — sample metadata (who was measured)
- **`.pgen`** — genotype data (the actual measurements, in compact binary)

Our tutorial dataset has 4 genetic variants measured across 4 individuals. Let's start by examining what variants and samples are present.

### Variant Metadata

Each variant (a position in the genome where people differ) has a chromosome, position, identifier, and the two alleles observed:

```sql
SELECT * FROM read_pvar('test/data/pgen_example.pvar');
```

| CHROM | POS   | ID  | REF | ALT |
|-------|-------|-----|-----|-----|
| 1     | 10000 | rs1 | A   | G   |
| 1     | 20000 | rs2 | C   | T   |
| 1     | 30000 | rs3 | G   | A   |
| 2     | 15000 | rs4 | T   | C   |

Three variants sit on chromosome 1 and one on chromosome 2. **REF** is the reference allele (the common form) and **ALT** is the alternate allele.

!!! tip "What are rsIDs?"
    Variant identifiers like `rs1` and `rs2` are shorthand here. In real data, these are [RefSNP IDs](https://www.ncbi.nlm.nih.gov/snp/) (e.g., `rs7903146`) that uniquely identify known variants across studies.

Since this is standard SQL, you can summarize the data however you like:

```sql
SELECT CHROM, COUNT(*) AS n_variants
FROM read_pvar('test/data/pgen_example.pvar')
GROUP BY CHROM
ORDER BY CHROM;
```

| CHROM | n_variants |
|-------|------------|
| 1     | 3          |
| 2     | 1          |

### Sample Metadata

Now let's look at who was genotyped. We'll use a companion `.psam` file that includes family IDs and sex:

```sql
SELECT * FROM read_psam('test/data/pfile_example.psam');
```

| FID    | IID     | SEX  |
|--------|---------|------|
| FAM001 | SAMPLE1 | 1    |
| FAM001 | SAMPLE2 | 2    |
| FAM002 | SAMPLE3 | NULL |
| FAM002 | SAMPLE4 | 1    |

**IID** is the individual identifier. **FID** groups individuals into families. **SEX** uses PLINK conventions: 1 = male, 2 = female, NULL = unknown.

!!! note
    Different `.psam` files can have different columns. The minimal format has just `IID`. See [read_psam](functions/read_psam.md) for format details.

---

## 2. Reading Genotypes

Now for the core data — the actual genotype calls. A **genotype** is an individual's genetic state at a variant: how many copies of the alternate allele they carry.

### List Mode (One Row per Variant)

`read_pgen` returns one row per variant, with all genotypes packed into a list:

```sql
SELECT * FROM read_pgen('test/data/pgen_example.pgen');
```

| CHROM | POS   | ID  | REF | ALT | genotypes       |
|-------|-------|-----|-----|-----|-----------------|
| 1     | 10000 | rs1 | A   | G   | [0, 1, 2, NULL] |
| 1     | 20000 | rs2 | C   | T   | [1, 1, 0, 2]    |
| 1     | 30000 | rs3 | G   | A   | [2, NULL, 1, 0] |
| 2     | 15000 | rs4 | T   | C   | [0, 0, 1, 2]    |

Each number in the genotype list is an **alternate allele count**: 0 = homozygous reference (no copies of ALT), 1 = heterozygous (one copy), 2 = homozygous alternate (two copies), NULL = missing. See [Genotype Encoding](genotype-encoding.md) for the full reference.

For variant rs1, SAMPLE1 has genotype 0 (AA), SAMPLE2 has 1 (AG), SAMPLE3 has 2 (GG), and SAMPLE4 is missing.

### Tidy Mode (One Row per Observation)

For analysis with SQL joins and aggregations, **tidy mode** is often more convenient. It produces one row per variant-sample pair:

```sql
SELECT CHROM, POS, ID, IID, genotype
FROM read_pfile('test/data/pfile_example', tidy := true)
ORDER BY ID, IID;
```

| CHROM | POS   | ID  | IID     | genotype |
|-------|-------|-----|---------|----------|
| 1     | 10000 | rs1 | SAMPLE1 | 0        |
| 1     | 10000 | rs1 | SAMPLE2 | 1        |
| 1     | 10000 | rs1 | SAMPLE3 | 2        |
| 1     | 10000 | rs1 | SAMPLE4 | NULL     |
| 1     | 20000 | rs2 | SAMPLE1 | 1        |
| 1     | 20000 | rs2 | SAMPLE2 | 1        |
| 1     | 20000 | rs2 | SAMPLE3 | 0        |
| 1     | 20000 | rs2 | SAMPLE4 | 2        |
| 1     | 30000 | rs3 | SAMPLE1 | 2        |
| 1     | 30000 | rs3 | SAMPLE2 | NULL     |
| 1     | 30000 | rs3 | SAMPLE3 | 1        |
| 1     | 30000 | rs3 | SAMPLE4 | 0        |
| 2     | 15000 | rs4 | SAMPLE1 | 0        |
| 2     | 15000 | rs4 | SAMPLE2 | 0        |
| 2     | 15000 | rs4 | SAMPLE3 | 1        |
| 2     | 15000 | rs4 | SAMPLE4 | 2        |

!!! note "`read_pfile` vs `read_pgen`"
    `read_pfile` takes a **prefix** (e.g., `'test/data/pfile_example'`) and auto-discovers the `.pgen`, `.pvar`, and `.psam` files. `read_pgen` takes the `.pgen` path directly. See [read_pfile](functions/read_pfile.md) and [read_pgen](functions/read_pgen.md) for details.

With tidy data, standard SQL works naturally. For example, count the total alternate alleles carried by each individual:

```sql
SELECT IID, SUM(genotype) AS total_alt_alleles
FROM read_pfile('test/data/pfile_example', tidy := true)
WHERE genotype IS NOT NULL
GROUP BY IID
ORDER BY IID;
```

| IID     | total_alt_alleles |
|---------|-------------------|
| SAMPLE1 | 3                 |
| SAMPLE2 | 2                 |
| SAMPLE3 | 4                 |
| SAMPLE4 | 4                 |

Or find all heterozygous calls (genotype = 1):

```sql
SELECT ID, IID
FROM read_pfile('test/data/pfile_example', tidy := true)
WHERE genotype = 1
ORDER BY ID, IID;
```

| ID  | IID     |
|-----|---------|
| rs1 | SAMPLE2 |
| rs2 | SAMPLE1 |
| rs2 | SAMPLE2 |
| rs3 | SAMPLE3 |
| rs4 | SAMPLE3 |

---

## 3. Quality Control

Before running any analysis, you need to check data quality. PlinkingDuck provides dedicated functions for the three pillars of genotype QC: **missingness**, **allele frequency**, and **Hardy-Weinberg equilibrium**. These use optimized paths in pgenlib that are faster than computing from raw genotypes.

### Missingness

**Missingness** is the fraction of genotype calls that failed. High missingness for a variant suggests a technical problem with that assay; high missingness for a sample suggests a problem with that individual's DNA.

#### Per-Variant Missingness

```sql
SELECT ID, MISSING_CT, OBS_CT, F_MISS
FROM plink_missing('test/data/pgen_example.pgen');
```

| ID  | MISSING_CT | OBS_CT | F_MISS |
|-----|------------|--------|--------|
| rs1 | 1          | 3      | 0.25   |
| rs2 | 0          | 4      | 0.0    |
| rs3 | 1          | 3      | 0.25   |
| rs4 | 0          | 4      | 0.0    |

Variants rs1 and rs3 each have one missing call (25% missingness). In a real study, you would typically remove variants exceeding a threshold (e.g., 5%):

```sql
SELECT ID, F_MISS
FROM plink_missing('test/data/pgen_example.pgen')
WHERE F_MISS > 0.05;
```

#### Per-Sample Missingness

Switch to sample mode to check individual-level data quality:

```sql
SELECT IID, MISSING_CT, OBS_CT, F_MISS
FROM plink_missing('test/data/pgen_example.pgen',
    mode := 'sample', psam := 'test/data/pfile_example.psam');
```

| IID     | MISSING_CT | OBS_CT | F_MISS |
|---------|------------|--------|--------|
| SAMPLE1 | 0          | 4      | 0.0    |
| SAMPLE2 | 1          | 3      | 0.25   |
| SAMPLE3 | 0          | 4      | 0.0    |
| SAMPLE4 | 1          | 3      | 0.25   |

SAMPLE2 and SAMPLE4 each have one missing genotype.

### Allele Frequency

**Allele frequency** is how common each allele is in the population. The alternate allele frequency (ALT_FREQ) ranges from 0 to 1. Variants with very low frequency (rare variants) may need special handling.

```sql
SELECT ID, ALT_FREQ, OBS_CT
FROM plink_freq('test/data/pgen_example.pgen');
```

| ID  | ALT_FREQ | OBS_CT |
|-----|----------|--------|
| rs1 | 0.5      | 6      |
| rs2 | 0.5      | 8      |
| rs3 | 0.5      | 6      |
| rs4 | 0.375    | 8      |

!!! tip "OBS_CT is allele count, not sample count"
    `OBS_CT` is the number of observed **alleles** (2 per non-missing sample, since humans are diploid). With 3 non-missing samples, OBS_CT = 6.

Use `counts := true` to see the underlying genotype counts:

```sql
SELECT ID, HOM_REF_CT, HET_CT, HOM_ALT_CT, MISSING_CT
FROM plink_freq('test/data/pgen_example.pgen', counts := true);
```

| ID  | HOM_REF_CT | HET_CT | HOM_ALT_CT | MISSING_CT |
|-----|------------|--------|------------|------------|
| rs1 | 1          | 1      | 1          | 1          |
| rs2 | 1          | 2      | 1          | 0          |
| rs3 | 1          | 1      | 1          | 1          |
| rs4 | 2          | 1      | 1          | 0          |

### Hardy-Weinberg Equilibrium

**Hardy-Weinberg equilibrium (HWE)** predicts the expected genotype frequencies given the allele frequencies. Departures from HWE can indicate genotyping errors, population structure, or selection. The HWE exact test produces a p-value — small values suggest the variant deviates from equilibrium.

```sql
SELECT ID, HOM_REF_CT, HET_CT, HOM_ALT_CT,
       ROUND(P_HWE, 4) AS P_HWE
FROM plink_hardy('test/data/pgen_example.pgen');
```

| ID  | HOM_REF_CT | HET_CT | HOM_ALT_CT | P_HWE  |
|-----|------------|--------|------------|--------|
| rs1 | 1          | 1      | 1          | 1.0    |
| rs2 | 1          | 2      | 1          | 1.0    |
| rs3 | 1          | 1      | 1          | 1.0    |
| rs4 | 2          | 1      | 1          | 0.4286 |

All variants are in equilibrium (p-values well above any reasonable threshold). In practice, variants with P_HWE < 1e-6 are often excluded.

### Putting It Together: A QC Summary

The power of SQL shines when you combine multiple analyses. This CTE joins frequency, missingness, and HWE results into a single QC summary:

```sql
WITH freq AS (
    SELECT ID, ALT_FREQ
    FROM plink_freq('test/data/pgen_example.pgen')
),
miss AS (
    SELECT ID, F_MISS
    FROM plink_missing('test/data/pgen_example.pgen')
),
hardy AS (
    SELECT ID, P_HWE
    FROM plink_hardy('test/data/pgen_example.pgen')
)
SELECT f.ID, f.ALT_FREQ, m.F_MISS, h.P_HWE
FROM freq f
JOIN miss m ON f.ID = m.ID
JOIN hardy h ON f.ID = h.ID
ORDER BY f.ID;
```

| ID  | ALT_FREQ | F_MISS | P_HWE             |
|-----|----------|--------|---------------------|
| rs1 | 0.5      | 0.25   | 1.0                 |
| rs2 | 0.5      | 0.0    | 1.0                 |
| rs3 | 0.5      | 0.25   | 1.0                 |
| rs4 | 0.375    | 0.0    | 0.4285714285714286  |

From here you could filter to only variants passing QC:

```sql
-- Variants passing all QC thresholds
-- (not run here since our tiny dataset passes everything)
WHERE f.ALT_FREQ > 0.01          -- MAF > 1%
  AND m.F_MISS < 0.05            -- missingness < 5%
  AND h.P_HWE > 1e-6             -- HWE p-value > 1e-6
```

!!! note
    For a deeper dive into QC workflows, see the [Quality Control Guide](guides/quality-control.md).

---

## 4. Linkage Disequilibrium

**Linkage disequilibrium (LD)** measures the statistical correlation between genotypes at different variants. Variants in high LD tend to be inherited together. LD is fundamental to understanding genome structure and interpreting association results.

### Pairwise LD

To measure LD between two specific variants, use pairwise mode:

```sql
SELECT ID_A, ID_B, ROUND(R2, 4) AS R2, D_PRIME, OBS_CT
FROM plink_ld('test/data/pgen_example.pgen',
    variant1 := 'rs1', variant2 := 'rs2');
```

| ID_A | ID_B | R2     | D_PRIME | OBS_CT |
|------|------|--------|---------|--------|
| rs1  | rs2  | 0.75   | 0.5     | 3      |

**r²** is the squared correlation between genotypes (0 = independent, 1 = perfectly correlated). **D'** is a normalized measure of allelic association. `OBS_CT` is the number of samples with non-missing genotypes at both variants.

### Windowed LD

To scan all variant pairs within a genomic window, omit the variant parameters. Set `r2_threshold` to control which pairs are reported:

```sql
SELECT ID_A, ID_B, ROUND(R2, 4) AS R2,
       ROUND(D_PRIME, 4) AS D_PRIME, OBS_CT
FROM plink_ld('test/data/pgen_example.pgen',
    r2_threshold := 0.0);
```

| ID_A | ID_B | R2     | D_PRIME | OBS_CT |
|------|------|--------|---------|--------|
| rs1  | rs2  | 0.75   | 0.5     | 3      |
| rs1  | rs3  | 1.0    | 1.0     | 2      |
| rs2  | rs3  | 0.25   | 0.3333  | 3      |

By default, only same-chromosome pairs within 1000 kb are tested. Variants rs1 and rs3 show perfect LD (r² = 1.0), though with only 2 shared observations.

### Cross-Chromosome LD

To include pairs across chromosomes (useful for checking for population structure artifacts), set `inter_chr := true`:

```sql
SELECT ID_A, ID_B, ROUND(R2, 4) AS R2, OBS_CT
FROM plink_ld('test/data/pgen_example.pgen',
    r2_threshold := 0.0, inter_chr := true);
```

| ID_A | ID_B | R2     | OBS_CT |
|------|------|--------|--------|
| rs1  | rs2  | 0.75   | 3      |
| rs1  | rs3  | 1.0    | 2      |
| rs1  | rs4  | 0.75   | 3      |
| rs2  | rs3  | 0.25   | 3      |
| rs2  | rs4  | 0.1818 | 4      |
| rs3  | rs4  | 1.0    | 3      |

!!! warning
    High LD between variants on different chromosomes is unexpected in real data and usually signals population stratification or genotyping artifacts.

---

## 5. Polygenic Risk Scoring

A **polygenic risk score (PRS)** summarizes an individual's genetic predisposition by weighting each variant's genotype by its effect size from a prior study. PlinkingDuck's `plink_score` computes this efficiently, handling missing data with mean imputation by default.

### Positional Weights

The simplest approach provides one weight per variant, in order:

```sql
SELECT IID, SCORE_SUM, SCORE_AVG
FROM plink_score('test/data/pgen_example.pgen',
    psam := 'test/data/pfile_example.psam',
    weights := [0.5, -0.3, 1.2, 0.8]);
```

| IID     | SCORE_SUM | SCORE_AVG |
|---------|-----------|-----------|
| SAMPLE1 | 2.1       | 0.2625    |
| SAMPLE2 | 1.4       | 0.175     |
| SAMPLE3 | 3.0       | 0.375     |
| SAMPLE4 | 1.5       | 0.1875    |

Each score is `sum(weight × dosage)` across all variants. `SCORE_AVG` normalizes by the number of observed alleles.

### ID-Keyed Weights

For real applications, you typically have a scoring file that maps variant IDs and alleles to weights. ID-keyed mode handles this — variants not found in the data are silently skipped:

```sql
SELECT IID, ROUND(SCORE_SUM, 2) AS SCORE_SUM,
       ROUND(SCORE_AVG, 4) AS SCORE_AVG
FROM plink_score('test/data/pgen_example.pgen',
    psam := 'test/data/pfile_example.psam',
    weights := [
        {'id': 'rs1', 'allele': 'G', 'weight': 0.5},
        {'id': 'rs2', 'allele': 'T', 'weight': -0.3},
        {'id': 'rs4', 'allele': 'C', 'weight': 0.8}
    ]);
```

| IID     | SCORE_SUM | SCORE_AVG |
|---------|-----------|-----------|
| SAMPLE1 | -0.3      | -0.05     |
| SAMPLE2 | 0.2       | 0.0333    |
| SAMPLE3 | 1.8       | 0.3       |
| SAMPLE4 | 1.5       | 0.25      |

The `allele` field specifies which allele the weight applies to. If it matches ALT, the dosage is used directly; if it matches REF, the dosage is flipped.

!!! tip "Missing data handling"
    By default, missing genotypes are imputed with the variant's mean dosage. Use `no_mean_imputation := true` to skip missing variants instead (each sample's `ALLELE_CT` will reflect only non-missing variants). See [plink_score](functions/plink_score.md) for details.

---

## 6. Working at Scale

The examples above used a tiny 4-variant dataset. Real genomic data has millions of variants and thousands of samples. PlinkingDuck is designed for this — parallel scanning, projection pushdown, and region filtering all help keep queries fast.

Let's demonstrate with a larger dataset: 3,000 variants across 8 samples on 3 chromosomes.

```sql
SELECT
    (SELECT COUNT(*) FROM read_pvar('test/data/large_example.pvar')) AS n_variants,
    (SELECT COUNT(*) FROM read_psam('test/data/large_example.psam')) AS n_samples;
```

| n_variants | n_samples |
|------------|-----------|
| 3000       | 8         |

```sql
SELECT CHROM, COUNT(*) AS n_variants
FROM read_pvar('test/data/large_example.pvar')
GROUP BY CHROM
ORDER BY CHROM;
```

| CHROM | n_variants |
|-------|------------|
| 1     | 1000       |
| 2     | 1000       |
| 3     | 1000       |

### Region Filtering

When you only need variants in a specific genomic region, the `region` parameter restricts processing **before** any computation begins. This is much more efficient than filtering with `WHERE` after the fact:

```sql
-- Efficient: only reads variants in the region
SELECT COUNT(*) AS n_variants
FROM plink_freq('test/data/large_example.pgen',
    region := '1:1-50000');
```

| n_variants |
|------------|
| 500        |

Region filtering works with all analysis functions (`plink_freq`, `plink_hardy`, `plink_missing`, `plink_ld`, `plink_score`) and `read_pfile`.

```sql
-- Combine region filtering with aggregation
SELECT COUNT(*) AS n_variants, ROUND(AVG(ALT_FREQ), 4) AS mean_freq
FROM plink_freq('test/data/large_example.pgen',
    region := '1:1-50000');
```

| n_variants | mean_freq |
|------------|-----------|
| 500        | 0.5       |

!!! tip "Performance tips"
    - Use `region` instead of `WHERE CHROM = ...` — it avoids reading genotype data outside the region
    - Select only the columns you need — [projection pushdown](guides/performance.md) skips unused computation
    - Use dedicated analysis functions (`plink_freq`, `plink_hardy`) instead of computing from raw genotypes — they use optimized pgenlib paths

    See the [Performance Guide](guides/performance.md) for the full set of optimization strategies.

---

## Next Steps

You've now seen the complete PlinkingDuck workflow: loading data, exploring genotypes, running QC, analyzing LD, computing PRS scores, and working efficiently at scale. Here are some directions to explore next:

- **[Function Reference](functions/index.md)** — full parameter and output documentation for every function
- **[Quality Control Guide](guides/quality-control.md)** — detailed QC workflows with recommended thresholds
- **[Performance Guide](guides/performance.md)** — optimization strategies for large datasets
- **[Common Parameters](common-parameters.md)** — sample subsetting, region filtering, and companion file options
- **[Genotype Encoding](genotype-encoding.md)** — how genotype values map to allele counts
