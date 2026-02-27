# plink_ld

Compute pairwise linkage disequilibrium statistics (r², D').

## Synopsis

```sql
-- Pairwise mode: LD between two specific variants
plink_ld(path VARCHAR, variant1 := ..., variant2 := ...) -> TABLE

-- Windowed mode: LD for all variant pairs within a window
plink_ld(path VARCHAR [, window_kb := ..., r2_threshold := ...,
         region := ..., samples := ..., inter_chr := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `path` | `VARCHAR` | *(required)* | Path to the `.pgen` file |
| `pvar` | `VARCHAR` | Auto-discovered | Path to `.pvar` or `.bim` file |
| `psam` | `VARCHAR` | Auto-discovered | Path to `.psam` or `.fam` file |
| `variant1` | `VARCHAR` | *(none)* | First variant ID (pairwise mode) |
| `variant2` | `VARCHAR` | *(none)* | Second variant ID (pairwise mode) |
| `window_kb` | `INTEGER` | `1000` | Window size in kilobases (windowed mode) |
| `r2_threshold` | `DOUBLE` | `0.2` | Minimum r² to report (windowed mode) |
| `samples` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All | Subset to specific samples |
| `region` | `VARCHAR` | All | Filter to genomic region (`chr:start-end`) |
| `inter_chr` | `BOOLEAN` | `false` | Include cross-chromosome pairs (windowed mode) |

See [Common Parameters](../common-parameters.md) for details on `pvar`, `psam`, `samples`, and `region`.

## Modes

### Pairwise Mode

When both `variant1` and `variant2` are specified, the function computes LD for that single pair and returns exactly one row. Both variant IDs must exist in the `.pvar` file.

### Windowed Mode

When neither `variant1` nor `variant2` is specified, the function scans all variant pairs within a sliding window. For each anchor variant, it tests all subsequent variants within `window_kb` kilobases on the same chromosome and reports pairs where r² >= `r2_threshold`.

With `inter_chr := true`, cross-chromosome pairs are also tested (no distance filter applies across chromosomes).

## Output Columns

| Column | Type | Description |
|--------|------|-------------|
| `CHROM_A` | `VARCHAR` | Chromosome of first variant |
| `POS_A` | `INTEGER` | Position of first variant |
| `ID_A` | `VARCHAR` | ID of first variant |
| `CHROM_B` | `VARCHAR` | Chromosome of second variant |
| `POS_B` | `INTEGER` | Position of second variant |
| `ID_B` | `VARCHAR` | ID of second variant |
| `R2` | `DOUBLE` | Squared correlation coefficient (NULL if invalid) |
| `D_PRIME` | `DOUBLE` | Normalized linkage disequilibrium coefficient (NULL if invalid) |
| `OBS_CT` | `INTEGER` | Number of samples with non-missing genotypes at both variants |

R2 and D_PRIME are NULL when the LD computation is invalid (monomorphic variants, fewer than 2 shared observations).

## Description

### LD Statistics

**r²** is the squared Pearson correlation between genotype dosages at two loci. Values range from 0 (no LD) to 1 (perfect LD).

**D'** is computed using the composite LD estimator (Weir 1979) from genotype-level statistics. D is the covariance between genotype dosages divided by 4, and D' = D / D_max where D_max depends on allele frequencies and the sign of D.

!!! note
    Because this uses genotype-level (not haplotype-level) statistics, D' can exceed 1.0 when samples deviate from Hardy-Weinberg equilibrium.

### Invalid Results

LD statistics are NULL when:

- Either variant is monomorphic (zero variance)
- Fewer than 2 samples have non-missing genotypes at both loci

### Parallel Scan

Windowed mode supports multi-threaded scanning where each thread claims anchor variants independently.

## Examples

```sql
-- Pairwise LD between two variants
SELECT * FROM plink_ld('data/example.pgen',
    variant1 := 'rs1', variant2 := 'rs2');
```

```sql
-- Windowed LD with default parameters (1000kb window, r² >= 0.2)
SELECT * FROM plink_ld('data/example.pgen');
```

```sql
-- Narrow window with lower threshold
SELECT * FROM plink_ld('data/example.pgen',
    window_kb := 250, r2_threshold := 0.05);
```

```sql
-- LD in a specific region
SELECT * FROM plink_ld('data/example.pgen',
    region := '1:1-1000000');
```

```sql
-- Include cross-chromosome LD
SELECT * FROM plink_ld('data/example.pgen',
    inter_chr := true, r2_threshold := 0.8);
```

## See Also

- [plink_freq](plink_freq.md) -- allele frequencies
- [plink_hardy](plink_hardy.md) -- HWE test
