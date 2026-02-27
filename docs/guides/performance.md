# Performance

PlinkingDuck is designed for efficient access to large PLINK datasets. This guide explains the performance features and how to take advantage of them.

## Projection Pushdown

All PlinkingDuck functions support **projection pushdown**: if a column is not referenced in your query, the work to compute it is skipped.

This is most impactful for genotype data. When you only need variant metadata, genotype decoding is skipped entirely:

```sql
-- Fast: no genotype decoding
SELECT CHROM, POS, ID FROM read_pgen('large_dataset.pgen');

-- Slower: requires decoding genotypes for every variant
SELECT * FROM read_pgen('large_dataset.pgen');
```

For analysis functions, projection pushdown also works:

```sql
-- Fast: skips PgrGetCounts, only reads variant metadata
SELECT CHROM, POS, ID FROM plink_freq('large_dataset.pgen');

-- Normal: computes frequencies
SELECT ID, ALT_FREQ FROM plink_freq('large_dataset.pgen');
```

## Parallel Scanning

Functions that process `.pgen` files use multi-threaded parallel scanning:

| Function | Parallel | Notes |
|----------|----------|-------|
| `read_pvar` | No | Sequential text file reading |
| `read_psam` | No | Sequential text file reading; `.psam` files are small |
| `read_pgen` | Yes | Per-thread pgenlib readers |
| `read_pfile` (default) | Yes | Same as `read_pgen` |
| `read_pfile` (tidy) | No | Row expansion is single-threaded |
| `plink_freq` | Yes | Per-thread readers, atomic batch claiming |
| `plink_hardy` | Yes | Same pattern as `plink_freq` |
| `plink_missing` (variant) | Yes | Parallel missingness extraction |
| `plink_missing` (sample) | No | Must scan all variants for each sample |
| `plink_ld` (windowed) | Yes | Per-thread anchor claiming |
| `plink_ld` (pairwise) | No | Single pair computation |
| `plink_score` | No | Sequential variant scoring |

Each parallel function uses atomic batch claiming: threads claim batches of variants from a shared counter, ensuring even work distribution without lock contention.

The maximum thread count scales with the number of variants (typically `min(variants/500 + 1, 16)`).

## Region Filtering

For large datasets, `region` filtering restricts processing to a genomic region before any computation begins:

```sql
-- Only processes variants on chr22 between positions 1-50M
SELECT * FROM plink_freq('biobank.pgen',
    region := '22:1-50000000');
```

Region filtering works by scanning the pre-loaded variant metadata (sorted by chromosome and position) to find the matching variant index range. The .pgen reader then starts at the first matching variant and stops after the last.

This avoids reading genotype data for variants outside the region.

## Sample Subsetting

When using the `samples` parameter, only the specified samples are processed:

```sql
SELECT * FROM plink_freq('biobank.pgen',
    samples := [0, 100, 200]);
```

Sample subsetting is handled at the pgenlib level using `PgrGetCounts` / `PgrGet` with a sample bitmask. For functions using `PgrGetCounts` (`plink_freq`, `plink_hardy`), this is very efficient since the counting path directly respects the bitmask without extracting individual genotypes.

## Tips for Large Datasets

### Use region filtering for targeted analyses

Instead of scanning all variants and filtering with `WHERE`:

```sql
-- Less efficient: scans all variants, filters in SQL
SELECT * FROM plink_freq('biobank.pgen')
WHERE CHROM = '22' AND POS BETWEEN 1 AND 50000000;
```

Use the `region` parameter:

```sql
-- More efficient: only reads variants in the region
SELECT * FROM plink_freq('biobank.pgen',
    region := '22:1-50000000');
```

### Select only the columns you need

Projection pushdown means unused columns are free to ignore:

```sql
-- Faster: only decodes what's needed
SELECT ID, ALT_FREQ FROM plink_freq('biobank.pgen');

-- Slower: computes everything including counts
SELECT * FROM plink_freq('biobank.pgen', counts := true);
```

### Use analysis functions instead of raw genotypes

When you need summary statistics, use the dedicated functions rather than computing from raw genotypes. They use optimized pgenlib paths:

```sql
-- Efficient: uses PgrGetCounts (no decompression)
SELECT ID, ALT_FREQ FROM plink_freq('data.pgen');

-- Less efficient: decompresses genotypes, then aggregates in SQL
SELECT ID, AVG(genotype) / 2.0 AS alt_freq
FROM read_pfile('data', tidy := true)
GROUP BY ID;
```

### Subset samples when possible

If your analysis only involves a subset of samples, specify them upfront:

```sql
-- Efficient: pgenlib only processes 3 samples per variant
SELECT * FROM plink_freq('data.pgen',
    samples := ['CASE1', 'CASE2', 'CASE3']);
```
