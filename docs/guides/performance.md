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
| `read_plink_vcf` | No | Sequential text VCF parsing |
| `read_pvar` | No | Sequential text file reading |
| `read_psam` | No | Sequential text file reading |
| `read_pgen` | Yes | Per-thread pgenlib readers |
| `read_pfile` (default) | Yes | Same as `read_pgen` |
| `read_pfile` (genotype orient) | No | Row expansion is single-threaded |
| `read_pfile` (sample orient) | No | Pre-reads all genotypes at bind time |
| `plink_freq` | Yes | Per-thread readers, atomic batch claiming |
| `plink_hardy` | Yes | Same pattern as `plink_freq` |
| `plink_missing` (variant) | Yes | Parallel missingness extraction |
| `plink_missing` (sample) | No | Parallel per-sample accumulation |
| `plink_ld` (windowed) | Yes | Per-thread anchor claiming |
| `plink_ld` (pairwise) | No | Single pair computation |
| `plink_score` | Yes | Parallel per-sample accumulation |
| `plink_glm` | Yes | Per-thread regression |

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

## Filter Pushdown

The `af_range`, `ac_range`, `include_genotypes`, and `genotype_range` parameters let you filter variants or genotypes before full decompression:

```sql
-- Skip variants outside MAF range (uses PgrGetCounts, no decompression)
SELECT * FROM read_pfile('biobank',
    af_range := {min: 0.01, max: 0.5});

-- Only keep heterozygous calls
SELECT * FROM read_pgen('biobank.pgen',
    include_genotypes := ['het']);
```

`af_range` and `ac_range` use pgenlib's `PgrGetCounts` to count genotypes without decompressing, then skip variants that fail the filter. `include_genotypes` (and its numeric alias `genotype_range`) also uses `PgrGetCounts` as a pre-check — if a variant has no matching genotypes, it is skipped entirely without decompression.

### Fast carrier lookups (subject IDs matching a genotype)

To list the subjects matching a genotype at a variant, combine `orient := 'sample'`
with `include_genotypes` so the scan emits **only the matching subjects** instead
of one row per sample. At biobank scale (millions of individuals) this is the
difference between materializing every subject and materializing just the carriers:

```sql
-- Subjects carrying an ALT allele at rs12345 (het or hom-alt)
SELECT IID
FROM read_pfile('biobank', orient := 'sample', genotypes := 'counts',
                variants := ['rs12345'], include_genotypes := ['het', 'hom_alt']);
```

The row filter is evaluated during the scan, so only carrier rows are produced —
a `WHERE` clause on the output cannot achieve this because the per-subject counts
are computed inside the table function. The same works in `genotype` orient
(one row per matching variant × subject), and the filter applies whether or not
the genotype column itself is selected.

## Per-sample counts & carrier finding

`orient := 'sample'` with `genotypes := 'counts'` or `'stats'` produces a per-sample
summary across a variant range — for counting how many variants each subject carries.
Both modes take a **streaming** path: a variant-parallel pass accumulates per-sample
category counts (`hom_ref`/`het`/`hom_alt`/`missing`) directly, **without materializing
the variants × samples genotype matrix**. So:

- Memory is `O(samples)`, not `O(variants × samples)` — a 5 000-variant × 7 M-sample
  range is tens of MB, not a ~35 GB matrix. These modes are **not** bounded by
  `plinking_max_matrix_elements` (the array/list/struct/columns modes still are).
- It scales over threads (the accumulation is parallel over the variant range).
- `stats` uses the **same** accumulation as `counts` — its `n`/`af`/`maf`/
  `missing_rate`/`carrier_count`/`het_rate` fields are derived from the same
  per-sample counters at no extra decode cost.

```sql
-- Per-subject carrier burden over a gene, streamed:
SELECT IID, genotypes.het + genotypes.hom_alt AS n_carrier_variants
FROM read_pfile('biobank', orient := 'sample',
    region := '17:43000000-43125000', genotypes := 'counts')
WHERE genotypes.het + genotypes.hom_alt > 0;
```

### Sparse (difflist) path for rare variants

For rare-variant carrier finding, most subjects are `hom_ref` at any given variant.
The `.pgen` difflist stores only the non-`hom_ref` genotypes of a sparse variant, so
the counts can be accumulated by touching only the carriers — skipping the `hom_ref`
majority entirely. Enable it with:

```sql
SET plinking_sample_counts_sparse = true;   -- default false
```

Measured **~4× (moderate-rare, MAF < 1%)** to **~6× (ultra-rare, MAF < 0.05%)** faster
on 100 K samples × 30 K variants; the win grows with both rarity and sample count.
The path auto-falls-back to the dense loop for common variants, so it is never worse
on common data, and it produces **identical** results either way. It applies to
`counts` and `stats` (the streaming aggregate modes) in sample orient. Left off by
default pending broad validation — A/B it on your data and enable per session.

## Configuration

| Setting | Default | Effect |
|---------|---------|--------|
| `plinking_max_threads` | `0` (cap 16) | Cap threads for all parallel scans |
| `plinking_max_matrix_elements` | `16 G` | Ceiling for the `orient := 'sample'` genotype-matrix pre-read (array/list/struct/columns; **not** counts/stats, which stream) |
| `plinking_sample_counts_sparse` | `false` | Use the sparse difflist path for sample-orient `counts`/`stats` (see above) |
| `plinking_use_parquet_companions` | `true` | Prefer `.pvar.parquet` / `.psam.parquet` companions |

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
FROM read_pfile('data', orient := 'genotype')
GROUP BY ID;
```

### Subset samples when possible

If your analysis only involves a subset of samples, specify them upfront:

```sql
-- Efficient: pgenlib only processes 3 samples per variant
SELECT * FROM plink_freq('data.pgen',
    samples := ['CASE1', 'CASE2', 'CASE3']);
```
