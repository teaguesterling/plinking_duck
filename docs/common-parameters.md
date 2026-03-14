# Common Parameters

Several parameters are shared across multiple PlinkingDuck functions. This page documents their behavior in one place.

## `pvar`

| Type | Default |
|------|---------|
| `VARCHAR` | Auto-discovered |

Explicit path to a `.pvar` or `.bim` file containing variant metadata. When omitted, the function replaces the `.pgen` extension with `.pvar` (then `.bim`) and uses the first file that exists.

```sql
-- Auto-discover (looks for data.pvar, then data.bim)
SELECT * FROM plink_freq('data.pgen');

-- Explicit path
SELECT * FROM plink_freq('data.pgen', pvar := '/other/path/variants.pvar');
```

**Used by:** `read_pgen`, `read_pfile`, `plink_freq`, `plink_hardy`, `plink_missing`, `plink_ld`, `plink_score`, `plink_glm`

## `psam`

| Type | Default |
|------|---------|
| `VARCHAR` | Auto-discovered |

Explicit path to a `.psam` or `.fam` file containing sample metadata. When omitted, the function replaces the `.pgen` extension with `.psam` (then `.fam`) and uses the first file that exists.

The `.psam` file is **optional** for most functions -- they can operate in index-only mode using just the sample count from the `.pgen` header. It is **required** for:

- `read_pfile` (needs sample IDs for genotype orient mode output)
- `plink_score` (needs sample IDs for output)
- `plink_missing` in sample mode (needs sample IDs for output)
- Any function when `samples` uses `LIST(VARCHAR)` (needs IDs to resolve names)

```sql
-- Explicit sample file
SELECT * FROM read_pgen('data.pgen', psam := '/other/path/samples.fam');
```

**Used by:** `read_pgen`, `read_pfile`, `plink_freq`, `plink_hardy`, `plink_missing`, `plink_ld`, `plink_score`, `plink_glm`

## `samples`

| Type | Default |
|------|---------|
| `LIST(VARCHAR)` or `LIST(INTEGER)` | All samples |

Subset the analysis to specific samples. Accepts either:

- **`LIST(VARCHAR)`** -- sample IDs (matched against IID in the `.psam` file)
- **`LIST(INTEGER)`** -- 0-based sample indices into the `.pgen` file order

```sql
-- By sample name (requires .psam)
SELECT * FROM plink_freq('data.pgen', samples := ['SAMPLE1', 'SAMPLE3']);

-- By 0-based index (works without .psam)
SELECT * FROM plink_freq('data.pgen', samples := [0, 2, 5]);
```

The list must not be empty and must not contain duplicates.

When no `.psam` file is available, only `LIST(INTEGER)` is accepted.

**Used by:** `read_pgen`, `read_pfile`, `plink_freq`, `plink_hardy`, `plink_missing`, `plink_ld`, `plink_score`, `plink_glm`

## `region`

| Type | Default |
|------|---------|
| `VARCHAR` | All variants |

Restrict the analysis to variants within a genomic region. Format: `chr:start-end` where positions are 1-based and inclusive.

```sql
-- Variants on chromosome 22 between positions 1 and 50,000,000
SELECT * FROM plink_freq('data.pgen', region := '22:1-50000000');
```

The function assumes variants are sorted by (CHROM, POS) as required by the PLINK format specification. The region filter finds all variants matching the chromosome and falling within the position range.

Both `start` and `end` must be specified (chromosome-only filtering is not supported in this parameter; use a `WHERE` clause instead).

**Used by:** `read_pfile`, `plink_freq`, `plink_hardy`, `plink_missing`, `plink_ld`, `plink_score`, `plink_glm`

!!! note
    `read_pgen` and `read_pvar` do not have a `region` parameter. Use a SQL `WHERE` clause to filter their output by position.

## `af_range`

| Type | Default |
|------|---------|
| `STRUCT(min DOUBLE, max DOUBLE)` | All variants |

Filter variants by alternate allele frequency. Uses `PgrGetCounts` for fast pre-filtering without genotype decompression. The frequency denominator is `2 * non_missing_samples`.

```sql
-- Keep only common variants (MAF > 1%)
SELECT * FROM read_pfile('data', af_range := {min: 0.01, max: 0.99});

-- Rare variants only
SELECT * FROM plink_freq('data.pgen', af_range := {min: 0.0, max: 0.01});
```

Can be combined with `ac_range`.

**Used by:** `read_pgen`, `read_pfile`

## `ac_range`

| Type | Default |
|------|---------|
| `STRUCT(min INTEGER, max INTEGER)` | All variants |

Filter variants by alternate allele count. Uses `PgrGetCounts` for fast pre-filtering.

```sql
-- Variants with at least 2 alt alleles
SELECT * FROM read_pfile('data', ac_range := {min: 2, max: 1000000});
```

Can be combined with `af_range`.

**Used by:** `read_pgen`, `read_pfile`

## `genotype_range`

| Type | Default |
|------|---------|
| `STRUCT(min TINYINT, max TINYINT)` | All genotype values |

Filter individual genotype values. Genotypes outside the range are set to NULL. The domain is [0, 2]. Uses `PgrGetCounts` for pre-decompression optimization: variants where no values fall in range are skipped entirely.

```sql
-- Only heterozygous calls (genotype = 1)
SELECT * FROM read_pgen('data.pgen', genotype_range := {min: 1, max: 1});

-- Het and hom-alt calls only
SELECT * FROM read_pfile('data', genotype_range := {min: 1, max: 2});
```

Incompatible with `dosages := true`.

**Used by:** `read_pgen`, `read_pfile`
