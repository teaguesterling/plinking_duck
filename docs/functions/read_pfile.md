# read_pfile

Read a complete PLINK 2 fileset (`.pgen` + `.pvar` + `.psam`) from a common prefix. Supports multiple orient modes, sample subsetting, variant filtering, and region filtering.

## Synopsis

```sql
read_pfile(prefix VARCHAR | LIST(VARCHAR) [, pgen := ..., pvar := ..., psam := ...,
           orient := ..., genotypes := ..., phased := ..., dosages := ...,
           samples := ..., variants := ..., region := ...,
           af_range := ..., ac_range := ...,
           include_genotypes := ..., genotype_range := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `prefix` | `VARCHAR` or `LIST(VARCHAR)` | *(required)* | Common file prefix (e.g., `'data/cohort'`), or a list of prefixes to read as one table (see [Multi-File Input](#multi-file-input)) |
| `pgen` | `VARCHAR` | `prefix.pgen` | Explicit `.pgen` path |
| `pvar` | `VARCHAR` | `prefix.pvar` | Explicit `.pvar` or `.bim` path |
| `psam` | `VARCHAR` | `prefix.psam` | Explicit `.psam` or `.fam` path |
| `orient` | `VARCHAR` | `'variant'` | Output orientation: `'variant'`, `'genotype'`, or `'sample'` |
| `combine_samples` | `VARCHAR` | `'implicit'` | Multi-file only: how shards' samples combine — `'implicit'` (positional trust) or `'identical'` (verify IIDs). See [Multi-File Input](#multi-file-input) |
| `genotypes` | `VARCHAR` | `'array'` | Genotype output: `'array'`, `'list'`, `'columns'`, `'struct'`, or the aggregate modes `'counts'` / `'stats'` (per-sample summaries in sample orient) |
| `phased` | `BOOLEAN` | `false` | Output phased haplotype pairs instead of dosage counts |
| `dosages` | `BOOLEAN` | `false` | Output dosage values instead of hard calls |
| `samples` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All | Filter to specific samples |
| `variants` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All | Filter to specific variants |
| `region` | `VARCHAR` | All | Filter to genomic region (`chr:start-end`) |
| `af_range` | `STRUCT(min DOUBLE, max DOUBLE)` | All | Filter variants by allele frequency range |
| `ac_range` | `STRUCT(min INTEGER, max INTEGER)` | All | Filter variants by allele count range |
| `include_genotypes` | `LIST(VARCHAR)` | All | Filter by hardcall category: `'hom_ref'`, `'het'`, `'hom_alt'`, `'missing'` (any subset; canonical genotype filter) |
| `genotype_range` | `STRUCT(min TINYINT, max TINYINT [, include_missing BOOLEAN])` | All | Numeric alias of `include_genotypes` for contiguous ranges over [0, 2] |

See [Common Parameters](../common-parameters.md) for details on `samples`, `region`, and filter parameters.

### Variant Filtering

The `variants` parameter accepts:

- **`LIST(VARCHAR)`** -- variant IDs (matched against ID in `.pvar`)
- **`LIST(INTEGER)`** -- 0-based variant indices

When both `region` and `variants` are specified, the result is the intersection (variants must match both filters).

## Output Columns

### Variant Orient Mode (Default)

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | `VARCHAR` | Chromosome |
| `POS` | `INTEGER` | Base-pair position |
| `ID` | `VARCHAR` | Variant identifier |
| `REF` | `VARCHAR` | Reference allele |
| `ALT` | `VARCHAR` | Alternate allele |
| `genotypes` | `ARRAY(TINYINT, N)` | Genotype calls, one per sample |

### Genotype Orient Mode (`orient := 'genotype'`)

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | `VARCHAR` | Chromosome |
| `POS` | `INTEGER` | Base-pair position |
| `ID` | `VARCHAR` | Variant identifier |
| `REF` | `VARCHAR` | Reference allele |
| `ALT` | `VARCHAR` | Alternate allele |
| `FID` | `VARCHAR` | Family ID (from `.psam`) |
| `IID` | `VARCHAR` | Individual ID (from `.psam`) |
| `SEX` | `INTEGER` | Sex (from `.psam`) |
| *additional* | varies | Any extra `.psam` columns |
| `genotype` | `TINYINT` | Single genotype value |

In genotype orient mode, each row represents one variant-sample combination. The sample metadata columns come from the `.psam` file and vary depending on what columns are present.

### Sample Orient Mode (`orient := 'sample'`)

| Column | Type | Description |
|--------|------|-------------|
| `FID` | `VARCHAR` | Family ID (from `.psam`) |
| `IID` | `VARCHAR` | Individual ID (from `.psam`) |
| `SEX` | `INTEGER` | Sex (from `.psam`) |
| *additional* | varies | Any extra `.psam` columns |
| `genotypes` | `ARRAY(TINYINT, M)` | Genotype calls, one per variant |

In sample orient mode, each row represents one sample. The genotypes array contains one value per variant (M = effective variant count after filtering). A multi-file list concatenates every shard's variants along this axis (M = total across shards).

**Aggregate output** (`genotypes := 'counts'` or `'stats'`) replaces the array with a per-sample summary STRUCT — how many `hom_ref` / `het` / `hom_alt` / `missing` genotypes that sample has across the variant range (`stats` adds `n`, `af`, `maf`, `missing_rate`, `carrier_count`, `het_rate`). These modes **stream** (they do not materialize the variants × samples matrix), so they are not bounded by `plinking_max_matrix_elements` and are the fast path for per-individual **carrier counting** over large variant ranges. See [Performance](../guides/optimizations.md#per-sample-counts-carrier-finding).

The array/list/struct/columns modes pre-read the genotype matrix at bind time (bounded by `plinking_max_matrix_elements`, default 16 G elements).

## Description

`read_pfile` is the most feature-rich reader function. It combines genotype, variant, and sample data with flexible filtering.

### File Discovery

By default, `read_pfile` constructs file paths by appending extensions to the prefix: `prefix.pgen`, `prefix.pvar`, `prefix.psam`. Each can be overridden individually.

Relative prefixes honor DuckDB's **`file_search_path`** setting (like `read_csv`): `SET file_search_path='/data/cohort'; SELECT * FROM read_pfile('chr22')` resolves `chr22.*` under that directory. Inputs containing a **glob** (`chr*.pgen`) or a registered **protocol** (`pathmacro:...`) are expanded via the VFS to concrete local paths before opening — a single glob/URL can fan out to a whole sharded fileset.

### Remote / cloud reads (`s3://`, `http`, …)

`read_pfile` (and `read_pgen` and every `plink_*` analysis function) can read a
`.pgen` **directly from any DuckDB filesystem** — `s3://`, `https://`, plus the
`.pvar`/`.psam` companions — with the `httpfs` extension loaded:

```sql
LOAD httpfs;
SELECT IID FROM read_pfile('s3://bucket/cohort', orient := 'sample',
    genotypes := 'counts', region := '17:43000000-43125000',
    include_genotypes := ['het','hom_alt']);   -- carriers in BRCA1, over S3
```

Credentials/secrets and settings are taken from the query context automatically.
A **targeted query fetches only the bytes it needs** (the variant index + the
region's records, via HTTP range reads) — not the whole file — so carrier/range
queries over a large remote `.pgen` are efficient. `.pgen` I/O is governed by
`plinking_pgen_io`:

| `plinking_pgen_io` | behavior |
|---|---|
| `'auto'` *(default)* | remote/VFS paths read through the VFS; local paths use native `fopen` (zero overhead) |
| `'native'` | always native `fopen` — errors on a remote path |
| `'vfs'` | always through the VFS (even local) — for testing |
| `'localize'` | *reserved* — download remote to a temp then read locally (not yet implemented) |

**Caveats:** a **full scan** of a remote `.pgen` re-reads more than a targeted
query (and each parallel thread reads the index + its range independently) —
load the community `cache_httpfs` extension (a block cache over `httpfs`) or
reduce threads for full-scan-over-remote workloads. Split-index (`.pgi`) filesets
are not yet supported (embedded-index `.pgen` only).

See the [File Handling](../guides/file-handling.md) guide for the full treatment of
discovery, companions, multi-file/sharded reads, path resolution, and remote reads.

### Multi-File Input

The first argument may be a single prefix or a `LIST(VARCHAR)` of prefixes, mirroring the way `read_csv` and `read_json` accept a list. A list **row-concatenates** the variants from each file, in file order, into one logical table:

```sql
SELECT * FROM read_pfile(['/data/chr1', '/data/chr2', '/data/chr3']);
```

This targets the biobank-scale layout where a fileset is **sharded by variant** — many `.pgen`/`.pvar` files that all share one identical `.psam`. All files in a list must carry the **same sample set: the same IIDs in the same order**. Alignment is **positional** and trust-the-caller — there is no per-file sample-identity check by design (adding one would reintroduce the `.psam`-parse cost on the real workload). The only guard is a free one: each `.pgen` header reports its sample count, and a file whose sample count differs from the first is rejected as a hard error. Sample order is *not* verified; that remains the caller's contract.

Metadata is taken from the first file only (identical by contract). Row order across files is not guaranteed (same as single-file parallel scan) — add `ORDER BY` if you need a stable order.

**All orient modes** accept a multi-file list. `variant` and `genotype` orient row-concatenate variants at scan time; `orient := 'sample'` concatenates every shard's variants into one genotypes array per sample (array dimension = total variants across all shards).

**`combine_samples`** controls how the shards' sample sets are combined:

- `'implicit'` *(default)* — trust positional alignment; no per-shard IID check, no extra I/O. This is the contract above.
- `'identical'` — verify every shard carries the **same IIDs in the same order** (one `.psam` parse per shard); a mismatch is a hard error naming the differing sample. Use it when you aren't certain the shards are aligned.
- `'union'` / `'intersect'` / `'concatenate'` — reserved for real sample-set joins; **not yet implemented** (they error).

A **`psam` override is allowed** with a multi-file list — it names the single shared `.psam` used for every shard (the common layout of per-shard `.pgen`/`.pvar` with one `.psam` kept elsewhere). `pgen`/`pvar` overrides remain single-file only (they can't disambiguate per shard). The `.psam` sample count must still match every shard's `.pgen`.

`region :=` works across a list (each shard is filtered to the range; the union is returned). `variants := [...]` is **not yet supported** with a multi-file list (IDs/indices don't resolve globally across shards yet) — use `region :=`.

The list of prefixes can come from anywhere: a literal list, a **glob** (`read_pfile('data/chr*.pgen')`), or a **catalog macro** returning a `VARCHAR[]` of shard prefixes — e.g. the scalarfs `pathmacro:` protocol, which `read_pfile` resolves to local shard paths and reads directly: `read_pfile('pathmacro:cohort?gene=BRCA1')`.

### Genotype Output Formats

The `genotypes` parameter controls how genotype data is structured:

- **`'array'`** (default): fixed-size `ARRAY(TINYINT, N)` — fastest, requires compile-time size
- **`'list'`**: variable-length `LIST(TINYINT)` — use for large cohorts (>100K samples)
- **`'columns'`**: one scalar `TINYINT` column per sample (variant orient) or per variant (sample orient); column names from IIDs or variant IDs. Incompatible with `orient := 'genotype'`.

### Phased Output

With `phased := true`, genotype elements change from `TINYINT` to `ARRAY(TINYINT, 2)` representing `[allele1, allele2]` haplotype pairs. `[0, 0]` = hom ref, `[0, 1]` = REF|ALT, `[1, 0]` = ALT|REF, `[1, 1]` = hom alt, NULL = missing. Incompatible with `dosages := true`. See [Genotype Encoding](../genotype-encoding.md).

### Filter Pushdown

`af_range` and `ac_range` filter variants by allele frequency or count using fast genotype counting (no decompression). `include_genotypes` (and its numeric alias `genotype_range`) filters by hardcall category — in `variant` orient it sets non-matching values to NULL; in `genotype` and `sample` orient it drops non-matching rows, so a carrier query in `sample` orient materializes only the matching subjects. See [Common Parameters](../common-parameters.md).

### Projection Pushdown

Genotype decoding is skipped when genotype columns (`genotypes` or `genotype`) are not referenced in the query.

### Parallel Scan

Default mode uses multi-threaded parallel scanning. Genotype orient mode is single-threaded because it expands each variant into multiple rows (one per sample). Sample orient mode is single-threaded (all genotypes are pre-read at bind time).

## Examples

```sql
-- Read all variants in default mode
SELECT * FROM read_pfile('data/example');
```

```sql
-- Multi-file: read several variant-sharded filesets (identical .psam) as one table
SELECT COUNT(*) FROM read_pfile(['data/chr1', 'data/chr2', 'data/chr3']);
```

```sql
-- Genotype orient mode: one row per variant-sample pair
SELECT chrom, pos, iid, genotype
FROM read_pfile('data/example', orient := 'genotype');
```

```sql
-- Sample orient: one row per sample
SELECT IID, genotypes
FROM read_pfile('data/example', orient := 'sample');
```

```sql
-- Sample subset by name
SELECT * FROM read_pfile('data/example',
    orient := 'genotype',
    samples := ['SAMPLE1', 'SAMPLE3']);
```

```sql
-- Region filter
SELECT * FROM read_pfile('data/example',
    region := '22:1-50000000');
```

```sql
-- Variant filter by ID
SELECT * FROM read_pfile('data/example',
    orient := 'genotype',
    variants := ['rs1', 'rs2']);
```

```sql
-- Combined filters
SELECT iid, genotype
FROM read_pfile('data/example',
    orient := 'genotype',
    region := '1:10000-50000',
    samples := ['SAMPLE1']);
```

```sql
-- Phased haplotype output
SELECT ID, genotypes
FROM read_pfile('data/example', phased := true);
```

```sql
-- Filter by allele frequency (common variants only)
SELECT * FROM read_pfile('data/example',
    af_range := {min: 0.01, max: 0.5});
```

```sql
-- Columns pivot mode
SELECT * FROM read_pfile('data/example',
    genotypes := 'columns');
```

### Multi-file, sharded filesets

```sql
-- Sample orient across a variant-sharded fileset: each sample's genotypes
-- array spans every shard's variants (one row per sample).
SELECT IID, genotypes
FROM read_pfile(['data/chr1', 'data/chr2', 'data/chr3'], orient := 'sample');
```

```sql
-- A glob fans out to every matching shard.
SELECT count(*) FROM read_pfile('data/chr*.pgen');
```

```sql
-- Verify the shards really share the same samples (else hard error).
SELECT count(*) FROM read_pfile(['data/chr1', 'data/chr2'],
    orient := 'sample', combine_samples := 'identical');
```

```sql
-- Per-shard .pgen/.pvar with one shared .psam kept elsewhere.
SELECT count(*) FROM read_pfile(['data/chr1', 'data/chr2'],
    psam := 'data/cohort.psam');
```

```sql
-- Resolve relative prefixes via file_search_path (like read_csv).
SET file_search_path = '/data/cohort';
SELECT count(*) FROM read_pfile(['chr1', 'chr2', 'chr3']);
```

### Per-individual carrier counting

```sql
-- How many rare-variant carriers does each subject have in a region?
-- Streams per-sample counts (no genotype matrix), scales to millions of samples.
SELECT IID, genotypes.het + genotypes.hom_alt AS n_carrier_variants
FROM read_pfile('data/cohort', orient := 'sample',
    region := '17:43000000-43125000',      -- e.g. BRCA1
    genotypes := 'counts')
WHERE genotypes.het + genotypes.hom_alt > 0
ORDER BY n_carrier_variants DESC;
```

```sql
-- Same, but let the row filter drop non-carriers for you.
SELECT IID FROM read_pfile('data/cohort', orient := 'sample',
    genotypes := 'counts',
    include_genotypes := ['het', 'hom_alt']);   -- keeps only carriers

-- For biobank-scale rare-variant ranges, enable the sparse (difflist) path:
SET plinking_sample_counts_sparse = true;       -- 4-6x faster on rare variants
```

## See Also

- [read_pgen](read_pgen.md) -- lower-level genotype reader
- [read_pvar](read_pvar.md) -- variant metadata only
- [read_psam](read_psam.md) -- sample metadata only
- [Genotype Encoding](../genotype-encoding.md) -- encoding and phased output reference
