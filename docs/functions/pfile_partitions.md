# pfile_partitions

Compute the ordinal-range partitions of a PLINK 2 fileset along the **variants** dimension. Each partition is a half-open `[lo, hi)` range in variant ordinal space, suitable for driving chunked or parallel reads of `read_pfile` / `read_pvar` (pair it with `variant_range :=` on `read_pfile`).

## Synopsis

```sql
pfile_partitions(prefix VARCHAR | LIST(VARCHAR) [, chunk_size := N]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `prefix` | `VARCHAR` or `LIST(VARCHAR)` | *(required)* | Common file prefix (e.g., `'data/cohort'`), or a list of prefixes. A glob or registered protocol is expanded via the VFS, like `read_pfile` |
| `chunk_size` | `BIGINT` | *(derived)* | Partition size in variants (must be a positive integer). When omitted, a size is derived from the variant count (see [Default Chunk Size](#default-chunk-size)) |

Only the `.pgen` header is read — the variant count comes from `raw_variant_ct`, so no `.pvar`/`.psam` parse and no variant-record scan is performed. Remote/VFS `.pgen` paths are opened through the VFS (a tiny range read for the header).

## Output Columns

| Column | Type | Description |
|--------|------|-------------|
| `source` | `VARCHAR` | The resolved `.pgen` path this partition reads from (absolute/resolved, not the input prefix) |
| `dim` | `VARCHAR` | The partitioned dimension; always `'variants'` (the samples dimension is unchunked by design) |
| `partition_id` | `BIGINT` | 0-based partition index **within the source** |
| `lo` | `BIGINT` | Inclusive low bound, in variant ordinal space |
| `hi` | `BIGINT` | Exclusive high bound, in variant ordinal space |

`lo`/`hi` live in the **same 0-based variant ordinal space** as `read_pfile` / `read_pvar`'s emitted [`__ordinal`](read_pvar.md#named-parameters) (and as `read_pfile`'s [`variant_range :=`](read_pfile.md#variant-range-ordinal-space) parameter).

## Description

`pfile_partitions` tiles `[0, n_variants)` into disjoint, ordered, covering half-open ranges. For a single source the partitions satisfy:

- **disjoint** — no two ranges overlap;
- **ordered** — emitted in ascending `lo` (`partition_id` increasing);
- **covering** — together they cover exactly `[0, n_variants)`, the same count `read_pfile` uses for its ordinal space.

These guarantees are **per source**. For a multi-prefix list, partitions are emitted per prefix and both `partition_id` and the `[lo, hi)` ordinal space **restart at 0 for each source** (matching `read_pfile`'s per-source `__ordinal`). The `source` column is what lets a consumer map per-source partitions into a global chunk space when concatenating shards, and to read each shard separately with its own `variant_range`.

A zero-variant prefix contributes **no** partition rows (it covers `[0, 0)`).

### Default Chunk Size

When `chunk_size` is omitted, the size is `min(n_variants / 16, 50000)` rounded **half-up** to the nearest 1000, but never 0 — a file too small for that formula (fewer than 8,000 variants, where the rounded size would be 0) becomes a single partition covering `[0, n_variants)`.

> This derived default is specific to `pfile_partitions`. Consumers that mint their own chunks (e.g. howard's `register_pgen`) pass an explicit `chunk_size :=`, so they are never subject to this default — and their own chunk-size heuristic need not match it exactly.

## Examples

```sql
-- Default partitioning of a fileset
SELECT * FROM pfile_partitions('data/cohort');
```

```sql
-- Fixed partition size of 10,000 variants
SELECT partition_id, lo, hi
FROM pfile_partitions('data/cohort', chunk_size := 10000)
ORDER BY partition_id;
```

```sql
-- Drive a chunked read: one variant_range scan per partition
SELECT p.partition_id, v.*
FROM pfile_partitions('data/cohort', chunk_size := 10000) p,
     LATERAL read_pfile('data/cohort', variant_range := [p.lo, p.hi]) v;
```

```sql
-- Per-source partitions across a variant-sharded fileset (ordinals restart per source)
SELECT source, partition_id, lo, hi
FROM pfile_partitions(['data/chr1', 'data/chr2', 'data/chr3'])
ORDER BY source, partition_id;
```

With a 4-variant fixture, `pfile_partitions(prefix, chunk_size := 2)` yields two partitions — `(partition_id=0, lo=0, hi=2)` and `(partition_id=1, lo=2, hi=4)` — while the default (n too small for the formula) yields one partition `(0, 0, 4)`.

## See Also

- [read_pfile](read_pfile.md) — read a fileset; pair with `variant_range := [lo, hi]` to read one partition
- [read_pvar](read_pvar.md) — variant metadata with an emitted `__ordinal` in the same ordinal space
- [plinking_duck_capabilities](plinking_duck_capabilities.md) — detect whether a build provides `pfile_partitions`
- [Optimizations](../guides/optimizations.md#ordinal-range-pushdown-and-partitioning) — the ordinal-range pushdown this feeds
