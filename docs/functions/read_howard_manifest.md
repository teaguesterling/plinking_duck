# read_howard_manifest

Parse, verify, and prune a `.howard.json` shard catalog, returning **one row per surviving shard**.

This is the parse + verify + prune surface plinking_duck owns for `.howard` **reading**. The **write**
side is the sharded pfile writer (`COPY (…) TO 'out' (FORMAT pfile, shard_by '<col>')`, see
[Writing filesets](write_pfile.md)), which emits exactly the `out.howard.json` this function reads —
closing the write→read loop. The format is howard's `.howard` catalog design (variant-shards over a
shared `.psam`, or sample-shards over a shared `.pvar`); this function is Phase A of consuming it and
does **not** itself scan genotypes or `ATTACH` — it resolves which shards a query needs and vouches
for them.

For the catalog format itself — the JSON schema, the boundary law, and the trust model — see the
[**`.howard` catalog format reference**](../formats/howard.md).

## Synopsis

```sql
read_howard_manifest(path VARCHAR
    [, region := VARCHAR,
       variant_index := ANY,
       sample_index := ANY,
       verify_all := BOOLEAN]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `path` | `VARCHAR` | *(required)* | Path to a multi-file `.howard.json` catalog. Every path inside it resolves **relative to the catalog's own directory** (an absolute path is honored as-is). |
| `region` | `VARCHAR` | *(none)* | Genomic prune, `'contig:lo-hi'` or a bare `'contig'`. Keeps only shards whose declared `range` could overlap. Valid **only** when `concat_dim = 'variants'`. POS bounds are **inclusive** `[lo, hi]`. |
| `variant_index` | `ANY` (2-element list `[lo, hi]`) | *(none)* | Index prune on the variants axis. Keeps shards whose `partition` overlaps the **inclusive** index predicate `[lo, hi]`. Valid **only** when `concat_dim = 'variants'`. |
| `sample_index` | `ANY` (2-element list `[lo, hi]`) | *(none)* | Index prune on the samples axis, inclusive `[lo, hi]`. Valid **only** when `concat_dim = 'samples'`. |
| `verify_all` | `BOOLEAN` | `false` | When `true`, verify **every** shard (counts + declared range) *before* pruning, closing the dropped-shard hole; otherwise only the survivors are verified. |

!!! warning "`variant_index` is INCLUSIVE `[lo, hi]` — not `read_pfile`'s `variant_range`"
    This function's index param is named `variant_index` (and `sample_index`), deliberately distinct
    from `read_pfile`'s `variant_range := [lo, hi)`, which is **half-open**. `variant_index` /
    `sample_index` are **inclusive** `[lo, hi]`, to match `WHERE variant_index BETWEEN lo AND hi` and
    howard's phase-3 index prune. A shard's `partition`
    is half-open `[start, stop)`, so a shard (holding ordinals `start..stop-1`) survives iff
    `start <= hi AND stop > lo`. Example: `[1000, 1000]` selects only the shard whose partition
    contains ordinal 1000 (`[1000, 2000)`), not `[0, 1000)`.

## Pruning model

Pruning is on **monotonic** coordinates only, and a predicate is valid on exactly one axis (an
off-axis predicate is an error, never a silent whole-cohort return):

- **Index prune** (`variant_index` / `sample_index`) reads the existing `partition` ordinal range —
  no `range` entry needed, and it is the prune key on **both** axes.
- **Genomic prune** (`region`, variants only) reads the per-shard `range`: a shard is a candidate iff
  the contig is in its `range.CHROM` list **and** the query POS window intersects its inclusive
  `range.POS` extent. The test is **conservative**: a shard is dropped only when its declared `range`
  *proves* disjointness. A shard with **no `range`** (or a partial one) is **always** a candidate —
  it is never wrongly dropped.

!!! note "Survivor ordinals stay global — no rebasing"

    When a prune drops shards, the surviving rows keep their **original global**
    `partition_start` / `partition_stop` and `source_ordinal` — this function does **not** rebase the
    subset to `[0, subset_total)`. (The Python `howard` loader rebases a pruned sub-cohort to a
    self-contained `[0, …)` index and warns when it does; the SQL manifest does not.) So after a prune
    the returned ordinals are still positions in the *full* catalog, not a renumbered subset.

## Single-fileset catalogs

A `.pfile.json` / single-`pfile` catalog (one fileset, no sharding) returns **one row** describing the
whole fileset: `concat_dim` is `NULL`, `partition_start = 0`, `partition_stop = <probed variant
count>`, `source_ordinal = 0`, resolved `prefix`/`pvar`/`psam`, and `chrom`/`pos_min`/`pos_max` read
from the single `.pvar`. A single fileset is always the whole candidate, so a prune predicate is
accepted and ignored (whole-load) rather than raising a wrong-axis error. A `chunks` object, if
present, is validated (`dim` must be `variants`; `chunk_size` a positive integer); its `chunk_size`
is surfaced in the `chunk_size` column (NULL on multi-shard rows).

An **omitted** companion (`pvar`/`psam` not declared, on any catalog form) resolves to `NULL` in the
output column, and verification lets `read_pfile`'s own discovery find it — including a
`.pvar.parquet` / `.psam.parquet` sidecar or a legacy `.bim` / `.fam`.

## Verification (the trust model)

`.howard` bounds are *assertions about foreign files plinking_duck did not write*, so they are checked
eagerly at open rather than trusted:

- **Count check** — each surviving shard's `.pgen` header `raw_variant_ct` / `raw_sample_ct` is
  asserted against its companion `.pvar` / `.psam` row counts (the shared companion for the shared
  axis), **and** the declared `partition` width (`partition_stop - partition_start`) must equal the
  real companion row count. Any mismatch raises — so `verified` is never a lie.
- **Range check** (variants only) — a surviving shard's declared `range` must **cover** its real
  `.pvar` extent (declared ⊇ actual): every contig the `.pvar` holds must be in `range.CHROM`, and
  the real POS extent must lie within `range.POS`. An **under-claim** raises (it would wrongly drop
  the shard); an **over-claim** is harmless. A shard with no `range` is skipped.

By default only the **survivors** of a prune are verified — a staleness detector for the drops
actually made. It cannot see an under-claim on a shard that was *itself* dropped (never opened).
`verify_all := true` verifies **every** shard before pruning, which closes that hole (at the cost of
reading every `.pvar`); use it for a hand-written catalog whose `range` values are not
writer-constructed.

## Output Columns

| Column | Type | Description |
|--------|------|-------------|
| `dataset` | `VARCHAR` | The catalog's `dataset` name. |
| `concat_dim` | `VARCHAR` | `'variants'` or `'samples'` — the axis the shards concatenate on; **`NULL` for a single-fileset catalog** (the signal a consumer routes a whole-fileset load on). |
| `prefix` | `VARCHAR` | Resolved fileset path (no extension). |
| `pvar` | `VARCHAR` | Resolved `.pvar` (the shard's own for variant-shards; the shared one for sample-shards). |
| `psam` | `VARCHAR` | Resolved `.psam` (the shared one for variant-shards; the shard's own for sample-shards). |
| `partition_start` | `BIGINT` | Global half-open `[start, stop)` ordinal on the concat axis — included. |
| `partition_stop` | `BIGINT` | Global half-open ordinal — excluded. |
| `source_ordinal` | `BIGINT` | The shard's descriptive source ordinal. |
| `chrom` | `LIST(VARCHAR)` | `range.CHROM` (the contigs the shard covers); `NULL` when the shard declares no `range`. |
| `pos_min` | `BIGINT` | Inclusive `range.POS.min`; `NULL` when absent. |
| `pos_max` | `BIGINT` | Inclusive `range.POS.max`; `NULL` when absent. |
| `verified` | `BOOLEAN` | Whether this shard's companions were checked during this call. Currently always `true` on every returned row (a verification failure raises); the column exists for a future degrade-gracefully mode. |
| `chunk_size` | `BIGINT` | A **single-fileset** catalog's declared `chunks.chunk_size` (consumed by howard's single-file loader). `NULL` for multi-shard rows — there the shards themselves are the chunks. |

## Error taxonomy

Errors carry distinct, classifiable message prefixes so a caller can tell failure classes apart:

| Prefix | Class | Exception |
|--------|-------|-----------|
| `read_howard_manifest: catalog parse error: …` | Malformed JSON or schema (bad version, missing field, both `pfile` and multi-file keys, duplicate `partition.start`, out-of-range integer) | `InvalidInputException` |
| `read_howard_manifest: file not found: …` | A referenced catalog / `.pgen` / `.pvar` / `.psam` is missing | `IOException` |
| `read_howard_manifest: verification failed (count mismatch): …` | A `.pgen` header count disagrees with its companion | `InvalidInputException` |
| `read_howard_manifest: verification failed (companion unreadable): …` | A companion exists but could not be read (corrupt `.pgen` header, wrong-format `.pvar`/`.psam`) | `InvalidInputException` |
| `read_howard_manifest: verification failed (range under-claim): …` | A declared `range` does not cover the real `.pvar` extent | `InvalidInputException` |
| `read_howard_manifest: invalid predicate: …` | Wrong-axis or malformed prune argument | `InvalidInputException` |

## Examples

```sql
-- One row per shard, with each shard's resolved paths, partition tiling, and range extent.
SELECT * FROM read_howard_manifest('cohort.howard.json');

-- Genomic prune: only the shards whose range overlaps chr2:1-5Mb survive (and are verified).
SELECT prefix, pvar FROM read_howard_manifest('cohort.howard.json', region := '2:1-5000000');

-- Index prune on the variants axis (inclusive): shards touching variant ordinals 1000..2000.
SELECT prefix FROM read_howard_manifest('cohort.howard.json', variant_index := [1000, 2000]);

-- Strict mode for a hand-written catalog: verify every shard's range before pruning.
SELECT count(*) FROM read_howard_manifest('cohort.howard.json', region := '7:1-1000000', verify_all := true);
```

## See Also

- [`.howard` catalog format](../formats/howard.md) — the format this function reads (schema, boundary law, trust model)
- [Writing filesets (`COPY … TO`)](write_pfile.md#sharded-output-shard_by) — the sharded writer that emits the `.howard.json` (the write→read partner)
- [`read_pfile`](read_pfile.md) — read a surviving shard's genotypes back
