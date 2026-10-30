# read_pfile_catalog

Parse, verify, and prune a `.pfile.json` shard catalog, returning **one row per surviving shard**.
For a version-2 matrix catalog, it returns **one row per surviving cell** (see
[Matrix catalogs](#matrix-catalogs-version-2)).

This is the parse + verify + prune surface plinking_duck owns for `pfile_catalog` **reading**. The **write**
side is the sharded pfile writer (`COPY (…) TO 'out' (FORMAT pfile, shard_by '<col>')`, see
[Writing filesets](write_pfile.md)), which emits exactly the `out.pfile.json` this function reads —
closing the write→read loop. The format is howard's `pfile_catalog` catalog design (variant-shards over a
shared `.psam`, or sample-shards over a shared `.pvar`); this function is Phase A of consuming it and
does **not** itself scan genotypes or `ATTACH` — it resolves which shards a query needs and vouches
for them.

For the catalog format itself — the JSON schema, the boundary law, and the trust model — see the
[**`pfile_catalog` catalog format reference**](../formats/pfile_catalog.md).

## Synopsis

```sql
read_pfile_catalog(path VARCHAR
    [, region := VARCHAR,
       variant_index := ANY,
       sample_index := ANY,
       keys := MAP(VARCHAR, VARCHAR[]),
       samples := ANY,
       verify_all := BOOLEAN]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `path` | `VARCHAR` | *(required)* | Path to a multi-file `.pfile.json` catalog. Every path inside it resolves **relative to the catalog's own directory** (an absolute path is honored as-is). |
| `region` | `VARCHAR` | *(none)* | Genomic prune, `'contig:lo-hi'` or a bare `'contig'`. Keeps only shards whose declared `range` could overlap. Valid **only** when `concat_dim = 'variants'`. POS bounds are **inclusive** `[lo, hi]`. |
| `variant_index` | `ANY` (2-element list `[lo, hi]`) | *(none)* | Index prune on the variants axis. Keeps shards whose `partition` overlaps the **inclusive** index predicate `[lo, hi]`. Valid **only** when `concat_dim = 'variants'`. |
| `sample_index` | `ANY` (2-element list `[lo, hi]`) | *(none)* | Index prune on the samples axis, inclusive `[lo, hi]`. Valid **only** when `concat_dim = 'samples'`. |
| `keys` | `MAP(VARCHAR, VARCHAR[])` | *(none)* | **Version 2 only.** Prune by key. See [Matrix catalogs](#matrix-catalogs-version-2). A version-1 catalog rejects it. |
| `samples` | `ANY` (list of `BIGINT`) | *(none)* | **Version 2 only.** Keep cells whose sample-block holds one of the listed sample indices. Pass 0-based integer indices; a string (IID) element raises `invalid predicate`, because IID lookup is not supported yet. A version-1 catalog rejects it. |
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
    subset to `[0, subset_total)`. (The Python `pfile_catalog` loader rebases a pruned sub-cohort to a
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

A companion may be **declared** (an explicit `pvar`/`psam` path — which can itself be a parquet file,
e.g. `foo.pvar.parquet`) or **omitted**. An **omitted** companion (on any catalog form) is **resolved
the way `read_pfile` discovers it** — the `.pvar.parquet` / `.psam.parquet` sidecar **first** (when
parquet companions are enabled), then the native text `.pvar` / `.bim` (or `.psam` / `.fam`) — and the
**resolved path is emitted** in the output column (it is **no longer `NULL`**), so a consumer reads the
exact file the manifest verified rather than re-discovering it. **Parquet companions are supported**
throughout: both as a declared `.pvar.parquet` / `.psam.parquet` path and via sidecar discovery, and
the manifest **verifies** them — the count check (through `read_pfile`) and the width/extent checks
(`read_pvar` for the `.pvar` extent, `read_psam` for the `.psam` count, on the resolved path) read a
parquet companion natively. A declared companion
that is missing, or an omitted companion whose default file is **truly absent**, still raises
`file not found` (see [Error taxonomy](#error-taxonomy)).

## Verification (the trust model)

`pfile_catalog` bounds are *assertions about foreign files plinking_duck did not write*, so they are checked
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
| `pvar` | `VARCHAR` | **Resolved** `.pvar` path — the exact file `read_pfile` would open: the shard's own for variant-shards, the shared one for sample-shards, or, when omitted, the discovered text `.pvar` **or `.pvar.parquet` sidecar**. Never `NULL` for a verified shard. |
| `psam` | `VARCHAR` | **Resolved** `.psam` path — the shared one for variant-shards, the shard's own for sample-shards, or, when omitted, the discovered text `.psam` **or `.psam.parquet` sidecar**. Never `NULL` for a verified shard. |
| `partition_start` | `BIGINT` | Global half-open `[start, stop)` ordinal on the concat axis — included. |
| `partition_stop` | `BIGINT` | Global half-open ordinal — excluded. |
| `source_ordinal` | `BIGINT` | The shard's descriptive source ordinal. |
| `chrom` | `LIST(VARCHAR)` | `range.CHROM` (the contigs the shard covers); `NULL` when the shard declares no `range`. |
| `pos_min` | `BIGINT` | Inclusive `range.POS.min`; `NULL` when absent. |
| `pos_max` | `BIGINT` | Inclusive `range.POS.max`; `NULL` when absent. |
| `verified` | `BOOLEAN` | Whether this shard's companions were checked during this call. Currently always `true` on every returned row (a verification failure raises); the column exists for a future degrade-gracefully mode. |
| `chunk_size` | `BIGINT` | A **single-fileset** catalog's declared `chunks.chunk_size` (consumed by howard's single-file loader). `NULL` for multi-shard rows — there the shards themselves are the chunks. |

## Matrix catalogs (version 2)

A catalog with `"pfile_catalog": 2` is a grid of cells: variant-blocks × sample-blocks.
The reader routes on that version number. The catalog format is described in
[Version 2: the matrix form](../formats/pfile_catalog.md#version-2-the-matrix-form).

The function returns **one row per surviving cell**. The row order is variants-major.
The output columns differ from version 1:

| Column | Type | Description |
|--------|------|-------------|
| `dataset` | `VARCHAR` | The catalog's `dataset` name. |
| `prefix` | `VARCHAR` | The cell's `.pgen` path, no extension. Absolute and normalized. |
| `pvar` | `VARCHAR` | The `.pvar` shared by every cell in this cell's variant-block. |
| `psam` | `VARCHAR` | The `.psam` shared by every cell in this cell's sample-block. |
| `variant_start`, `variant_stop` | `BIGINT` | The variant-block's half-open ordinal span `[start, stop)`. |
| `sample_start`, `sample_stop` | `BIGINT` | The sample-block's half-open ordinal span `[start, stop)`. |
| `source_ordinal` | `BIGINT` | Position of the cell among the **present** cells, variants-major, counted before any prune. A prune does not renumber it. |
| `keys` | `MAP(VARCHAR, VARCHAR)` | Every key of both axes and its value for this cell, for example `{CHROM=1, CHUNK=001, GROUP=A, SLICE=0, age=young}`. |
| `verified` | `BOOLEAN` | Always `true`. A verification failure raises. |

### Prune parameters

A cell survives only if it passes **every** predicate you give. The predicates multiply
across the two axes. A prune that keeps nothing returns 0 rows. It is not an error.

| Parameter | Axis | Meaning |
|-----------|------|---------|
| `keys := {'GROUP': ['A', 'B']}` | either | Enum membership only: the cell's value must be in the list. A key of type `ordinal_interval` or `value_bin` is an `invalid predicate` error; use `region`, `variant_index` or `samples`. |
| `region := '1:10000-20000'` | variants | The cell's `CHROM` must equal the contig. If the variants axis has an `ordinal_interval` over `POS`, the block's inclusive `[min, max]` must also overlap the window. A bare contig (`'1'`) keeps every POS. |
| `variant_index := [lo, hi]` | variants | Inclusive `[lo, hi]` against the variant-block's half-open `[start, stop)`. |
| `samples := [i, j, …]` | samples | A **list of points**, not a range. The cell stays if any listed index lies in the sample-block's `[start, stop)`. |

Interval rules, by key type:

- `ordinal_interval` over `POS`: the block bound is inclusive `[min, max]`. The query
  `[lo, hi]` is inclusive. They overlap if they share at least one value.
- `ordinal_interval` over `sample_index`: the block bound is its half-open
  `[start, stop)`.
- `value_bin`: the block bound is half-open `[lo, hi)`. A query of `[40, 40]` overlaps
  `[40, 60)` but not `[20, 40)`.

A predicate that cannot apply raises `invalid predicate` and never returns the whole grid:

- `keys` names a key the catalog does not declare, or a key with a `NULL` value or
  `NULL` list element.
- `keys` on an interval key does not hold exactly two integers, or `lo > hi`.
- `region` on a catalog with no variants-axis contig key: an `enum` key marked `"over":"CHROM"`, or else one named `CHROM`.
- `sample_index` on a version-2 catalog. Use `samples` instead.
- `keys` or `samples` on a version-1 catalog.

`variant_index` and `region` keep the meaning they have for version 1: inclusive
`[lo, hi]`, and an inclusive POS window.

### Holes

A cell that the catalog lists in `holes` is absent. The reader skips it **before** it
resolves any path. It never opens the cell, so a hole may point at files that do not
exist. A hole never appears in the output. A predicate that selects only holes returns
0 rows. It does not raise `file not found`. `verify_all := true` also skips holes.

### Verification

The trust model is the same as for version 1: bounds are checked, not assumed. A
failure raises at bind. For each verified cell the reader checks:

- the cell's `.pgen` header counts against its `.pvar` and `.psam` row counts;
- the variant span width against the real `.pvar` row count, and the sample span width
  against the real `.psam` row count (each distinct companion is read once);
- each contig in the `.pvar` against the block's contig label, when the variants axis
  has a contig `enum` key (marked `"over":"CHROM"`, or else named `CHROM`);
- the block's declared POS bound against the real `.pvar` POS extent, when the variants
  axis has an `ordinal_interval` over `POS`. A bound that does not cover the real extent
  is a `range under-claim`.

The `verified` column does **not** cover `value_bin` key bounds. A `value_bin` bound
(for example an age bin) has no stored column to check against, so `verify_all` cannot
catch a stale bin. Only the `CHROM` and `POS` bounds above are checked against data.

The errors use the same prefixes as in [Error taxonomy](#error-taxonomy).

Verification is **lazy on survivors**. By default only the cells that survive the prune
are verified. A bad bound on a cell that the prune drops is not seen. This is the same
gap as in version 1. `verify_all := true` verifies **every present cell before** the
prune. Use it for a hand-written catalog.

### Not supported yet

- `ATTACH` of a matrix catalog.
- Reading genotypes through the catalog (the genotype cascade). Read each surviving
  cell with [`read_pfile`](read_pfile.md), using its `prefix`, `pvar`, and `psam`.
- Analysis functions (`plink_*`) over the grid.
- A writer for version 2. The sharded writer emits version 1.

### Matrix examples

```sql
-- All present cells of the grid, with the shared companions of each.
SELECT keys, prefix, pvar, psam FROM read_pfile_catalog('compound.pfile.json');

-- Enum membership on the samples axis: GROUP A only.
SELECT count(*) FROM read_pfile_catalog('compound.pfile.json', keys := {'GROUP': ['A']});

-- Genomic prune: contig 1, POS 15000. Keeps CHROM=1/CHUNK=001 only (not CHROM=2).
SELECT DISTINCT keys['CHROM'], keys['CHUNK']
FROM read_pfile_catalog('compound.pfile.json', region := '1:15000-15000');

-- Sample points 1 and 3: the sample-blocks [1,2) and [3,4).
SELECT count(*) FROM read_pfile_catalog('compound.pfile.json', samples := [1, 3]);

-- Both axes at once: the product of a region and a sample-axis key.
SELECT count(*) FROM read_pfile_catalog('compound.pfile.json',
    region := '1:10000-20000', keys := {'GROUP': ['A']});
```

## Error taxonomy

Errors carry distinct, classifiable message prefixes so a caller can tell failure classes apart:

| Prefix | Class | Exception |
|--------|-------|-----------|
| `read_pfile_catalog: catalog parse error: …` | Malformed JSON or schema (bad version, missing field, both `pfile` and multi-file keys, duplicate `partition.start`, out-of-range integer; for version 2: a bad key schema, a tiling gap or overlap, a bad template, a bad hole) | `InvalidInputException` |
| `read_pfile_catalog: file not found: …` | A referenced catalog / `.pgen` / `.pvar` / `.psam` is missing | `IOException` |
| `read_pfile_catalog: verification failed (count mismatch): …` | A `.pgen` header count disagrees with its companion | `InvalidInputException` |
| `read_pfile_catalog: verification failed (companion unreadable): …` | A companion exists but could not be read (corrupt `.pgen` header, wrong-format `.pvar`/`.psam`) | `InvalidInputException` |
| `read_pfile_catalog: verification failed (range under-claim): …` | A declared `range` does not cover the real `.pvar` extent | `InvalidInputException` |
| `read_pfile_catalog: invalid predicate: …` | Wrong-axis or malformed prune argument | `InvalidInputException` |

## Examples

```sql
-- One row per shard, with each shard's resolved paths, partition tiling, and range extent.
SELECT * FROM read_pfile_catalog('cohort.pfile.json');

-- Genomic prune: only the shards whose range overlaps chr2:1-5Mb survive (and are verified).
SELECT prefix, pvar FROM read_pfile_catalog('cohort.pfile.json', region := '2:1-5000000');

-- Index prune on the variants axis (inclusive): shards touching variant ordinals 1000..2000.
SELECT prefix FROM read_pfile_catalog('cohort.pfile.json', variant_index := [1000, 2000]);

-- Strict mode for a hand-written catalog: verify every shard's range before pruning.
SELECT count(*) FROM read_pfile_catalog('cohort.pfile.json', region := '7:1-1000000', verify_all := true);
```

## See Also

- [`pfile_catalog` catalog format](../formats/pfile_catalog.md) — the format this function reads (schema, boundary law, trust model)
- [Writing filesets (`COPY … TO`)](write_pfile.md#sharded-output-shard_by) — the sharded writer that emits the `.pfile.json` (the write→read partner)
- [`read_pfile`](read_pfile.md) — read a surviving shard's genotypes back
