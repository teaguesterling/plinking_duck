# The `pfile_catalog` catalog format

`pfile_catalog` is a small, human-editable **JSON catalog** layered on top of **existing**
PLINK 2 filesets. It stores no genotype data of its own — it references `.pgen` /
`.pvar` / `.psam` files that already exist and records just enough metadata to
assemble them into one logical dataset: which files make up a cohort, how they
concatenate, and (for pruning) the genomic extent of each shard.

PlinkingDuck both **writes** `pfile_catalog` (the sharded [`COPY … (FORMAT pfile, shard_by …)`](../functions/write_pfile.md#sharded-output-shard_by)
writer emits one as its output manifest) and **reads** it back
([`read_pfile_catalog`](../functions/read_pfile_catalog.md) parses, verifies,
and prunes it) — closing the write→read loop. The Python
[howard](https://github.com/teaguesterling/howard) library also reads `pfile_catalog` to
register a cohort as a diar dataset; the format is its design, and this page
documents the subset PlinkingDuck implements.

!!! note "Why a catalog, not a new store"

    VCF-Zarr, Hail MatrixTable, and TileDB all **re-store** the genotype data in a
    new on-disk form. `pfile_catalog` deliberately does not: the canonical artifacts stay
    the foreign `.pgen`/`.pvar`/`.psam` files, and the catalog is a thin
    serialization pointing at them. JSON was chosen because DuckDB reads it
    natively (`read_json`), so the catalog is itself a queryable relation — the
    overlap test that prunes shards is a plain SQL query over it.

## What it is for

A single `read_pfile('cohort', pvar := …, psam := …)` already handles one fileset
with companion overrides. `pfile_catalog` earns its keep on the **multi-file shard case**
that SQL has no syntax for: several `.pgen` that **share** one `.pvar` or one
`.psam`, concatenated along an axis, each annotated with the global ordinal range
it occupies and the genomic range it covers. That annotation is what lets a query
skip shards it cannot touch without opening their `.pgen` at all (see
[Partition and range](#partition-and-range)).

The single-fileset case is in scope too, as declarative capture — it is barely more
than a written-down `read_pfile(…)`, but it gives `pfile_catalog` one shape for both.

## One extension, two shapes

A pfile catalog always uses the **`.pfile.json`** extension, whether it indexes a
single fileset or a sharded cohort. The shape is discriminated by the catalog's own
top-level keys, not by the file name — the reader routes on the presence of a
`pfile` key vs `files`/`concat_dim` (and rejects a catalog that declares both):

| Shape | How it reads | Use |
|-------|--------------|-----|
| **multi-file** | has `concat_dim` + `shared` + `files` | A **multi-file** cohort (one dataset assembled from several `.pgen`), e.g. `cohort.pfile.json`. This is what the sharded writer emits. |
| **single-fileset** | has `pfile` | A **single-fileset** sidecar sitting beside one fileset, e.g. `cohort.pfile.json` next to `cohort.{pgen,pvar,psam}`. |

A catalog describes **exactly one dataset**. Grouping several datasets into one
catalog file is out of scope.

Both shapes above are **version 1** and are 1-D: shards concatenate along one axis.
A grid of shards (variants × samples) is **version 2**. See
[Version 2: the matrix form](#version-2-the-matrix-form).

## JSON schema

A catalog is either the **multi-file** form (`concat_dim` + `shared` + `files`) or
the **single-fileset** form (`pfile`), never both.

| Field | Where | Meaning |
|-------|-------|---------|
| `pfile_catalog` | top | Format version integer. `1` for the forms on this page; `2` for the [matrix form](#version-2-the-matrix-form). Any other value is rejected. |
| `dataset` | top | The dataset name. |
| `concat_dim` | multi-file | `"variants"` (shards share one `.psam`) or `"samples"` (shards share one `.pvar`) — the axis the shards concatenate on. |
| `shared` | multi-file | The reused companion: `{ "psam": … }` for variant-shards, `{ "pvar": … }` for sample-shards. |
| `files` | multi-file | Ordered array of shards; each `{ prefix, <own-companion>, partition, range? }`. |
| `prefix` | per file | Fileset path with **no** extension (`<prefix>.pgen`/`.pvar`/`.psam`). |
| `pvar` / `psam` | per file | The shard's **own** non-shared companion (the shared one comes from `shared`). May be a **text** (`.pvar`/`.psam`) or **parquet** (`.pvar.parquet`/`.psam.parquet`) file. Optional — an omitted companion falls back to the reader's own discovery (see [Path resolution](#path-resolution)). |
| `partition` | per file | `{ start, stop, source_ordinal }` — the shard's **global half-open** ordinal span on the concat axis (see [below](#partition-and-range)). Required. |
| `range` | per file | `{ "CHROM": […], "POS": { "min", "max" } }` — the genomic extent, the pruning stat. **Strongly recommended** on every shard (it is the point of a shard catalog) but optional. |
| `pfile` | single-file | `{ prefix, pvar?, psam? }` for one fileset. |
| `chunks` | single-file | `{ dim, chunk_size }`; `dim` must be `"variants"`. Tells the loader how to chunk the single fileset. |

**Companion values** (`pvar`, `psam`, `shared.psam`, `shared.pvar`) are normally a
**string path**. A companion may instead be an **object**
`{ "path": …, "path_is_relative": …, "columns": { … } }` when a relative-path
opt-out or a column name-map is needed.

### Path resolution

Every `prefix` and companion path resolves **relative to the catalog file's own
directory** by default — so a catalog moves together with its files: copy the
`.pfile.json` and the filesets as a set and the references still resolve. An
**absolute path is honored as-is.** On a companion **object**,
`"path_is_relative": false` opts that one path out of relative resolution.

!!! note "Written relative, returned absolute"

    PlinkingDuck's writer emits **relative** paths (a shard `prefix` of
    `cohort/CHROM=1/data`, a `shared.psam` of `cohort.psam`). The reader resolves
    them against the catalog directory and returns **absolute, normalized** paths,
    so a consumer can use them regardless of its own working directory. `path_is_relative`
    is honored on **companion objects** only; a `prefix` is always taken through
    relative resolution (an absolute `prefix` is still honored). This holds for a
    parquet companion too: a declared `foo.pvar.parquet`, or a discovered sidecar, is
    returned as its absolute, normalized path.

### Companion kinds and omitted-companion discovery

A companion may be a **text** PLINK file (`.pvar`/`.bim`, `.psam`/`.fam`) or a **parquet**
sidecar (`.pvar.parquet`/`.psam.parquet`) — declared explicitly as that path, or **omitted**
and discovered. An omitted companion is resolved **exactly the way `read_pfile` discovers it**:
the `.pvar.parquet` / `.psam.parquet` sidecar **first** (parquet-first, when parquet companions
are enabled), then the native text `.pvar`/`.bim` (or `.psam`/`.fam`). The manifest then **surfaces
the resolved path** in its `pvar`/`psam` output columns — so a consumer reads the *same* file the
manifest verified, instead of re-discovering it (and so a bare `.pvar.parquet` sidecar is not missed
by a text-only `read_pvar('<prefix>.pvar')`). Verification is parquet-aware: the count and extent
checks read the resolved companion natively. A companion that resolves to **nothing** (no declared
file, no discoverable default) raises `file not found`.

### Example — a variant-sharded cohort (`cohort.pfile.json`)

Three `.pgen`, 8 shared samples, one shared `.psam`, concatenated along `variants`,
one contig per shard. This is the **actual output** of
`COPY (…) TO 'cohort' (FORMAT pfile, shard_by 'CHROM')` over a 3000-variant,
3-contig fileset:

```json
{
  "pfile_catalog": 1,
  "dataset": "cohort",
  "concat_dim": "variants",
  "shared": { "psam": "cohort.psam" },
  "files": [
    {
      "prefix": "cohort/CHROM=1/data",
      "pvar": "cohort/CHROM=1/data.pvar",
      "partition": { "start": 0, "stop": 1000, "source_ordinal": 0 },
      "range": { "CHROM": ["1"], "POS": { "min": 100, "max": 100000 } }
    },
    {
      "prefix": "cohort/CHROM=2/data",
      "pvar": "cohort/CHROM=2/data.pvar",
      "partition": { "start": 1000, "stop": 2000, "source_ordinal": 1 },
      "range": { "CHROM": ["2"], "POS": { "min": 100, "max": 100000 } }
    },
    {
      "prefix": "cohort/CHROM=3/data",
      "pvar": "cohort/CHROM=3/data.pvar",
      "partition": { "start": 2000, "stop": 3000, "source_ordinal": 2 },
      "range": { "CHROM": ["3"], "POS": { "min": 100, "max": 100000 } }
    }
  ]
}
```

The `partition.start`/`stop` are global `variants`-axis ordinals and tile
`[0, 3000)` with no gap or overlap. There is **no per-file `psam`** — every shard
reuses the one `shared.psam`. The `samples` axis is identical across all shards
(guaranteed by the shared `.psam` and the header invariant below).

### Example — a sample-sharded cohort (mirror)

The symmetric case: same variants, different samples, one shared `.pvar`,
concatenated along `samples`. Each shard carries its own `.psam`:

```json
{
  "pfile_catalog": 1,
  "dataset": "cohort",
  "concat_dim": "samples",
  "shared": { "pvar": "data/cohort.pvar" },
  "files": [
    {
      "prefix": "data/cohort.batchA",
      "psam": "data/cohort.batchA.psam",
      "partition": { "start": 0, "stop": 5000, "source_ordinal": 0 }
    },
    {
      "prefix": "data/cohort.batchB",
      "psam": "data/cohort.batchB.psam",
      "partition": { "start": 5000, "stop": 12000, "source_ordinal": 1 }
    }
  ]
}
```

Here `partition` ranges are global `samples`-axis ordinals (the `sample_index`
extent of each shard). A genomic `range` prunes nothing on this axis — the variants
are identical across sample-shards — so sample-shards prune on `sample_index` via
`partition` instead, and typically carry no `range`.

!!! note "PlinkingDuck's writer emits variant-shards only"

    The sharded `COPY … shard_by` writer produces the **variant-sharded** form
    (a shared `.psam`, `concat_dim: "variants"`). The sample-sharded form is a valid
    catalog that `read_pfile_catalog` reads, but writing it is a separate
    follow-on.

### Example — a single fileset (`cohort.pfile.json`)

```json
{
  "pfile_catalog": 1,
  "dataset": "wes_single",
  "pfile": { "prefix": "large_example" },
  "chunks": { "dim": "variants", "chunk_size": 1000 }
}
```

A single fileset declares no `partition`/`range` — the reader probes the fileset's
`.pvar` to fill them (`partition` = `[0, variant_count)`, `source_ordinal` = 0, plus
the contig set and POS extent). `pvar`/`psam` may be given as overrides (text or parquet);
omitted, they fall back to `read_pfile`'s own discovery (parquet sidecar first, then text), and the
resolved path is surfaced in the output — see [Companion kinds and omitted-companion discovery](#companion-kinds-and-omitted-companion-discovery).

## The shared-companion invariant

A `.pgen` header carries `raw_variant_ct` and `raw_sample_ct`, and PLINK 2 requires
the companions to match them exactly. So sharing **one** companion across several
`.pgen` is legal in exactly two orientations, and `concat_dim` names which:

- **Shared `.psam`** ⇒ every `.pgen` has the same `sample_ct` ⇒ the files are
  **variant-shards** (same samples, different variants), each with its own `.pvar`
  ⇒ `concat_dim: "variants"`.
- **Shared `.pvar`** ⇒ every `.pgen` has the same `variant_ct` ⇒ the files are
  **sample-shards** (same variants, different samples), each with its own `.psam`
  ⇒ `concat_dim: "samples"`.

There is no third orientation. A `shared.psam` with `concat_dim: "samples"` (or
vice versa) is a contradiction and is **rejected at parse** — as is a per-file
companion of the *shared* type.

## Partition and range

These two per-shard fields serve different jobs. Getting their boundary conventions
right is what makes pruning correct, so they are stated explicitly rather than left
to convention:

!!! abstract "Boundary law"

    - **`partition { start, stop }` is half-open `[start, stop)`** — `start`
      included, `stop` excluded. A shard holds ordinals `start … stop-1`.
    - **`range.POS { min, max }` is inclusive `[min, max]`** — both endpoints are
      real coordinates the shard contains.
    - **`range.CHROM` is a membership list** — a contig is covered or it is not.

**`partition`** is the shard's **global ordinal span** on the concat axis — *not* a
file-local range. Three variant-shards of 1000 variants each are `[0, 1000)`,
`[1000, 2000)`, `[2000, 3000)`, tiling `[0, total)` with no gap or overlap. This is
the required correctness field (it places each shard on the axis) and it doubles as
the **index prune** key: a predicate on `variant_index` / `sample_index` reads
`partition` directly, needing no `range`. The index is monotonic by construction, so
this prune works on **both** axes.

**`range`** is the **genomic extent** for region pruning — the one thing ordinals
cannot express ("which contig does this shard hold"). It applies to variant-shards
only. A region predicate keeps a shard iff its contig is in the shard's `CHROM` list
**and** the query's POS window intersects the shard's inclusive `POS` extent.

`range` is **strongly recommended** on every shard — it is the entire point of a
shard catalog — but formally optional so a minimal hand-written catalog still loads.
A shard **without** `range` is correct but un-prunable: it is **always** a candidate
and is **never** wrongly dropped (a conservative keep).

!!! warning "Only `CHROM` and `POS` are read from `range`"

    The format's design envelope allows other `range` keys (e.g. `CPRA`, or the
    monotonic index keys). PlinkingDuck's reader parses **only `CHROM` and `POS`**;
    any other key in a `range` object is silently ignored. The writer likewise emits
    only `CHROM` (always) and `POS` (when positions are present).

## Trust model

`pfile_catalog` bounds are **assertions about foreign files** PlinkingDuck did not write —
a `.pgen`/`.pvar` can be swapped, truncated, or re-sorted out from under the catalog.
A stale bound would then cause a shard to be **wrongly skipped** (a correctness bug,
not just a slow query). So, unlike Hail's constructive partition invariants, these
bounds are **verified, never assumed**:

- **Count check** — each shard's `.pgen` header `raw_variant_ct` / `raw_sample_ct`
  is asserted against its companion row counts, and the declared `partition` width
  (`stop - start`) must equal the real companion count. A mismatch raises at open.
- **Range check** (variants only) — a declared `range` must **cover** the shard's
  real `.pvar` extent (declared ⊇ actual). An **under-claim** raises (it would
  wrongly drop the shard); an over-claim is merely less effective pruning.

By default only the **survivors** of a prune are verified. `verify_all := true`
checks every shard before pruning, closing the hole where an under-claiming range on
a shard that was *itself* dropped would never be opened. See
[`read_pfile_catalog`](../functions/read_pfile_catalog.md#verification-the-trust-model)
for the full behavior.

## Producing a `pfile_catalog`

The sharded pfile writer emits one as its output manifest:

```sql
COPY (SELECT * FROM read_pfile('cohort', genotypes := 'struct') ORDER BY CHROM, POS)
TO 'cohort' (FORMAT pfile, shard_by 'CHROM');
```

This writes one fileset per contiguous `CHROM` value under a hive layout
(`cohort/CHROM=<value>/data.{pgen,pvar}`), **one** shared `cohort.psam`, and
`cohort.pfile.json` **last** — the catalog is the completion marker, so a
mid-run failure leaves no catalog and the shard set reads as
incomplete-by-construction. The writer fills each shard's `range` for free as it
streams the data. See [Writing filesets](../functions/write_pfile.md#sharded-output-shard_by)
for the full layout, contract, and limitations.

## Consuming a `pfile_catalog`

In SQL, [`read_pfile_catalog`](../functions/read_pfile_catalog.md) parses,
verifies, and prunes the catalog, returning one row per surviving shard:

```sql
-- Every shard, with resolved paths, partition tiling, and range extent.
SELECT * FROM read_pfile_catalog('cohort.pfile.json');

-- Genomic prune: only the shards whose range overlaps chr2:1-5Mb (and verify them).
SELECT prefix, pvar FROM read_pfile_catalog('cohort.pfile.json', region := '2:1-5000000');
```

This is the **manifest** surface (Phase A): it resolves *which* shards a query needs
and vouches for them. It does not itself scan genotypes or `ATTACH` — a consumer
reads each surviving shard with [`read_pfile`](../functions/read_pfile.md), using the
resolved `prefix` and the shared companion the manifest reports. The Python
`pfile_catalog` library consumes the same catalog to register the cohort as a single
concatenated diar dataset.

## What PlinkingDuck validates — and what it leaves to the loader

PlinkingDuck's reader is a *manifest* over the catalog, not the full dataset loader.
It deliberately checks less than the end-to-end howard loader does, so a catalog that
passes `read_pfile_catalog` is not guaranteed to be a complete, gapless cohort:

- **Per-shard checks only.** The reader validates each shard's `partition` width
  (`stop > start`), rejects a duplicate `partition.start`, and sorts shards by
  `start`. It does **not** check the cross-shard `[0, total)` tiling (no-gap /
  no-overlap contiguity) — that stays in the howard loader, which owns registration.
- **No index rebasing.** When a prune drops shards, `read_pfile_catalog` returns
  the survivors with their **original global** `partition_start`/`stop` and
  `source_ordinal` — it does not rebase the subset to `[0, subset_total)`. (The
  howard Python loader does rebase, and warns when it does; the SQL manifest does
  not.) If you read the manifest after a prune, the ordinals are still global.
- **`columns`, `record_count`, `file_size_bytes` are not used.** A `columns`
  name-map on a companion object is accepted and ignored; `record_count` /
  `file_size_bytes` staleness stamps are not parsed at all.
- **Companions may be text or parquet, and an omitted one is resolved — not left to the consumer.**
  A `.pvar`/`.psam` may be declared as text or as a `.pvar.parquet`/`.psam.parquet` sidecar, and an
  omitted companion is resolved by `read_pfile`'s discovery order (**parquet-first**, then native).
  The manifest both **verifies** the resolved companion (parquet-aware) and **surfaces its resolved
  path** in the output, so the consumer reads the exact file the manifest checked. What is still left
  to the loader is the companion's *content* beyond the row count and (for `.pvar`) the CHROM/POS
  extent — PlinkingDuck does not validate every field.

## Version 2: the matrix form

Version 1 tiles one axis. Version 2 tiles **two** axes at once: a grid of `.pgen` files
where variants are split one way and samples are split the other way.

### Version discriminator

A catalog is version 2 when its top-level `"pfile_catalog"` is the integer `2`. It must
then have an `axes` block, a `template` object, and a `dataset` name. Version 1 is
unchanged. Any other version number is rejected. A version-2 catalog must **not** also
carry `concat_dim`, `files`, `shared`, or `pfile`. The reader rejects the mix.

### Data model

- A matrix has two **axes**, named `variants` and `samples`. Both are required. The
  format has two axes today.
- Each axis is cut into ordered **blocks**. The blocks tile the axis `[0, length)` with
  no gap and no overlap.
- A **cell** is one variant-block × one sample-block. It is one `.pgen` fileset.
- Companions are shared **along** each axis. All cells in a variant-block share that
  block's `.pvar`. All cells in a sample-block share that block's `.psam`. A grid of
  M variant-blocks and N sample-blocks has M `.pvar`, N `.psam`, and M×N `.pgen`.
- Block spans are half-open `[start, stop)`, as in the [boundary law](#partition-and-range).

### Grammar

Each axis has a **key schema** (`keys`) and a **tiling** (`blocks`). The catalog
defines its own key names. Only the key **types** are fixed.

| Key type | Extra field | Per-block value | Bound the block carries |
|----------|-------------|-----------------|-------------------------|
| `enum` | none | A label string. | None. |
| `ordinal_interval` | `over`: `"POS"` on `variants`, `"sample_index"` on `samples` | A label string. | Over `POS`: `bounds.<key>` = `{ "min", "max" }`, **inclusive**. Over `sample_index`: none. The block's own `[start, stop)` is the bound. |
| `value_bin` | `over`: the name of the binned value (for example `"age"`) | A label string. | `bounds.<key>` = `{ "lo", "hi" }`, **half-open** `[lo, hi)`. |

Bound values are integers.

| Field | Where | Meaning |
|-------|-------|---------|
| `axes[].name` | axis | `"variants"` or `"samples"`. Each appears once. |
| `axes[].keys[]` | axis | Ordered key schema: `{ name, type, over? }`. At least one key. A key name is unique across **both** axes. |
| `axes[].blocks[]` | axis | The tiling, in ordinal order: `{ keys, start, stop, bounds? }`. |
| `blocks[].keys` | block | An object with a **string** value for **every** key the axis declares. The same combination of values may not appear twice on an axis. |
| `blocks[].start`, `stop` | block | Integer ordinals, half-open. The first block starts at 0. Each next block starts where the last one stopped. `stop` must be greater than `start`. |
| `blocks[].bounds` | block | Required for each `value_bin` key and each `ordinal_interval` key over `POS`. `hi` (or `max`) may not be less than `lo` (or `min`). |
| `template` | top | `{ "pfile", "pvar", "psam" }`. Path templates (see below). |
| `holes` | top | Optional. Absent cells (see [Sparse grids](#sparse-grids-holes)). |

`ordinal_interval` over the wrong quantity for its axis (for example `POS` on `samples`)
is a parse error. So is an unknown key type.

!!! note "Short form for one `enum` key"

    An axis with exactly one `enum` key may omit `blocks`. It lists
    `"values": [{ "value", "start", "stop" }, …]` inside the key instead. The reader
    turns each entry into one block. An axis with more than one key needs `blocks`.
    An axis may not use both forms.

### Path templates

`template.pfile`, `template.pvar`, and `template.psam` are paths with `{KEY}`
placeholders. For each cell, the reader replaces `{KEY}` with that cell's value for the
key. Paths resolve relative to the catalog directory, as in
[Path resolution](#path-resolution). `template.pfile` has no extension.

- `template.pfile` must use at least one key of **each** axis. Otherwise two cells
  resolve to the same `.pgen`.
- `template.pvar` may use **only** `variants` keys, and at least one.
- `template.psam` may use **only** `samples` keys, and at least one.
- A placeholder that names an undeclared key, or an unclosed `{`, is a parse error.
- Two cells that resolve to the same `.pgen`, or two blocks that resolve to the same
  `.pvar` or `.psam`, are a parse error.

Omitted-companion discovery (see [above](#companion-kinds-and-omitted-companion-discovery))
does not apply to version 2. The templates always name the companions.

### Example — a 3 × 4 grid with compound keys

Variants are keyed by `CHROM` (an `enum`) and `CHUNK` (an `ordinal_interval` over
`POS`). Samples are keyed by `GROUP` (an `enum`), `SLICE` (an `ordinal_interval` over
`sample_index`), and `age` (a `value_bin`). This is the catalog the
`read_pfile_catalog_matrix_compound` test builds. It has 3 variant-blocks and 4
sample-blocks, so 12 cells, 3 `.pvar`, and 4 `.psam`.

```json
{
  "pfile_catalog": 2,
  "dataset": "compound",
  "axes": [
    {
      "name": "variants",
      "keys": [
        { "name": "CHROM", "type": "enum" },
        { "name": "CHUNK", "type": "ordinal_interval", "over": "POS" }
      ],
      "blocks": [
        { "keys": { "CHROM": "1", "CHUNK": "001" }, "start": 0, "stop": 2,
          "bounds": { "CHUNK": { "min": 10000, "max": 20000 } } },
        { "keys": { "CHROM": "1", "CHUNK": "002" }, "start": 2, "stop": 3,
          "bounds": { "CHUNK": { "min": 30000, "max": 30000 } } },
        { "keys": { "CHROM": "2", "CHUNK": "001" }, "start": 3, "stop": 4,
          "bounds": { "CHUNK": { "min": 15000, "max": 15000 } } }
      ]
    },
    {
      "name": "samples",
      "keys": [
        { "name": "GROUP", "type": "enum" },
        { "name": "SLICE", "type": "ordinal_interval", "over": "sample_index" },
        { "name": "age", "type": "value_bin", "over": "age" }
      ],
      "blocks": [
        { "keys": { "GROUP": "A", "SLICE": "0", "age": "young" }, "start": 0, "stop": 1,
          "bounds": { "age": { "lo": 20, "hi": 40 } } },
        { "keys": { "GROUP": "A", "SLICE": "1", "age": "old" }, "start": 1, "stop": 2,
          "bounds": { "age": { "lo": 40, "hi": 60 } } },
        { "keys": { "GROUP": "B", "SLICE": "0", "age": "young" }, "start": 2, "stop": 3,
          "bounds": { "age": { "lo": 20, "hi": 40 } } },
        { "keys": { "GROUP": "B", "SLICE": "1", "age": "old" }, "start": 3, "stop": 4,
          "bounds": { "age": { "lo": 40, "hi": 60 } } }
      ]
    }
  ],
  "template": {
    "pfile": "c_{CHROM}_{CHUNK}_{GROUP}_{SLICE}",
    "pvar": "v_{CHROM}_{CHUNK}.pvar",
    "psam": "s_{GROUP}_{SLICE}.psam"
  }
}
```

`template.pvar` uses `CHROM` and `CHUNK`, so each variant-block gets its own `.pvar`.
`template.psam` uses `GROUP` and `SLICE`, so each sample-block gets its own `.psam`.
The `age` key labels the sample-blocks and drives the `age` prune. It is not in a
template, which is allowed.

### Example — build a catalog from existing exports (`glob` + JSON)

You often have the filesets on disk already, but no catalog. You can build one in
SQL. Fix the **axes** by hand — you know your data has a `CHROM` axis and a `GROUP`
axis. Discover the **blocks** (each value and its ordinal span) from the files with
`glob`, so you never hand-type the chromosomes or the groups.

This example assumes the shared-companion layout, relative to the catalog directory:

```
cohort/
  variants/CHROM=1.pvar, CHROM=2.pvar, …      one shared .pvar per variant-block
  samples/GROUP=A.psam,  GROUP=B.psam,  …      one shared .psam per sample-block
  cells/CHROM=1/GROUP=A/cell.pgen, …           one .pgen per cell
```

One query assembles the catalog. The two axes and their key names are literals. The blocks
of each axis come from the files on disk, and a running count turns each block's size into
its contiguous `[start, stop)` ordinal span. The array order you emit **is** the axis order
(the format tiles the blocks as listed), so each axis is ordered on purpose: the variant
axis by an **explicit karyotypic contig order**, the sample axis by file order (groups have
no canonical order).

```sql
COPY (
  WITH
  -- The contig order is EXPLICIT — one row per contig, in the order you want them tiled
  -- (karyotypic here). Use your own contig names (a 'chr' prefix is fine, be consistent). A
  -- contig absent from this table fails loud below, so none is silently dropped or misordered.
  contig_order(value, rank) AS (VALUES
    ('1',1),('2',2),('3',3),('4',4),('5',5),('6',6),('7',7),('8',8),('9',9),('10',10),
    ('11',11),('12',12),('13',13),('14',14),('15',15),('16',16),('17',17),('18',18),
    ('19',19),('20',20),('21',21),('22',22),('X',23),('Y',24),('XY',25),('MT',26)),
  -- Variant blocks: read_pvar gives each contig and its count; contig_order fixes the order.
  -- The LEFT JOIN + coalesce(…, error(…)) raises on any contig missing from contig_order.
  variant_counts AS (
    SELECT CHROM AS value, count(*) AS n
    FROM read_pvar('cohort/variants/*.pvar') GROUP BY CHROM),
  variant_ranked AS (
    SELECT c.value, c.n,
           coalesce(o.rank,
             error('pfile_catalog: contig ' || c.value || ' is missing from contig_order')) AS rank
    FROM variant_counts c LEFT JOIN contig_order o USING (value)),
  variant_vals AS (
    SELECT value,
           row_number() OVER w - 1      AS idx,    -- block position along the axis
           (sum(n) OVER w)::BIGINT - n   AS start,  -- running ordinal: where this block begins
           (sum(n) OVER w)::BIGINT       AS stop
    FROM variant_ranked
    WINDOW w AS (ORDER BY rank ROWS UNBOUNDED PRECEDING)),
  -- Sample blocks: GROUP is in the companion path; glob fixes the order, read_csv counts.
  sample_files AS (
    SELECT file, regexp_extract(file, 'GROUP=([^/.]+)', 1) AS value
    FROM glob('cohort/samples/*.psam')),
  sample_counts AS (
    SELECT regexp_extract(filename, 'GROUP=([^/.]+)', 1) AS value, count(*) AS n
    FROM read_csv('cohort/samples/*.psam', filename := true, delim := '\t',
                  header := true, all_varchar := true)
    GROUP BY value),
  sample_vals AS (
    SELECT f.value,
           row_number() OVER w - 1         AS idx,
           (sum(c.n) OVER w)::BIGINT - c.n  AS start,
           (sum(c.n) OVER w)::BIGINT        AS stop
    FROM sample_files f JOIN sample_counts c USING (value)
    WINDOW w AS (ORDER BY f.file ROWS UNBOUNDED PRECEDING)),
  v AS (SELECT list({'value': value, 'start': start, 'stop': stop} ORDER BY idx) AS arr FROM variant_vals),
  s AS (SELECT list({'value': value, 'start': start, 'stop': stop} ORDER BY idx) AS arr FROM sample_vals)
  SELECT
    2        AS pfile_catalog,
    'cohort' AS dataset,
    [ {'name': 'variants', 'keys': [{'name': 'CHROM', 'type': 'enum', 'values': v.arr}]},
      {'name': 'samples',  'keys': [{'name': 'GROUP', 'type': 'enum', 'values': s.arr}]} ] AS axes,
    {'pfile': 'cells/CHROM={CHROM}/GROUP={GROUP}/cell',
     'pvar':  'variants/CHROM={CHROM}.pvar',
     'psam':  'samples/GROUP={GROUP}.psam'} AS template
  FROM v, s
) TO 'cohort/cohort.pfile.json' (FORMAT json);
```

`COPY … (FORMAT json)` writes each result row as one JSON object, keyed by the column
names — so this single row is the whole catalog object, with no quoting to undo. `idx` is
each block's position on its axis and the running `sum` is its ordinal; `sum` returns a
128-bit total, so the spans are cast to `BIGINT`, which the parser reads. The result is a
normal catalog:

```sql
SELECT keys['CHROM'] AS chrom, keys['GROUP'] AS grp,
       variant_stop - variant_start AS n_variants,
       sample_stop - sample_start  AS n_samples, verified
FROM read_pfile_catalog('cohort/cohort.pfile.json');
```

Three notes:

- **Contig order is the author's job, and lexicographic order is wrong for it.** A variant
  axis keyed by `CHROM` must tile in karyotypic order (`1, 2, …, 22, X, Y, XY, MT`). File or
  lexicographic order puts `CHROM=10` before `CHROM=2`, so the variant ordinal stops tracking
  genomic position. The example fixes the order with an explicit `contig_order` table — one row
  per contig, written in your own contig names and in the order you want. A contig absent from
  the table fails loud (the `coalesce(…, error(…))` guard), so none is silently dropped or
  sorted to an undefined position. The sample axis has no canonical order, so file order is
  fine there. `verified` checks each block's width against its own `.pvar` and that the blocks
  tile with no gap or overlap; it does **not** check that the order is genomic. That is on you.
- This is the single-`enum`-key form, where each block is one `value`. An
  `ordinal_interval` key (a `CHUNK` over `POS`) also needs per-block `bounds`. Compute
  them with `min(POS)` and `max(POS)` per block and emit the `blocks` form instead (see
  the [compound example](#example--a-3--4-grid-with-compound-keys)).
- If your exports keep a companion in every cell rather than one per axis, point the
  template at any one cell's companion per block. The `.pvar` is the same across a
  variant-block's cells, and the `.psam` is the same across a sample-block's cells.

### Sparse grids (`holes`)

A grid need not be full. The optional top-level `holes` array lists **absent** cells.
Each entry is an object of `key: "value"` pairs over declared keys of either axis. A
hole removes every cell whose block labels match **all** of its pairs. A partial hole
removes a whole row or column of the grid:

```json
"holes": [ { "CHROM": "2", "GROUP": "C" } ]
```

The example above (a 2 × 3 grid with `CHROM` ∈ {1, 2} and `GROUP` ∈ {A, B, C}) has 5
cells, not 6. `{ "GROUP": "C" }` alone would remove both `GROUP=C` cells.

A hole is never opened, and its path is never resolved. A hole may therefore point at
files that do not exist. The reader rejects:

- a hole that matches no cell (a typo must not leave the cell present);
- a hole that names an undeclared key;
- an empty hole `{}`;
- a hole value that is not a string;
- a `holes` value that is not an array.

### Reading and pruning a matrix

[`read_pfile_catalog`](../functions/read_pfile_catalog.md) returns **one row per
surviving cell**. Its output columns, prune parameters, and verification are described
on the function page. In short:

- `keys` prunes by key, on `enum` keys only (membership). An interval or bin key is
  rejected; use the positional parameters below.
- `region`, `variant_index`, and `samples` prune by position.
- A cell survives only if it passes **every** predicate. The predicates multiply
  across the two axes.
- A predicate that names an undeclared key, or the wrong axis, raises an error. It is
  never ignored.

### Inspecting a matrix without reading it

Two functions look at the catalog file alone and open no data file:

- [`read_pfile_catalog_metadata`](../functions/read_pfile_catalog_metadata.md) returns the
  **literal declared schema**: one row per declared key, with its axis, `type`, `over` and
  position.
- [`read_pfile_catalog_layout`](../functions/read_pfile_catalog_layout.md) returns the
  **interpreted chunk map**: one row per grid cell, **holes included** (`present = false`),
  with the 2-D spans and the template-resolved (unchecked) paths.

### Trust model for a matrix

Version 2 keeps the version 1 rule: bounds are **verified, never assumed**. For each
verified cell:

- The `.pgen` header counts must match the cell's `.pvar` and `.psam` row counts.
- The variant span width (`variant_stop - variant_start`) must equal the real `.pvar`
  row count. The sample span width must equal the real `.psam` row count.
- If the variants axis has an `enum` key named `CHROM`, each contig the `.pvar` holds
  must equal the block's `CHROM` label.
- If the variants axis has an `ordinal_interval` over `POS`, the block's `{ min, max }`
  must cover the real `.pvar` POS extent. An under-claim raises. An over-claim only
  prunes less.

By default only the cells that survive a prune are verified. `verify_all := true`
verifies every present cell first. Holes are never verified. See
[Verification](../functions/read_pfile_catalog.md#matrix-catalogs-version-2).

### What version 2 does not do yet

- **No `ATTACH`.** A matrix catalog cannot be attached as a database.
- **No genotype cascade.** `read_pfile_catalog` returns cell metadata and resolved
  paths. It does not scan genotypes across cells. A consumer reads each surviving cell
  with [`read_pfile`](../functions/read_pfile.md).
- **No analysis over the grid.** The `plink_*` analysis functions do not take a matrix
  catalog.
- **No writer.** The sharded `COPY … shard_by` writer emits version 1 only. A matrix
  catalog is written by hand or by other tooling.
- **Two axes only.** The format names `variants` and `samples`. Other axes are not
  supported.
- **`region` needs a `CHROM` key.** `region` matches its contig against a variants-axis
  `enum` key named exactly `CHROM`. A catalog with no such key rejects `region`. All
  other key names are free.

## See also

- [`read_pfile_catalog`](../functions/read_pfile_catalog.md) — read, verify, and prune a catalog in SQL
- [`read_pfile_catalog_metadata`](../functions/read_pfile_catalog_metadata.md) — the catalog's literal declared schema
- [`read_pfile_catalog_layout`](../functions/read_pfile_catalog_layout.md) — the catalog's interpreted chunk map (holes included, no file I/O)
- [Writing filesets (`COPY … TO`)](../functions/write_pfile.md#sharded-output-shard_by) — emit a catalog with the sharded writer
- [`read_pfile`](../functions/read_pfile.md) — read a single fileset (how a surviving shard is scanned)
- [File Formats](../file-formats.md) — the underlying `.pgen` / `.pvar` / `.psam` fileset
