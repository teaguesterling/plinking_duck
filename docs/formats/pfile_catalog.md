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

## JSON schema

A catalog is either the **multi-file** form (`concat_dim` + `shared` + `files`) or
the **single-fileset** form (`pfile`), never both.

| Field | Where | Meaning |
|-------|-------|---------|
| `pfile_catalog` | top | Format version integer. Must be `1`. |
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

## See also

- [`read_pfile_catalog`](../functions/read_pfile_catalog.md) — read, verify, and prune a catalog in SQL
- [Writing filesets (`COPY … TO`)](../functions/write_pfile.md#sharded-output-shard_by) — emit a catalog with the sharded writer
- [`read_pfile`](../functions/read_pfile.md) — read a single fileset (how a surviving shard is scanned)
- [File Formats](../file-formats.md) — the underlying `.pgen` / `.pvar` / `.psam` fileset
