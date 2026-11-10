# Writing filesets (`COPY … TO`)

Write PLINK 2 filesets back out of SQL. PlinkingDuck registers four DuckDB **COPY
formats** — one bundle and three standalone writers:

| `FORMAT` | Writes | Input |
|----------|--------|-------|
| `pfile` | a whole fileset: `prefix.pvar` + `prefix.pgen` + `prefix.psam` | one row per variant; genotypes as a per-sample `LIST`/`ARRAY` or `STRUCT` |
| `pvar` | a text `.pvar` | variant-metadata columns |
| `pgen` | a binary `.pgen` (biallelic hardcalls) | one row per variant; genotypes column |
| `psam` | a minimal `#IID`-only `.psam` | sample IIDs (out-of-band or per row) |

These are **writers**, invoked through `COPY (query) TO 'target' (FORMAT …)`, not
table functions — so they do not appear in the [function reference
table](index.md). The headline round-trips the reader:

```sql
COPY (SELECT * FROM read_pfile('in', genotypes := 'struct'))
TO 'out' (FORMAT pfile);
```

!!! warning "COPY option syntax is `KEY value`, not `:=`"

    DuckDB **COPY** options are space-separated `KEY value` pairs
    (`FORMAT pfile, samples ['A','B']`). The `:=` form used by the *reader* table
    functions (`read_pfile('x', genotypes := 'struct')`) is the table-function
    argument operator — the COPY parser rejects it with `Parser Error: syntax error
    at or near ":="`. Struct keys are barewords; struct, scalar, and list option
    values are string/number/list literals.

## Scope (what this release writes)

This is the **v0.10.0 Phases A–C** writer. It writes **biallelic diploid
hardcalls** from a **variant-oriented, position-sorted** stream, to **local**
paths. See [Limitations](#limitations-and-deferrals) for everything that is
deliberately out of scope — read it before relying on the writer; the deferrals
are real, not polish.

Later releases add **sharding** on top of that base: [`shard_by`](#sharded-output-shard_by)
along the variant axis, [`shard_samples`](#sample-sharded-output-shard_samples) along
the sample axis, and [both at once](#matrix-output-shard_by-and-shard_samples) as a matrix.
Each emits a [`.pfile.json` catalog](../formats/pfile_catalog.md).

## `FORMAT pfile` — the bundle

```sql
COPY (SELECT … ORDER BY CHROM, POS) TO 'prefix' (FORMAT pfile
  [, variant_fields {chrom: 'CHROM', pos: 'POS', id: 'ID', ref: 'REF', alt: 'ALT'}]
  [, genotype_fields 'genotypes']
  [, samples ['IID1', 'IID2', …]]
  [, variant_count N]
  [, pvar 'path.pvar'] [, pgen 'path.pgen'] [, psam 'path.psam']);
```

One `COPY` fans out to all three companion files. A trailing **`.pfile`** on the
target is stripped to the prefix, so `TO 'out'` and `TO 'out.pfile'` both write
`out.pvar` / `out.pgen` / `out.psam`.

### Options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `variant_fields` | `STRUCT` | each role ← same-named column | Map the `.pvar` roles (`chrom`,`pos`,`id`,`ref`,`alt`) to input columns. Keys are barewords, values are the input column names as strings. Any role left out defaults to its same-named column (case-insensitive). |
| `genotype_fields` | `VARCHAR` | `'genotypes'` | Name of the input genotypes column. |
| `samples` | `LIST(VARCHAR)` | *(from STRUCT field names)* | Explicit sample IIDs when the genotypes column is a plain `LIST`/`ARRAY` (see [Sample identity](#sample-identity-contract)). |
| `variant_count` | `INTEGER` | `16777216` (2²⁴) | Upper bound on the variant count (a pgenlib memory knob; see [below](#the-variant_count-bound)). |
| `pvar` / `pgen` / `psam` | `VARCHAR` | `prefix.{pvar,pgen,psam}` | Per-file path overrides (symmetric with the reader's path overrides). |
| `shard_by` | `VARCHAR` | *(none)* | Shard along the **variant** axis by a key column. See [Sharded output](#sharded-output-shard_by). |
| `shard_samples` | `STRUCT` | *(none)* | Shard along the **sample** axis into named IID groups. See [Sample-sharded output](#sample-sharded-output-shard_samples). |

Every input column not consumed by a `variant_fields` role (and not the genotypes
column) is carried through into the `.pvar` as an extra column under its own name,
in input order — mirroring `read_pvar`'s "extra columns" output, so a
`read_pfile` → `COPY` → `read_pfile` round-trip preserves QUAL/FILTER/INFO/etc.

### Sample identity contract

A real `.psam` is written **only from explicit sample identity** — the writer
never fabricates IIDs:

1. **`genotypes := 'struct'`** — the STRUCT field *names* are the IIDs (field
   order = `.pgen` sample order = `.psam` IID order). This is the natural source
   from `read_pfile(…, genotypes := 'struct')`.
2. **`samples ['id1', …]`** — an explicit IID list, for a plain `LIST`/`ARRAY`
   genotypes column that carries no identity of its own.

If **neither** is present, `COPY … (FORMAT pfile)` raises a clean bind error and
writes **no files** — a plain `TINYINT[N]` genotypes stream has no sample identity
to record. If **both** are present they must agree (same IIDs, same order); a
`samples` list that disagrees with the STRUCT field names is an error, so an
explicit `samples` is never silently ignored. When the genotypes column fixes the
sample count (a fixed-size `ARRAY`, or a `STRUCT`), a `samples` list of a different
length is a bind-time error.

!!! note "Per-sample metadata is IID-only in this release"

    The `.psam` carries **only `#IID`**. `sample_fields` for anything beyond `iid`
    (`sex`, phenotype, covariates) raises an explicit *"not yet supported"* error —
    requested metadata is never silently dropped. A variant-oriented stream has no
    per-sample axis to carry that metadata; writing it is a follow-on.

### Local-only

Because `.pgen` is written through a raw `FILE*` that bypasses DuckDB's VFS, the
**whole fileset must be local**. A remote or scheme target (`s3://`, `http://`,
`pathmacro://`) is rejected at bind — including a remote *prefix* with a local
`pgen` override (which would otherwise write `.pvar`/`.psam` to the remote store).
To land a fileset on a remote store, write it locally and copy it up separately.

### Atomic overwrite

The three files are written to `.staging` siblings in the same directory and
renamed onto their final paths only when all three are complete. So a failure at
any point **before** the final promote leaves a pre-existing fileset
**byte-identical** — a failed overwrite never clobbers the prior data, and no
`.staging` residue is left behind.

The promote itself (three renames) is **not atomic as a group**: a failure between
the renames can leave a new `.pvar`/`.psam` beside the old `.pgen`. In almost all
such cases the new and old filesets have different variant counts, which
`read_pfile` catches as a variant-count mismatch — a clean error, never a silent
mis-serve. The one uncaught edge is a `.pgen`-rename failure on a *same-variant-
count* overwrite, which would serve the new metadata over the old genotypes
undetected; that requires an intra-directory `rename(2)` to fail on a just-synced
file, which is effectively impossible on a local filesystem. (This is strictly
safer than truncating the real targets up front.)

### Input must be position-sorted

The `.pvar` writer enforces, in-stream, that variants are **position-sorted and
contiguous per chromosome** (`POS` non-decreasing within a chromosome; a
chromosome never reappears after another). Unsorted input raises an error pointing
you at `ORDER BY CHROM, POS`. This matches what PLINK tooling and the reader's
region pushdown assume, and a round-trip test alone cannot catch disorder.

### Sharded output (`shard_by`)

```sql
COPY (SELECT … ORDER BY CHROM, POS) TO 'out' (FORMAT pfile, shard_by 'CHROM');
```

`shard_by '<col>'` writes **one pfile per contiguous value of `<col>`** — the
common case is sharding a cohort by chromosome — plus **one shared `.psam`** and a
**[`.pfile.json` catalog](../formats/pfile_catalog.md)** that
[`read_pfile_catalog`](read_pfile_catalog.md) and the Python
[howard](https://github.com/teaguesterling/howard) loader read back as a single
concatenated dataset. The layout is:

```
out/<col>=<value>/data.pgen      # one fileset per shard (hive-style dir)
out/<col>=<value>/data.pvar
out.psam                         # ONE shared .psam (same samples across all shards)
out.pfile.json                  # the catalog — written LAST, the completion marker
```

The single ordered `COPY` rolls shards in the sink: on a `<col>` transition it
finalizes and promotes the open shard's `.pgen`/`.pvar`, then opens the next. This
is `PARTITION_BY` **semantics** implemented by hand — native DuckDB `PARTITION_BY`
is unusable for `.pgen` (it forces an unordered parallel sink that interleaves
rows, can split one key across files via `max_open_files` eviction, and strips the
partition column).

- **Input must be sorted so each key value is contiguous** — `ORDER BY <col>, POS`.
  A key value that reappears after another is rejected (the same contiguity
  contract the `.pvar` monotonicity enforces, lifted to the shard key).
- **Variant-shards over a shared `.psam`.** All shards share the one root `.psam`
  (same samples, different variants), so the catalog's `concat_dim` is `variants`.
  Reading a shard directly back with `read_pfile` needs the shared `.psam`:
  `read_pfile('out/CHROM=1/data', psam := 'out.psam')`.
- **`variant_count N` bounds *each shard*** (it re-arms per shard — per-chromosome,
  not genome-wide), so the default is almost always enough even for whole-genome
  input.
- **The catalog is the completion marker.** It is emitted only after every shard
  and the shared `.psam` are durable. A mid-run failure leaves no `.pfile.json`
  and no `.staging` residue; shards that already closed before the failure remain
  on disk, but without the catalog the set reads as incomplete-by-construction.
- **Overwriting an existing target.** A prior run's `.pfile.json` is removed just
  before the first shard is promoted, because from then on it would describe a mix
  of new and old shards. A re-run that fails after that point leaves no catalog
  (the cohort reads as absent), never a stale one. A re-run that fails earlier
  leaves the prior cohort intact.
- **`shard_by` cannot be combined with `pvar`/`pgen`/`psam` path overrides** (the
  sharded layout generates each shard's paths), and the key must be a scalar
  variant column (not the genotypes column). The whole tree is **local-only**, as
  with the non-sharded bundle.

The emitted `.pfile.json` records, per shard, the global half-open `[start, stop)`
variant-ordinal `partition` (tiling `[0, total)`), the shard's own `pvar`, and a
`range` of the covered `CHROM`(s) + `POS` min/max (a pruning stat the writer fills
for free as it streams). The shared `.psam` is referenced once under `shared.psam`.

!!! warning "Re-running to an existing target can orphan stale shard directories"

    A shard-key value must be **path-safe** — it becomes a real directory name
    `out/<col>=<value>`. Values with a path separator, `=`, whitespace, a control
    character, a glob metacharacter (`* ? [ ] { }`, because the catalog reader treats
    a shard prefix as a path pattern), or the names `.`/`..`/empty are rejected with
    an actionable error (CHROM names, including contigs such as `GL000192.1`, and
    ordinary alphanumeric/`_`/`-` labels are fine).

    Re-running a sharded write to an **existing** target that produces a
    **different** set of keys leaves the prior run's orphaned shard directories on
    disk — the writer deliberately does **not** delete directories on the target
    (auto-deletion on a target is exactly where data-loss bugs live). The
    `.pfile.json` catalog governs what is loaded (the loader never globs), so this
    is **not** a correctness issue for catalog-driven reads. But a consumer that
    globs `out/<col>=*/…` *outside* the catalog would see the stale filesets, and a
    stale shard could pair with a changed shared `.psam`. If you rely on globbing
    rather than the catalog, write to a fresh target (or clear it first).

### Sample-sharded output (`shard_samples`)

```sql
COPY (SELECT * FROM read_pfile('cohort', genotypes := 'struct') ORDER BY CHROM, POS)
TO 'out' (FORMAT pfile, shard_samples {batch1: ['s1', 's2'], batch2: ['s3', 's4', 's5']});
```

`shard_samples` splits the **sample** axis. Each named group becomes one fileset
with its own `.pgen` and `.psam`. All the filesets share **one `.pvar`**, because
the variants are the same in every shard. The layout is:

```
out.pvar                    # ONE shared .pvar
out/<label>/data.pgen       # per shard
out/<label>/data.psam       # per shard: that group's IIDs
out.pfile.json             # the catalog (concat_dim "samples") — written LAST
```

[`read_pfile`](read_pfile.md) reads the shards back as one cohort with
`combine_samples := 'concatenate'`. `ATTACH 'out.pfile.json' (TYPE pfile_catalog)`
attaches the cohort. A shard has no `.pvar` of its own, so a direct read of one
shard names the shared `.pvar`: `read_pfile('out/batch1/data', pvar := 'out.pvar')`.

- **The groups must partition the input samples exactly.** Every sample must be
  in exactly one group. An IID that is not an input sample is an error.
- **Each group is a contiguous run of the input sample order, and the groups follow
  the input order.** A gap, a reordered group, or groups listed out of order is an
  error. This keeps the catalog's sample ordinals correct: concatenating the shards
  in catalog order gives the input order back.
- **A label becomes a directory name** and part of every catalog prefix, so it may
  contain only ASCII letters, digits, `_` and `-`. The catalog reader treats a
  prefix as a path pattern, so a glob character such as `[` or `?` is rejected. A
  label that is a SQL keyword must be quoted (`{"all": […]}`).
- **Sample identity** is resolved as for the plain bundle: STRUCT field names, or
  `samples [...]`. A bare `TINYINT[N]` list is rejected.
- **Variant-orient input only.** No transpose: a sample-orient stream fails on the
  missing `CHROM` role.
- **Combinable with `shard_by`** to shard both axes at once — see
  [Matrix output](#matrix-output-shard_by-and-shard_samples). **Not combinable** with
  `pvar`/`pgen`/`psam` path overrides (the layout generates every path).
- **Memory and file handles.** Unlike `shard_by`, which keeps one writer open at a
  time, every variant row goes to every shard, so all N `.pgen` writers are open
  at once. Each one reserves about 1.5 bytes per variant of `variant_count` (about
  25 MB at the default 2²⁴) and holds 2 file handles. So a COPY needs N times that
  memory and 2N handles. To reduce it, pass a smaller `variant_count` (an upper
  bound on the variants written), or use fewer groups. The COPY fails at bind,
  before any file is created, when the estimate exceeds the host's memory or the
  open-file limit (`ulimit -n`).
- **Failure leaves nothing.** Files open on the first row. Any failure, including an
  upstream query error mid-stream, removes every staging file and every directory
  the COPY created (only if empty). On success the writer promotes the `.pvar`,
  then every `.psam`, then every `.pgen`, then writes the catalog last. When it
  overwrites an existing target, it removes the prior catalog just before the
  first promotion. A failure during promotion then leaves no catalog, never a
  stale one, and a failure before promotion leaves the prior cohort intact.

The emitted `.pfile.json` records `shared.pvar`, and per shard its `prefix`, its
`psam`, and a global half-open `[start, stop)` **sample**-ordinal `partition` that
tiles `[0, total)`.

### Matrix output (`shard_by` and `shard_samples`)

Give **both** `shard_by` and `shard_samples` in one `COPY` to shard the two axes at
once. The result is a **matrix**: an M × N grid of `.pgen` cells — M variant-blocks by
N sample-blocks — described by a
[version-2 `pfile_catalog`](../formats/pfile_catalog.md#version-2-the-matrix-form).

```sql
COPY (SELECT * FROM read_pfile('cohort', genotypes := 'struct') ORDER BY CHROM, POS)
TO 'out' (FORMAT pfile, shard_by 'CHROM',
          shard_samples {batch1: ['s1', 's2'], batch2: ['s3', 's4', 's5']});
```

The layout is one `.pvar` per variant-block, one `.psam` per sample-block, and one
`.pgen` per cell:

```
out/CHROM=<v>/variants.pvar             # M, one per variant-block
out/GROUP=<label>/samples.psam          # N, one per sample-block
out/CHROM=<v>/GROUP=<label>/data.pgen   # M x N cells
out.pfile.json                          # the v2 catalog — written LAST
```

The grid is the `shard_by` **variant roll** with the `shard_samples` **sample
fan-out** inside each variant-block. The roller holds one variant-block open at a
time; inside it, N `.pgen` writers each take their sample slice. So **N** writers are
open at once, not M × N — the sample count, not the grid size, sets the footprint.
The N `.psam` are the same for every variant-block, so each is written once.

**The emitted catalog** (`pfile_catalog: 2`) keys the variant axis with **two**
compound keys and the sample axis with one:

- `CHROM` (an `enum`). It carries `"over": "CHROM"` **only** when the `shard_by`
  column is `CHROM` (case-insensitive), which lets `region` prune cells by contig. A
  non-`CHROM` shard key omits `over`, and `region` is then refused.
- `CHUNK` (an `ordinal_interval` over `POS`), with the per-block `POS` `{min, max}`
  bounds the writer fills as it streams. This lets `region` prune cells by position,
  not only by contig.
- `GROUP` (an `enum`) on the sample axis — one label per sample-block.

The grid is **dense**: every variant-block pairs with every sample-block, so the
catalog has no `holes`. Read it back the usual ways:

- [`read_pfile_catalog`](read_pfile_catalog.md) — one row per cell, with prune and
  verify. `region := '2:…'` prunes on both the contig and the `POS` bounds.
- `ATTACH 'out.pfile.json' (TYPE pfile_catalog)` — the structural relations
  (`c.variants`, `c.samples`, `c.cells`, `c.variant_blocks`, `c.sample_blocks`) and
  the variant-oriented `c.genotypes`. Each cell is read through its template-resolved
  `.pvar` and `.psam`, so a cell needs no companion of its own.
- One cell, directly, with the block's shared companions:
  `read_pfile('out/CHROM=1/GROUP=batch1/data', pvar := 'out/CHROM=1/variants.pvar', psam := 'out/GROUP=batch1/samples.psam')`.

**Reserved key names.** The writer names the variant `POS` key `CHUNK` and the sample
key `GROUP`. A `shard_by` column of either name (case-insensitive) would put one key
name on both axes, so it is rejected at bind. The shard key and the target basename
each become a catalog template literal, so a name with a `{` or `}` brace is rejected
too.

**Guards and rejects** — all at bind, before any file is written:

- **The resource guard counts N**, the concurrent writer count — the same guard as
  [`shard_samples`](#sample-sharded-output-shard_samples). Too many groups for the
  host's memory or open-file limit fails the `COPY`, naming the mitigations (a smaller
  `variant_count`, or fewer groups).
- The same identity, orientation, contiguity, and label rules as `shard_samples`
  apply (the sample-axis code is shared): a bare `TINYINT[N]` list with no identity, a
  sample-oriented stream, a non-contiguous group, or a bad label is rejected.
- `pvar`/`pgen`/`psam` path overrides cannot be combined (the layout generates every
  path). Zero input rows write nothing.

**Failure leaves nothing.** Every file stages until the whole grid is complete, and
the catalog is promoted last. Any failure — upstream, in the sink, or a shard key that
reappears after another — removes every staging file and every directory the `COPY`
created. An overwrite replaces the whole grid and leaves the catalog valid; a mid-run
failure leaves no catalog, never a stale one.

!!! warning "The sample-oriented read of a matrix is not supported yet"

    `c.genotypes_by_sample` (the sample-oriented ATTACH relation) does **not** work on
    a writer-produced matrix. A multi-sample-block grid is guarded, and a
    multi-variant-block grid cannot supply a per-block `.pvar` to a sample-oriented
    read. Both need the two-dimensional read (P13-006). The variant-oriented
    `c.genotypes`, the structural relations, and `read_pfile_catalog` are correct. See
    [What version 2 does not do yet](../formats/pfile_catalog.md#what-version-2-does-not-do-yet).

## Standalone formats

### `FORMAT pvar`

```sql
COPY (SELECT … ORDER BY CHROM, POS) TO 'out.pvar' (FORMAT pvar
  [, variant_fields {…}] [, pvar 'override.pvar']);
```

Writes a tab-delimited text `.pvar` re-readable by [`read_pvar`](read_pvar.md).
Output headers are always canonical (`#CHROM POS ID REF ALT`) even when
`variant_fields` remaps a role to a differently-named input column. `NULL` cells
are written as `.` (and read back as `NULL`). A cell containing a tab, newline, or
carriage return is **rejected** (it cannot be represented in a tab-delimited
`.pvar` without corrupting the layout). The same `(CHROM, POS)` monotonicity check
as the bundle applies.

### `FORMAT pgen`

```sql
COPY (SELECT … genotypes … ORDER BY CHROM, POS) TO 'out.pgen' (FORMAT pgen
  [, genotype_fields 'genotypes'] [, variant_count N] [, pgen 'override.pgen']);
```

Writes a PLINK 2 `.pgen` of **biallelic hardcalls** from a variant-oriented
stream, re-readable by [`read_pgen`](read_pgen.md). The genotypes column is a
`LIST`/`ARRAY` of per-sample calls (`0`/`1`/`2`, and `-9` or SQL `NULL` =
missing), or a `STRUCT` whose field values (in field order) are the per-sample
calls. The sample count is the genotype list length and must be identical on every
row. `read_pgen` needs the `.pvar`/`.psam` companions to read the result back, so
standalone `.pgen` writing is for when you already have (or separately write)
those companions — otherwise use `FORMAT pfile`.

#### The `variant_count` bound

`pgenlib` pre-allocates an in-memory index proportional to a variant-count **upper
bound** given up front. The default (2²⁴ ≈ 16.7M variants) covers exome / chip /
single-chromosome scale at ~25 MB of reserved address space. For larger datasets,
pass `variant_count N` with any value `≥` the number of variants being written.
Exceeding the bound is caught as a **clean error** (never heap corruption), so the
default is safe to exceed — it is a memory knob, not a hard cap. A value above
~2×10⁹ (a likely typo that would reserve multiple GB) is rejected at bind.

### `FORMAT psam`

```sql
COPY (query) TO 'out.psam' (FORMAT psam
  [, samples ['id1', …]] [, sample_fields {iid: '<col>'}] [, psam 'override.psam']);
```

Writes a minimal `#IID`-only `.psam` re-readable by [`read_psam`](read_psam.md).
The IID source is resolved in **priority order**:

1. a **STRUCT genotypes column** (its field names are the IIDs) — written once,
   out of band, independent of the row stream;
2. a **`samples ['id', …]`** list — written once, out of band;
3. otherwise a **per-row IID column** (default `IID`, remappable with
   `sample_fields {iid: '<col>'}`) — one `.psam` line per input row, for a
   sample-oriented input such as `SELECT IID FROM read_psam('x')`.

Duplicate IIDs are rejected (`read_psam` cannot read a dup-IID file). As with the
bundle, `sample_fields` roles beyond `iid` raise a *"not yet supported"* error.

## Examples

```sql
-- Headline: round-trip a fileset (struct field names become the .psam IIDs).
COPY (SELECT * FROM read_pfile('data/cohort', genotypes := 'struct'))
TO 'data/cohort_copy' (FORMAT pfile);
```

```sql
-- Write a region subset as a new fileset (reader pushdown keeps it cheap).
COPY (SELECT * FROM read_pfile('data/cohort', genotypes := 'struct',
                               region := '17:43000000-43125000'))
TO 'data/brca1' (FORMAT pfile);
```

```sql
-- Plain LIST genotypes carry no IIDs — supply them with `samples`.
COPY (SELECT CHROM, POS, ID, REF, ALT, genotypes::TINYINT[] AS genotypes
      FROM read_pfile('data/cohort'))
TO 'data/relabelled' (FORMAT pfile, samples ['A', 'B', 'C', 'D']);
```

```sql
-- Just the variant metadata, to a text .pvar.
COPY (SELECT * FROM read_pvar('data/cohort.pvar'))
TO 'data/variants.pvar' (FORMAT pvar);
```

```sql
-- Remap differently-named variant columns back to canonical headers.
COPY (SELECT CHROM AS CONTIG, POS AS BP, ID AS SNP, REF AS A2, ALT AS A1
      FROM read_pvar('data/cohort.pvar'))
TO 'data/out.pvar' (FORMAT pvar,
  variant_fields {chrom: 'CONTIG', pos: 'BP', id: 'SNP', ref: 'A2', alt: 'A1'});
```

## Limitations and deferrals

The writer is deliberately scoped. These are honest cuts, not oversights:

- **Per-sample metadata beyond IID** (SEX, phenotype, covariates) is **not
  written** — the `.psam` is `#IID`-only, and `sample_fields` beyond `iid` errors.
  A variant-oriented stream has no per-sample axis to carry it.
- **Phased genotypes and dosages are rejected.** A phased (`TINYINT[2][N]`) or
  dosage (`DOUBLE[N]`) genotypes column raises an actionable error at bind. Only
  biallelic diploid hardcalls are written (Phase B.2 follow-on).
- **Multiallelic sites collapse to biallelic in the `.pgen`.** The genotype
  *values* written are faithful to what `read_pfile` produced (`PgrGet` already
  collapses a multiallelic site to biallelic 0/1/2/missing), but allele *identity*
  is lost: you can get a 2-allele `.pgen` beside a `.pvar` whose ALT still reads
  `A,C`. There is a guard — a comma in **ALT** is rejected — but it only fires when
  a column *literally named* `ALT` is in the stream. The default round-trip carries
  ALT and is correctly rejected; but if you remap the ALT role off a differently-
  named column (`variant_fields {alt: 'A1'}`), or write standalone `FORMAT pgen`
  from a genotypes-only projection with no ALT column, the guard cannot see the
  comma and the site is written collapsed. Multiallelic writing is deferred.
- **No `.pvar` `##` meta lines.** Header/`##contig`/`##INFO` comment lines are not
  reconstructed, so a round-trip is **value-equal, not byte-equal**.
- **Position-sorted, chromosome-contiguous input is required** (the monotonicity
  check), and the whole fileset is **local-only**.
- **A matrix has no sample-oriented read yet.** `shard_by` and `shard_samples`
  combine into a [matrix](#matrix-output-shard_by-and-shard_samples), but its
  `c.genotypes_by_sample` ATTACH relation is not supported — the variant-oriented
  `c.genotypes` and the other relations are. It needs the 2-D read (P13-006).
- **Parquet companions, `.pgi` split index, and PLINK-VCF export are deferred**
  (Phases D–E).

## See Also

- [`pfile_catalog` catalog format](../formats/pfile_catalog.md) — the catalog the sharded writer emits (schema, boundary law, trust model), including the [version-2 matrix](../formats/pfile_catalog.md#version-2-the-matrix-form)
- [read_pfile_catalog](read_pfile_catalog.md) — read the emitted `.pfile.json` back in SQL
- [read_pfile](read_pfile.md) — read a fileset back (the round-trip partner)
- [read_pvar](read_pvar.md) · [read_pgen](read_pgen.md) · [read_psam](read_psam.md)
- [File Handling](../guides/file-handling.md#writing-filesets) — where writing fits
  in the fileset lifecycle
- [Genotype Encoding](../genotype-encoding.md) — the hardcall / missing encoding
