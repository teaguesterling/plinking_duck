# ATTACH genotype views: tall `c.genotypes`, labelled wide views, an ID index

Date: 2026-11-19. Status: design, for Teague's approval.
Decisions by Teague in the session (AskUserQuestion) are marked **(Teague)**. Calls made by
Claude are marked **(Claude's call)**.

## Intent

An `ATTACH`-ed cohort (`ATTACH 'x' AS c (TYPE pfile_catalog)`) must make these reads easy,
with a plain `WHERE`, and fast through filter pushdown:

1. The genotypes in a region, for all samples.
2. The genotypes in a region, for a subset of samples.
3. One variant (or a few) by `ID`, for example `WHERE ID = 'rs123'`.

**Teague's common queries (2026-11-19).** These are the acceptance tests:

| # | Ask | Query |
|---|---|---|
| Q1 | All genotypes for rsid123 | `FROM c.genotypes WHERE ID = 'rsid123'` |
| Q2 | All non-ref genotypes in chr1:3000-4000 | `FROM c.genotypes WHERE CHROM = '1' AND POS BETWEEN 3000 AND 4000 AND genotype > 0` |
| Q3 | The IIDs with a hom-alt call in a list of rsids | `SELECT DISTINCT IID FROM c.genotypes WHERE ID IN (…) AND genotype = 2` |

A result row must say which sample (and which variant) each genotype belongs to. Today
`c.genotypes` returns one row per variant with an unlabelled genotypes ARRAY/LIST. The user
must know that element *i* is row *i* of `c.samples`.

## What exists today (measured on c1b5545's code, local build of 949e382)

`src/pfile_catalog_storage.cpp:312-338` builds the views:

| View | SQL behind it |
|---|---|
| `c.genotypes` | `read_pfile(catalog := …, freq := true)`. A matrix with holes: a `UNION ALL` of per-cell reads (`MatrixGenotypesSql`), with `genotypes := 'list'` when sample blocks differ in width. A bare fileset or a v1 list: `GenotypesSql(…)` with `freq := true`. |
| `c.genotypes_by_sample` | A matrix: `read_pfile(catalog := …, orient := 'sample')`. Else `GenotypesSql(…, by_sample = true)`. No freq columns. |

The engine already has every read shape this design needs:

- `orient := 'genotype'` (tall: one row per variant × sample, with `IID` and `genotype`).
- `genotypes := 'columns'` (one scalar column per IID, or per variant in sample orient).
- `genotypes := 'struct'`.
- A `WHERE` on `CHROM`/`POS` drops variant blocks. A `WHERE IID = …` / `IID IN (…)` drops
  sample blocks in genotype orient (PRs #96-#109).

Measurements. Data: 400,000 variants × 200 samples, written as a 20 × 2 matrix (`$S/att2/mx`);
and 200 variants × 1,000,000 samples (`$S/m1/m`). One run each, `srun -c 8` on longbottom.

| Query | Time |
|---|---|
| `c.genotypes WHERE POS BETWEEN 199000 AND 199999` (matrix) | 0.008 s, one variant block opened |
| `c.genotypes WHERE ID = 'snp199998'` (matrix) | 0.115 s: every `.pvar` is parsed, no block is dropped |
| `read_pfile('mx', orient := 'genotype') WHERE POS BETWEEN … AND IID IN (2 IIDs)` | 0.009 s |
| `read_pfile('mx', orient := 'sample', genotypes := 'columns') WHERE IID IN (2)`, 1 column | 0.776 s |
| 1M samples, 1 variant: LIST | 0.17 s |
| 1M samples, 1 variant: tall, `IID = 's1'` | 0.28 s |
| 1M samples, 1 variant: `genotypes := 'columns'` | 0.9 s |
| 1M samples, 1 variant (or 10): `genotypes := 'struct'` | 4.3 s, a fixed cost per query; peak memory 2.7 GB |

The values agreed across the four shapes (`s1` = 2 in each).

**"No freq columns" was not reproduced.** On the local build, `c.genotypes` has the six freq
columns for a matrix, a v1 sharded cohort and a bare fileset. The freq change (f0a145a, #113) is
in no release: `git tag --contains f0a145a` is empty, and the latest release is v0.17.0. So the
published extension has no freq columns. This is inferred; Teague has not yet said which build he
uses.

## Goals

- G1. `c.genotypes` is tall: one row per variant × sample. **(Teague)**
- G2. `c.genotypes_by_variant` is one row per variant with **one column per sample, named by
  IID**. **(Teague: wide columns, after the STRUCT measurement.)**
- G3. `c.genotypes_by_sample` is one row per sample with **one column per variant, named by
  ID**. **(Teague)**
- G4. `WHERE ID = '…'` / `ID IN (…)` open only the variant blocks that hold those IDs, and
  `WHERE IID` uses the same kind of index, through optional **index files** beside the file set
  or the catalog (§4). **(Teague)**
- G4b. A `WHERE` on `genotype` in the tall view (`= 2`, `> 0`, `IN (1, 2)`, `<> 0`, `IS NULL`)
  becomes the `include_genotypes` filter inside the scan, not a filter on decoded rows. Measured
  at 1M samples: Q2-shape 2.17 s with `WHERE genotype > 0` vs 1.05 s with
  `include_genotypes := ['het','hom_alt']`; Q3-shape 0.80 s vs 0.45 s. The test data is 75%
  non-ref; on real, mostly hom-ref data the gain is expected to be larger (inferred, not measured).
- G4c. List versions of the two wide views: `c.genotypes_by_variant_list` and
  `c.genotypes_by_sample_list`. **(Teague)**
- G4d. `ID` stays exactly what the `.pvar` holds, in every row view. A new virtual column
  `CPRA` (`CHROM:POS:REF:ALT`, VARCHAR) is built only when a query projects it, and a
  `WHERE CPRA = '…'` / `IN (…)` is parsed into `CHROM`/`POS` bounds plus an allele check.
  **(Teague, 2026-11-19.)** See §2c.
- G5. All views keep today's block drops: `CHROM`/`POS` drop variant blocks, `IID` drops
  sample blocks (where the read shape allows it).

## Non-goals

- A new query function or table macro. The views plus `WHERE` are the interface.
- Grid pruning inside `plink_freq` / `plink_missing` (#131's open ruling).
- A 2-D read of a matrix with holes, in this work. See §1b for what each view does on holes.
- Multiallelic support.

## Anti-goals

- A view whose cost depends on the cohort size even when the `WHERE` names one region. The
  wide views break this at bind (see Risks). That cost is accepted for the wide views only.
- Results that differ between the views for the same cells. Each view must equal a
  `read_pfile` of the un-sharded source, filtered the same way (the oracle used since #92).

## The design

### 1. View names and shapes (breaking)

| View | Row | Columns | SQL |
|---|---|---|---|
| `c.genotypes` | variant × sample | `CHROM, POS, ID, REF, ALT, [extra .pvar cols], IID, [FID, SEX, extra .psam cols]?, genotype, ALT_FREQ, OBS_CT, HOM_REF_CT, HET_CT, HOM_ALT_CT, MISSING_CT` | `read_pfile(catalog := …, orient := 'genotype', freq := true)` |
| `c.genotypes_by_variant` | variant | `CHROM, POS, ID, REF, ALT, …, <IID 1>, …, <IID n>, ALT_FREQ, …` | `read_pfile(catalog := …, genotypes := 'columns', freq := true)` |
| `c.genotypes_by_sample` | sample | `IID, [FID, SEX, …], <ID 1>, …, <ID m>` | `read_pfile(catalog := …, orient := 'sample', genotypes := 'columns')` |
| `c.genotypes_by_variant_list` | variant | today's `c.genotypes`: `CHROM, POS, ID, REF, ALT, …, genotypes, ALT_FREQ, …` | `read_pfile(catalog := …, freq := true)` |
| `c.genotypes_by_sample_list` | sample | today's `c.genotypes_by_sample`: `IID, …, genotypes` | `read_pfile(catalog := …, orient := 'sample')` |

In the list views, element *i* of `genotypes` is row *i* of `c.samples` (by variant) or of
`c.variants` (by sample). The docs say so. **(Claude's call: no ordinal column is added to
`c.variants` or `c.samples` in this work.)** Neither `c.samples` (a `UNION ALL` of `.psam`
reads) nor `c.variants` (a list `read_pvar`) has an `ORDER BY` today, so SQL does not promise
that order. PR 3 adds a test that element *i* of a list row equals the genotype of row *i* of
`c.samples` / `c.variants`, at 1 and 4 threads. If the order can drift, PR 3 adds an `ORDER BY`
on the block key and the in-block ordinal. **(Claude's call.)**

Which `.psam` columns the tall view carries is whatever genotype orient emits today. **(Claude's
call: no new columns.)**

This is a breaking change to `c.genotypes`. The release notes and `docs/formats/pfile_catalog.md`
say: "the old `c.genotypes` (one row per variant, a genotypes array) is now
`c.genotypes_by_variant_list`; the old `c.genotypes_by_sample` is now
`c.genotypes_by_sample_list`." **(Claude's call: no compatibility alias. Teague may overrule.)**

The non-catalog fallbacks (a bare fileset, a v1 list read) get the same shapes, by adding the
same `orient` / `genotypes` arguments to `GenotypesSql`.

### 1b. A matrix with holes

**The target (Teague):** a hole's cells are NULL. `genotypes_by_variant` and
`genotypes_by_sample`, and their list versions, give one row per variant (or per sample), with
NULL in each position that falls in a hole.

**This work (Teague: "can raise unsupported for now"):** on a matrix with holes, the four wide and
list views raise a "not supported yet" error. The message names `c.genotypes` as the read that
works. The error must come at bind, so `DESCRIBE` and `count(*)` get it too. This replaces
today's per-block fan-out of `c.genotypes` on a holed matrix. That fan-out is not kept under the
`_list` name, because each of its rows covers one sample block, so "element *i* is row *i* of
`c.samples`" would be false. **(Claude's call, from Teague's rule.)**

**The tall `c.genotypes` on a holed matrix** stays a `UNION ALL` of per-cell tall reads
(`MatrixGenotypesSql` with `orient := 'genotype'`). Each row names its own `IID`, so the rows are
correct. A hole gives no rows, not NULL rows. **(Claude's call: NULL rows for holes would be the
tall form of Teague's target; it is deferred with the 2-D-with-holes read, and the docs say
that a hole's rows are absent.)**

### 2. Freq columns in genotype orient

`read_pfile` rejects `freq := true` with `orient := 'genotype'` today ("its columns hold one value
per variant"). Change: allow it, and repeat the variant's six values on each of its sample rows.
The values are over the cohort's full sample set, not over the rows a `WHERE IID` keeps.
**(Claude's call.)** This keeps "rare variants in these samples" one query:
`SELECT * FROM c.genotypes WHERE ALT_FREQ < 0.01 AND IID IN (…)`.

A sample-block drop (`WHERE IID …`) must not change the freq values. So with freq projected, the
counts still come from every sample block. The builder decides how: read counts from all blocks,
or turn off the sample-block drop when a freq column is projected. The second option is simpler
and makes the query slower; the builder must say which one it picked, and measure it.

Sample orient keeps no freq columns. **(Claude's call: a per-sample row has no per-variant value.)**

### 2b. `WHERE genotype` pushdown (tall view)

In genotype orient, a pushed-down filter on the `genotype` column whose allowed set is a subset
of {0, 1, 2, missing} compiles to `GenotypeRangeFilter` (the type behind `include_genotypes`),
so a row the filter rejects is never emitted. The pushed filter still runs on the output
(today's residual-filter rule), so a wrong compile cannot add rows. A filter that cannot be
compiled stays a plain residual filter. A query that also passes `include_genotypes :=`
intersects the two. **(Claude's call.)**

Freq columns (§2) and `af_range`/`ac_range` must keep counting every call, as they do today with
`include_genotypes` (CLAUDE.md: counts report TRUE counts).

### 2c. The virtual `CPRA` column

- **Where:** `read_pfile` variant and genotype orient, so `c.genotypes`, `c.genotypes_by_variant`
  and `c.genotypes_by_variant_list` all have it, after `ALT` and before the extra `.pvar` columns
  **(Claude's call on position)**. It is projection-aware, like the freq columns: a query that
  does not name it pays nothing.
- **Format:** `CHROM:POS:REF:ALT`, the format that `ResolveCpraString`
  (`src/plink_common.cpp:2045`) already parses for `variants := [...]`. A multiallelic ALT
  stays one string (`1:3000:G:A,C`), because the parser splits only on `:`.
- **Pushdown:** a pushed `CPRA = const` or `CPRA IN (consts)` becomes `CHROM` + `POS` bounds,
  which reuse the existing variant-block drop and the #120 `.pgen` narrowing, plus a REF/ALT
  match. A `CHROM:POS` constant (two parts) matches the site with any alleles **(Claude's call,
  matching `ResolveCpraString`)**. An unparsable constant gives zero rows, not an error, because
  it is a `WHERE` value and not a parameter. The pushed filter still runs on the output rows, so
  a wrong parse can only lose rows. The tests must cover that case.
- **A `.pvar` that has its own `CPRA` column (Teague):** the view uses that column as is, in the
  same position. Each row that is read must have a `CPRA` equal to `CHROM:POS:REF:ALT`; a row
  that does not is an error that names the file and the row. "Parses as a CPRA" is not enough: a
  `chr1:…` value over a `CHROM` of `1`, or another allele normalization, would make the `POS`
  pushdown drop the wrong rows. Only the rows read are checked; a block that the pushdown drops
  is never checked. **(Claude's definition of "invalid".)**
- **No index:** `CPRA` lookups need no ID index; the `POS` bounds do the block drop. The ID index
  (§4) holds `ID` only.
- **Timing to measure in its PR:** `WHERE CPRA = …` on the 20 × 2 matrix against
  `WHERE POS = …` (8 ms measured). The two are expected to match (inferred).

### 3. Wide-column names and their limits

- `genotypes_by_variant`: column name = IID. This is `genotypes := 'columns'` today. A name that
  equals a variant column (`CHROM`, `POS`, `ID`, …, a freq column) must be an error, or renamed
  under a stated rule. The builder checks what `'columns'` does today and keeps it if it is
  defined.
- `genotypes_by_sample`: column name = variant ID. Today `read_pfile(orient := 'sample',
  genotypes := 'columns')` raises an error on a duplicate ID (`src/pfile_reader.cpp:2774`) and
  names an empty ID `CHROM:POS` (`:2741`). Real `.pvar` files often have many `.` IDs, so the
  view would raise that error on most cohorts. **(Teague: a view-only naming option, with the
  `CPRA` value as the fallback.)** A new `read_pfile` option, which only the view passes, gives
  an ID of `.`, an empty ID, or a duplicate ID the variant's `CPRA` value as its column name. If that is also a duplicate, it adds a suffix `_<n>`
  in variant order. Without the option, `read_pfile` keeps today's behaviour, and its tests do
  not change. The builder names the option, and documents it in `docs/functions/read_pfile.md`.
- **Risk, accepted by Teague:** the by-sample schema holds one column per variant in the
  cohort. Every query on the view, including `DESCRIBE`, must load every `.pvar` to bind. At
  100M variants that is very slow, or out of memory (inferred, not measured). The docs say so and
  point to `c.genotypes` with a region filter for large cohorts.

### 4. The index files (`.pvar.index.parquet`, `.psam.index.parquet`)

`CHROM`/`POS` and `CPRA` need no index: the `.pvar` is in position order, and a `CPRA` constant
becomes `POS` bounds (§2c). An index is needed only for the names that are not sorted: `ID` and
`IID`. **(Teague, 2026-11-19.)**

**Format (Teague).** One format, for a single file set and for a catalog:

| File | Columns | Sorted by |
|---|---|---|
| `<p>.pvar.index.parquet` | `ID VARCHAR, [_variant_block INTEGER,] _pvar_ordinal BIGINT` | `ID` |
| `<p>.psam.index.parquet` | `IID VARCHAR, [_sample_block INTEGER,] _psam_ordinal BIGINT` | `IID` |

- The block column is present only when more than one file is involved (a catalog). It is the
  block's position in the catalog's variant (or sample) block list. **(Teague: "a block/group
  column if relevant".)**
- `_pvar_ordinal` / `_psam_ordinal` is the row in that block's `.pvar` / `.psam` (0-based), so the
  column has one meaning in both forms. **(Claude's call.)**
- The writer writes each index with `(FORMAT parquet, PARQUET_VERSION V2, COMPRESSION zstd)`.
  **(Teague, 2026-11-19.)** A test checks `parquet_metadata(…).compression = 'ZSTD'` and
  `parquet_file_metadata(…).format_version = 2`. Measured on 1,000 sorted `rs<n>` IDs: V2 writes
  `DELTA_LENGTH_BYTE_ARRAY` / `DELTA_BINARY_PACKED` and 1,082 bytes, against `PLAIN` and 2,586
  bytes for V1.
- Row groups are small enough that one `ID =` lookup reads one row group (the builder picks the
  size and measures it). Parquet min/max statistics on the sorted key let DuckDB skip the other
  row groups.
- A `.` or empty `ID` is not indexed **(Claude's call: it is not a lookup key)**. A duplicate `ID`
  gets one index row per occurrence.

**Discovery (Claude's call).** By name, as the other companions are found: beside
`<p>.pfile.json` for a catalog, or beside `<p>.pgen` for a single file set. The catalog JSON gets
no new field. A catalog reader uses only the catalog-level index; an index beside one cell is not
read when that cell is read through its catalog. `plinking_use_parquet_companions = false` also
turns index discovery off **(Claude's call)**.

**Optional (Teague).** The writer writes an index only when a new option `index true` is set
**(Claude's call: off by default, because it costs a sort of every ID; the builder measures that
cost at 400k variants, and Teague may make it the default)**. It writes one index pair at the
output prefix: the two-column form for a single file set, the block-column form for a catalog. It
writes no index per block. A reader uses an index only when it finds one; without one, every
query gives the same rows, only slower.

**Old file sets and catalogs.** The docs give the `COPY` that builds each index from `read_pvar` /
`read_psam` (with the block number, for a catalog).

**The read side.** In `read_pfile` (so every view gets it), with an index present:

1. A pushed `ID = const` / `ID IN (consts)` queries `.pvar.index.parquet`, which returns the
   `(block, ordinal)` pairs. Only those variant blocks are kept, through the `CHROM`/`POS` block
   drop's mechanism, and the `.pgen` read narrows to the ordinals where the read path allows it
   (the #120 `WHERE ID` narrowing is the model).
2. A pushed `IID = const` / `IID IN (consts)` queries `.psam.index.parquet`, and keeps only those
   sample blocks (where `WHERE IID` already drops sample blocks today: genotype orient, and sample
   orient without a filter parameter). For a single file set it gives the sample ordinals, so the
   `.psam` need not be parsed for the lookup when the read path allows it (the builder states
   whether it does).
3. A name that is not in the index gives zero rows and opens no block.

The pushed filter still runs on the rows that are read, so a stale index cannot add rows. It can
lose rows; the docs say that the writer owns the index. **(Claude's call: no index verification
at read time, as the catalog's bounds rule. The builder decides whether a cheap check, such as
the index's row count, is enough.)**

The index is read through the VFS, so it works for a remote file set or catalog. A remote cold
`ID =` lookup should cost: the JSON, one index footer and row group, and the one block's header
and records. The PR measures this over the range server (`LAT_MS=40`) on the 20 x 2 matrix.

### 5. Pushdown of the three target reads

Each read below must be one query on a view, and must open only the blocks named:

| Read | Query | Blocks opened |
|---|---|---|
| Region, all samples | `FROM c.genotypes WHERE CHROM = '1' AND POS BETWEEN a AND b` | the variant blocks that overlap |
| Region, sample subset | `… AND IID IN ('x', 'y')` | the overlapping variant blocks × the sample blocks that hold x, y |
| Point by ID (Q1) | `FROM c.genotypes WHERE ID = 'rs123'` | the one variant block in the index × all sample blocks |
| Non-ref in a region (Q2) | `… WHERE CHROM = '1' AND POS BETWEEN … AND genotype > 0` | the overlapping variant blocks; hom-ref calls are not emitted |
| Hom-alt IIDs for IDs (Q3) | `SELECT DISTINCT IID … WHERE ID IN (…) AND genotype = 2` | the variant blocks in the index for those IDs |
| Region, wide | `SELECT ID, x, y FROM c.genotypes_by_variant WHERE POS BETWEEN …` | the overlapping variant blocks; sample blocks are dropped by projection if the read path allows it (the builder measures, and states if it does not) |

## Approaches considered

- **Shape of `genotypes_by_variant`:** STRUCT keyed by IID, wide columns, a LIST plus an IID
  list, or a MAP. Teague first picked STRUCT. The measurement showed 4.3 s per query at 1M
  samples (wide columns: 0.9 s), and he changed to wide columns.
- **Shape of `genotypes_by_sample`:** a LIST plus a `c.variants` ordinal (recommended, because it
  scales), wide columns named by ID, or wide columns under a size cap. Teague picked wide columns
  named by ID, with the schema-size risk stated.
- **ID lookup:** an ID index file (picked), a per-session ID → block map built on the first query
  (no format change, but the first query stays slow, and costly when remote), or both.
- **Index layout:** a catalog-level `ids.parquet` named in the JSON (first draft), one index per
  file (an `ID` lookup over N blocks opens N files: about N serial round trips remote, inferred),
  or one format with an optional block column. Teague picked the last, which also adds the IID
  index and covers a single file set.

## Delivery (Claude's call: four PRs, in this order)

1. `freq := true` in genotype orient, with the freq-correctness rule of §2. RED test first.
2. `WHERE genotype` pushdown (§2b), and the `CPRA` column with its pushdown (§2c). Two
   separate commits in one PR, or two PRs if the diff is large (builder's call). RED test first: a timing-free check that the filter reached
   the scan (an `EXPLAIN` or a decode counter), plus the row oracle.
3. The five views (§1, §3), docs, and the breaking-change note. RED test first: the new
   schemas and a row oracle against `read_pfile` on the un-sharded source.
4. The index files: writer option, discovery, the `ID` and `IID` read-side drops, and docs. RED
   test first: a block-open count (the cell-skip tests' method) for `ID =` and `IID =`, with and
   without an index.

## Tests (each PR)

- The oracle: `EXCEPT ALL` both ways against the un-sharded source, at 1 and 4 threads (as in
  #131), on `large_example` and `sexchr_example`, for a matrix, a v1 variant-sharded cohort, a v1
  sample-sharded cohort, a bare fileset, and a holed matrix where the view supports it.
- The block-drop evidence: the existing cell-skip counters, for each row of the §5 table.
- Freq: a tall row's freq equals `plink_freq` on the source, with and without `WHERE IID`.
- Q1, Q2 and Q3 verbatim, on every catalog kind, each against the source oracle, and each
  timed at 1M samples in the PR body (local, and Q1 also over the range server).
- `CPRA`: equals `CHROM || ':' || POS || ':' || REF || ':' || ALT` on every row; `WHERE CPRA`
  gives the same rows as the matching `CHROM`/`POS`/`REF`/`ALT` filter; a bad constant gives zero
  rows; a `.pvar` with a valid `CPRA` column is used as is, and one with a wrong value is an
  error on a read of that row.
- Names: an ID of `.`, a duplicate ID, and an IID that equals `POS`; and `read_pfile` without
  the naming option still raises the duplicate-ID error.
- Holes: each of the four wide and list views raises the "not supported" error on a holed matrix
  (at `DESCRIBE`, too); the tall view equals the source oracle minus the hole cells' rows.
- List order: element *i* against row *i* of `c.samples` / `c.variants`, at 1 and 4 threads.

## Left to the builder

- Row-group size of the index, and the measurement of an `ID =` lookup, local and over the
  test range server (`sb/rangeserver3.py`, `LAT_MS=40`).
- How freq stays full-cohort under a sample-block drop (§2).
- Whether projection drops sample blocks in `genotypes := 'columns'`.
- The exact column-name collision rules, after reading what `'columns'` does today.
