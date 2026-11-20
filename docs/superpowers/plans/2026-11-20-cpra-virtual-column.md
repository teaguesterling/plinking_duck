# The CPRA virtual column (views PR 2b) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `read_pfile` in variant and genotype orient has a virtual column `CPRA` (`CHROM:POS:REF:ALT`, VARCHAR). It is not in `SELECT *`. It is computed only when a query names it. A `WHERE CPRA = '…'` / `IN (…)` narrows the read like the same `CHROM` + `POS` filter.

**Architecture:**
- **Declare it.** `fn.get_virtual_columns` returns `{CPRA_COLUMN_ID: TableColumn("CPRA", VARCHAR)}`, with `CPRA_COLUMN_ID = VIRTUAL_COLUMN_START + 1`, in variant and genotype orient. Sample orient gets an empty map.
- **Emit it.** `PfileDefaultScan` and `PfileTidyScan` handle the id before any lookup indexed by column.
- **Push it.** `CollectCellPruneConstraints` turns a CPRA `=` or `IN` with parsable constants into `AllowChroms` + `NarrowPos`, so `PushedRegion` and the cell drop apply as they do for `CHROM`/`POS`.
- **Expose it in ATTACH.** `c.genotypes` is built as `SELECT *, CPRA FROM read_pfile(…)`, because `SELECT *` leaves virtual columns out.

**Spec:** `docs/superpowers/specs/2026-11-19-attach-genotype-views-design.md` §2c, approved 2026-11-19.

## Global Constraints and rulings
- §2c: the format is `CHROM:POS:REF:ALT`; a multiallelic ALT stays one string.
- §2c: an unparsable constant gives zero rows, not an error.
- §2c: the pushed filter still runs on every row.
- **Ruling (Claude): a virtual column, not a schema column.** A schema column would change every `SELECT *` from `read_pfile`, and `COPY` would write it into the `.pvar`. The spec's "after `ALT`" position does not apply to a virtual column.
- **Ruling (Claude): "a `.pvar` that has its own `CPRA` column" is moot.** `read_pfile` exposes no extra `.pvar` columns: the schema is fixed to CHROM, POS, ID, REF, ALT. So the file's own column is never read.
- **Ruling (Claude): a 2-part constant (`'1:3000'`) gives zero rows.** That is SQL equality. The spec's "matches the site" is dropped; `CHROM` + `POS` already express a site query.
- **Pushdown:** a constant must have exactly 4 `:`-separated parts with an integer POS to narrow anything. Any other constant narrows nothing; the residual filter then gives zero rows.
- Golden rule, trailers, the submodule, hooks, tools: as in the PR 2a plan. Stacked on `feat/genotype-where-pushdown` (#137).

## Review Focus
1. `SELECT *` is unchanged, and `DESCRIBE` does not list `CPRA`.
2. `count(*)` and empty projections are unchanged (`GetAnyColumn` still picks column 0).
3. In genotype orient, the CPRA id is above `geno_sample_col_start`. It must not be read as a `.psam` column.
4. A residual filter on CPRA is evaluated on the emitted value.
5. Sample orient rejects `CPRA` at bind.

### Task 1: RED
`test/sql/read_pfile_cpra.test`:
- Value: `CPRA = CHROM || ':' || POS || ':' || REF || ':' || ALT` on every row, in variant and genotype orient, at 1 and 4 threads.
- `DESCRIBE` has no `CPRA`. `count(*)` is unchanged.
- Sample orient: `SELECT CPRA` is a binder error.
- `WHERE CPRA = '1:10000:T:C'` gives var100. `IN ('1:10000:T:C', '1:20000:T:C')` gives 2 rows. Genotype orient: 8 rows.
- **Narrowing proof:** the drop3 companion (as in `read_pfile_where_region_pushdown.test`). `WHERE CPRA = '1:10000:T:C'` reads cleanly, where a full read fails with "variant count mismatch".
- `'1:10000'`, `'garbage'` and `'1:x:T:C'` give 0 rows and no error.
- `IN` across chr1 and chr2 gives the same rows as `CHROM`/`POS`/`REF`/`ALT` on the source.
- ATTACH: `SELECT CPRA, ID FROM c.genotypes WHERE CPRA = '1:10000:T:C'` gives var100.
- [ ] Write it; confirm it fails at the first `CPRA` reference; commit.

### Task 2: GREEN
- `src/pfile_reader.cpp`:
  - `PfileGetVirtualColumns`;
  - the scan emission, for the default scan and the tidy scan;
  - `CollectCellPruneConstraints` (`=`, `IN`);
  - the `need_genotypes` / psam-index guards for `CPRA_COLUMN_ID`.
- `src/pfile_catalog_storage.cpp`: `SELECT *, CPRA` in the variant/genotype-orient genotypes views.
- Docs: `docs/functions/read_pfile.md` (new "CPRA" subsection), `docs/formats/pfile_catalog.md` (views), `CLAUDE.md`.
- [ ] Implement; run the file, the suite and the format check; commit.
