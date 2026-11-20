# The ID index, read side (views PR 4a) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** a `WHERE ID = '…'` / `ID IN (…)` on `read_pfile`, or on any ATTACH view, uses `<p>.pvar.index.parquet` when the index exists.
- It opens only the variant blocks that hold those IDs.
- With a `.pvar.parquet`, it loads only those rows.

Measured gap: at 10M variants, `WHERE ID =` takes 3.4 s locally and 6 s remote. At 400k variants it takes 0.129 s, against 0.013 s for `WHERE CPRA`.

**Spec:** `docs/superpowers/specs/2026-11-19-attach-genotype-views-design.md` §4, approved 2026-11-19.

## Format (spec §4, verbatim)

| File | Columns | Sort |
|---|---|---|
| `<p>.pvar.index.parquet` | `ID VARCHAR, [_variant_block INTEGER,] _pvar_ordinal BIGINT` | `ID` |

- The block column is present only for a catalog: the block's position in the catalog's variant-block list.
- `_pvar_ordinal` is the row in that block's `.pvar`, counted from 0.
- `.` and empty IDs are not indexed.
- Discovery is by name: beside `<p>.pfile.json` for a catalog, beside `<p>.pgen` for a file set. `plinking_use_parquet_companions = false` turns it off.
- A stale index can lose rows, but never adds any: the `WHERE` still runs on every row.

## Rulings (Claude)
- **Split:** this PR is the read side. The writer option `index true` is PR 4b. The docs give the SQL recipe that builds an index from `read_pfile(…, with_ordinal := true)`, and the tests use it.
- **The `.psam` index's read side is deferred.** Bind already holds every IID for the `WHERE IID` block drop, so there is no gain to measure yet.
- **With an active `region :=`** (explicit or pushed), the ordinal load is skipped and the region load is used. The ID filter still applies to every row.
- **A text `.pvar`** keeps its full parse; only the block drop applies.
- **Constraints:** `ID = 'x'` and `ID IN ('x', …)` with VARCHAR constants only. Anything else gives no ID constraint.

## Design
- `CellPruneConstraints` gets `has_id` + `ids`, collected for `ID_COL` in `CollectCellPruneConstraints`.
- **In `PfilePushdownComplexFilter`,** when `lazy_variants && !loaded && has_id`:
  1. Find the index.
  2. Query `SELECT [_variant_block,] _pvar_ordinal FROM read_parquet(idx) WHERE ID IN (…)`.
  3. Build a keep-set of blocks. A catalog's `sources[i]` is variant block *i* before any drop.
  4. Store `id_ordinals` (a `vector<vector<uint32_t>>` parallel to `sources`, plus a flag).
  5. Apply the keep-set together with the existing CHROM/POS drop, so `sources`, `source_ranges` and `id_ordinals` stay parallel. Source 0 is kept if no block is hit.
- **`EnsureLazyVariants`:** for a source with ordinals, a parquet `.pvar` and no active region, it calls `LoadVariantMetadataFromParquetOrdinals(context, path, ordinals, raw_variant_ct, "read_pfile")`. That is a sibling of the region loader: `WHERE __vidx IN (…)`, a sparse index, and `BuildChromOffsets` over its rows.

## Review Focus
1. **An ID in the index with a stale ordinal:** the row is dropped by the WHERE (wrong ID), or an out-of-range error comes from ingestion. Never a wrong row.
2. **A duplicate ID** (several ordinals), and IDs across several blocks.
3. **Index absent, or companions off:** the behaviour is today's.
4. **ID + CHROM/POS + region together.**
5. **A prepared statement executed twice.**

### Task 1: RED
`test/sql/read_pfile_id_index.test`:
- **File set:** the drop3 proof. `{A}` has a `.pvar.parquet` missing one chr3 row and an index built by the recipe. `WHERE ID = 'var100'` and `ID IN ('var100', 'var200')` read cleanly; without the index they fail with "variant count mismatch".
- **Catalog:** `shard_by 'CHROM'`, with the index built per block. A garbage `.pgen` sits in the CHROM 2 and 3 cells. `WHERE ID = 'var100'` succeeds with the index and fails without it.
- **Rows:** each case equals the same `WHERE` on the source, including a duplicate ID and an ID that is absent (0 rows).
- **Companions off:** with `plinking_use_parquet_companions = false`, the full-read behaviour returns.
- **ATTACH:** `c.genotypes WHERE ID = 'var100'` is 8 rows, with the index used.
- [ ] Write the file; confirm the RED cases fail; commit.

### Task 2: GREEN
- Implement as in Design.
- Docs: the recipe and the rules in `docs/functions/read_pfile.md` and `docs/formats/pfile_catalog.md`, plus a `CLAUDE.md` bullet.
- Measure the 10M bare file set, local and remote: `WHERE ID =` with and without the index.
- [ ] Run the suite and the format check; commit.
