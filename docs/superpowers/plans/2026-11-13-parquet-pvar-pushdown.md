# Region-Pruned Parquet `.pvar` Companion — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> **Target path in repo:** `docs/superpowers/plans/2026-11-13-parquet-pvar-pushdown.md` (drafted in `perf-bench/`).

**Goal:** A gene-region read of a `pfile_catalog` fetches ~one parquet row-group of variant metadata instead of the whole block's `.pvar`, by emitting a POS-sorted parquet `.pvar` companion that carries an explicit `__vidx` column and pushing the region predicate into its scan.

**Architecture:** Two seams. (1) The pfile/sharded/matrix writers gain the ability to emit a parquet `.pvar` companion (`<prefix>.pvar.parquet`) alongside (or instead of) the text one, POS-sorted, zstd, with a `__vidx` BIGINT column = the 0-based original variant ordinal per file. (2) `read_pvar`'s external-source materialization (`pvar_reader.cpp` `PvarBind`, which today runs `SELECT * FROM '<parquet>'` and buffers every row) pushes the region/`variant_range` predicate into that query when a `__vidx` column is present, and takes each surviving row's ordinal from `__vidx` instead of a post-emit counter. The existing `filter_pushdown` ban (which protects ordinal stability) is lifted ONLY on the parquet-with-`__vidx` path.

**Tech Stack:** C++ DuckDB extension (v2.0.0 line, `make release` via vcpkg+ninja); pgenlib/plink2 vendored; tests are sqllogictest `.test` files under `test/sql/` (`require plinking_duck`, `require parquet`); run via `run_tests.py ./build/release/test/unittest --test-flags --require plinking_duck test/*`.

**Spec:** `perf-bench/SPEC-parquet-pvar-pushdown.md` (artifact 9mzDiCBVwztAZrbNAhgNWY). Measured justification: `perf-bench/FINDINGS.md`.

## Global Constraints
- Extension is **v2.0.0-only**; sharded reads of Teague's cohorts go through the **concat** path.
- Parquet companions are gated by the EXISTING setting **`plinking_use_parquet_companions`**; do not add a new setting unless a task proves one is needed.
- `__vidx` column: name literally `__vidx`, type `BIGINT`, value = the **0-based original variant ordinal within that file** (resets per file, matching `read_pvar`'s `__ordinal` / read_pfile's per-shard `vidx`).
- Parquet companion write properties: **POS-sorted**, `COMPRESSION zstd`, **DuckDB default row-group size** (measured: fine row-groups fatten the footer and cost net remote bytes — do NOT set a small `ROW_GROUP_SIZE`).
- The **text `.pvar` path is unchanged** — its `filter_pushdown` ban stays. Pushdown is a parquet+`__vidx`-only capability.
- A parquet companion **without** a `__vidx` column MUST fall back to today's materialize-all (no-pushdown) read — never misalign.
- Our repo: normal commit norms; end commits with the session's `Co-Authored-By` / `Claude-Session` lines.

## Review Focus
Edge cases the spec implies but a happy-path test won't catch — each gets a test in the owning task:
- **Multi-file / sharded block (per-file vidx reset):** a region read across shards maps each surviving row to the correct per-file `__vidx` and thus the correct `.pgen` record. *(Task 2)*
- **Region split across a row-group boundary:** variants straddling two row-groups — pruning keeps both; no variant dropped at the seam. *(Task 2)*
- **Parquet companion with no `__vidx` column:** falls back to materialize-all, correct rows, no crash, no misalignment. *(Task 3)*
- **CHROM match + POS type:** `region.chrom` matched against the companion's `CHROM`, POS predicate pushed as integer so row-groups actually prune (a string POS would silently not prune → correctness still holds via post-filter, but assert pruning on the int path). *(Task 2)*
- **Empty-window and whole-contig region:** a region with zero variants, and a bare-contig region, both return the right set with ordinals intact. *(Task 2)*

---

### Task 1: Writer emits a parquet `.pvar` companion with `__vidx`

**Files:**
- Modify: `src/pfile_writer.cpp` (the `.pvar` emission + the combined pfile/sharded writers; options parsed from `input.info.options`, see `BuildPvarBind`).
- Modify: `src/plinking_duck_extension.cpp` if a writer option is registered there.
- Test: `test/sql/write_pvar_parquet_companion.test`

**Interfaces:**
- Produces: `<prefix>.pvar.parquet` with columns `CHROM, POS, ID, REF, ALT, __vidx` (BIGINT), POS-sorted, zstd, default row-groups. `__vidx` is the per-file 0-based variant ordinal.
- Consumes: the same variant rows the text `.pvar` writer already buffers/streams.

- [ ] **Step 1: Write the failing test** `test/sql/write_pvar_parquet_companion.test`: `COPY (SELECT … FROM read_pfile('test/data/pfile_example')) TO '{TEST_DIR}/out' (FORMAT pfile, …, <companion-parquet option>)`; then assert `read_parquet('{TEST_DIR}/out.pvar.parquet')` returns the same CHROM/POS/ID/REF/ALT as `read_pvar('{TEST_DIR}/out.pvar')`, that `__vidx` is `0,1,2,…` in POS order, and that the file is POS-sorted (`SELECT POS = list_sort(POS) …` style check).
- [ ] **Step 2: Run it — expect FAIL** (option unknown / no `.pvar.parquet` written).
- [ ] **Step 3: Implement** the parquet emission: add the writer option (decide `pvar_format := 'parquet'|'both'` vs auto-when-`plinking_use_parquet_companions`; state the choice in the commit) and emit `<prefix>.pvar.parquet` via a `COPY … (FORMAT parquet, COMPRESSION zstd)` over the buffered variant rows with a `row_number()-1 AS __vidx` after `ORDER BY POS`. Reuse the five canonical role columns.
- [ ] **Step 4: Run it — expect PASS.**
- [ ] **Step 5: Commit** (`feat(writer): emit parquet .pvar companion with __vidx`).

### Task 2: Reader pushes the region into the parquet companion, ordinal from `__vidx`

**Files:**
- Modify: `src/pvar_reader.cpp` — `PvarBind` external-source branch (`is_external_source`, the `conn.Query("SELECT * FROM '<parquet>'")` at ~line 375) + the `with_ordinal`/`materialized_ordinals` handling.
- Modify: `src/pfile_reader.cpp` — the companion read site, to pass the active `RegionFilter` (`struct RegionFilter`, `ParseRegion`) down to the companion materialization.
- Test: `test/sql/read_pfile_parquet_pvar_pushdown.test` (the ordinal-identity acceptance test).

**Interfaces:**
- Consumes: Task 1's `<prefix>.pvar.parquet` with `__vidx`; the `region`/`variant_range` predicate.
- Produces: identical rows + ordinals to the text path, having scanned only the overlapping row-groups.

- [ ] **Step 1: Write the failing acceptance test** — build a fileset with BOTH a text `.pvar` and a parquet `.pvar` companion over a block with a non-trivial layout, then assert all three of these return byte-identical rows AND identical genotypes:
  (a) `read_pfile(prefix, region:='<gene>')` with `plinking_use_parquet_companions=false` (text),
  (b) same with `=true` (pruned parquet),
  (c) `read_pfile(prefix, region:='<gene>')` forced through an unpruned parquet read.
  Include the Review-Focus cases: a region spanning a row-group boundary; an empty window; a bare-contig region; and a 2-shard fixture where `__vidx` resets per file (assert the genotypes map to the right records).
- [ ] **Step 2: Run it — expect FAIL** (today the parquet companion is read in full and ordinals come from the emit counter; results still *correct* but the test also asserts the pruned path is taken — add a probe, e.g. a `read_pfile_catalog_layout`/EXPLAIN or a row-group-count assertion, that fails until pushdown exists).
- [ ] **Step 3: Implement** — in `PvarBind`, when the source is parquet AND has a `__vidx` column AND a region/variant_range is supplied, materialize via `SELECT *, __vidx FROM '<parquet>' WHERE CHROM = '<c>' AND POS BETWEEN <lo> AND <hi>` (DuckDB prunes row-groups), and populate `materialized_ordinals` from `__vidx`, not the emit counter. Thread the `RegionFilter` from `pfile_reader.cpp`'s companion read into the pvar bind.
- [ ] **Step 4: Run it — expect PASS** (identical rows/genotypes across (a)/(b)/(c); pruned path taken).
- [ ] **Step 5: Commit** (`feat(read): push region into parquet .pvar companion, ordinal from __vidx`).

### Task 3: No-`__vidx` fallback + catalog region-prune wiring + docs

**Files:**
- Modify: `src/pvar_reader.cpp` / `src/pfile_catalog_reader.cpp` (ensure a parquet companion lacking `__vidx` uses the materialize-all path; ensure catalog region-prune reads benefit).
- Test: `test/sql/read_pfile_parquet_pvar_no_vidx_fallback.test`
- Modify: `docs/formats/pfile_catalog.md` + `docs/formats/write_pfile.md` (document the parquet companion, `__vidx`, and the `plinking_use_parquet_companions` + region recommendation).

**Interfaces:**
- Consumes: Task 2's pushdown path.
- Produces: safe behavior on a hand-made parquet companion with no `__vidx`; documented surface.

- [ ] **Step 1: Write the failing fallback test** — a parquet companion WITHOUT `__vidx` (hand-written) read with a region: assert correct rows/ordinals via the materialize-all path (no crash, no misalignment).
- [ ] **Step 2: Run it — expect FAIL or crash.**
- [ ] **Step 3: Implement** the `__vidx`-present detection + fallback; confirm `ATTACH`/catalog region reads route through the same companion path.
- [ ] **Step 4: Run it — expect PASS.**
- [ ] **Step 5: Docs** — one section each in the two doc files; a line that the gene query should use `region`, parquet companions, and (on an s3-mount) `plinking_pgen_io='vfs'`.
- [ ] **Step 6: Commit** (`feat(read): __vidx fallback + catalog wiring; docs`).

## Self-Review
- **Spec coverage:** Step A (writer emits parquet) = Task 1; Step B (pushdown preserving vidx) = Task 2; approach A (explicit `__vidx`) is the Global Constraint + Task 1 output; non-goals (`.psam`, custom index, rsid pruning) are untouched. ✓
- **Interfaces:** Task 1 Produces `__vidx` (BIGINT, per-file 0-based) == Task 2 Consumes `__vidx`. ✓
- **Review Focus:** all five edge cases have tests in Task 2/Task 3. ✓
- **Proportion:** three tasks, test-first; bodies left to the implementer (the hard decision — buffer-vs-restream in the writer, and threading the RegionFilter — are named, not coded). ✓
- **The one risk:** ordinal correctness. Task 2's acceptance test is byte-and-genotype identity across text / pruned-parquet / unpruned-parquet, written RED first. Nothing merges until green.
