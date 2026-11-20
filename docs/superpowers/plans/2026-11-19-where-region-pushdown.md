# A WHERE on CHROM/POS narrows the variant load (Part 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A `WHERE CHROM = '…' [AND POS BETWEEN a AND b]` on `read_pfile` (and so on every ATTACH view) loads only the matching variants, as `region :=` does. Measured today at 10M variants: 3.42 s against 0.10 s.

**Architecture:**
- **Task 2: postpone the load.** Variant-orient and genotype-orient reads with no `variants :=` load their `.pvar` at InitGlobal (`lazy_variants`), not only catalog reads.
- **Task 3: push the region.** `PfilePushdownComplexFilter` turns a pushed single-contig `CHROM` (+ `POS`) constraint into the bind's `RegionFilter`. The load then narrows itself through the existing region path: `LoadSourceVariants` uses `LoadVariantMetadataFromParquetRegion` for a parquet companion, and `BuildEffectiveVariantList` follows.
- **The rows do not change,** because the pushed filter still runs on every output row.

**Tech Stack:** C++17 DuckDB extension (DuckDB v2.0), sqllogictest.

**Spec:** `docs/superpowers/specs/2026-11-19-where-region-pushdown-and-attach-cache-design.md`, Part 1. It was approved by Teague on 2026-11-19.

## Global Constraints

- Spec G2: "The rows do not change: the pushed filter still runs on the output rows."
- Non-goals: `WHERE ID` (the ID index, a later PR), sample orient (its schema needs the variant count), and a partial `.pgen` index read.
- **Ruling (Claude, from the code):** `RegionFilter` (`src/pfile_reader.cpp:86`) is one contig plus one interval.
  - Only a pushed constraint with exactly one `CHROM` value becomes a region; its `POS` bounds are optional.
  - Two or more contigs, or `POS` without `CHROM`, push nothing and load as today.
  - The spec's "a contig set … loads each of them" is not built; it is a follow-up.
- A text `.pvar` keeps today's full parse under a region too: `LoadSourceVariants` region-pushes only parquet. The `WHERE` must cost what `region :=` costs; parity is the goal, not speeding up the text path.
- This branch is stacked on PR #134 (`feat/freq-genotype-orient`), because both edit `PfilePushdownComplexFilter` and `PfileInitGlobal`.
- Golden rule: implementation OR tests per commit. The RED commit comes first. Commit trailer:
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_01NSRdqdsBejcBWhYiBrD8Ft`.
- Never stage `third_party/plink-ng`. Do not bypass git hooks.
- **Tools:** as in `docs/superpowers/plans/2026-11-19-freq-genotype-orient.md`, Global Constraints:
  - build: `srun … make release -j8`;
  - one test: `srun … build/release/test/unittest "test/sql/<f>.test"`;
  - suite: `unittest "test/*"`;
  - format: `clang-format -i`, then `duckdb/scripts/format.py --check`;
  - git: `mcp__jetsam__git` with `cwd=`;
  - reading code: the Read tool.

## Review Focus

1. **A `.pvar` error is now reported when the scan starts, not at `DESCRIBE`.** A user who ran `DESCRIBE` to validate a file set no longer gets the error there. The query still reports it. Task 2 pins both halves.
2. **An explicit `region :=` together with a `WHERE` on another contig, or a disjoint interval.** The result must be zero rows, with no error and no widened load. Task 3 pins it.
3. **A `CHROM` whose case or `chr` prefix differs from the file's** (`WHERE CHROM = 'chr1'` on a `1` file). The rows must equal the unnarrowed read, which is empty. Task 3 pins it.
4. **A `WHERE` that DuckDB rewrites:** `POS >= a AND POS <= b`, `POS = x`, `CHROM IN ('1')` with one value, and `POS + 0 BETWEEN …`, which must push nothing. Task 3 pins these.
5. **A prepared statement executed twice, and two threads.** The region is applied once per bind. Task 3 pins it.

---

### Task 1: RED — the narrowing proof, plus the row guards

**Files:**
- Create: `test/sql/read_pfile_where_region_pushdown.test` (`require plinking_duck`, `require parquet`).

**Interfaces:**
- Consumes: nothing. Behaviour only.

- [ ] **Step 1: Write the tests**

Fixtures:
- `{TEST_DIR}/A`: `COPY (SELECT * FROM read_pfile('test/data/large_example', genotypes := 'struct') ORDER BY CHROM, POS) TO '{TEST_DIR}/A' (FORMAT pfile, pvar_parquet)`.
- `{TEST_DIR}/drop3.pvar.parquet`: the rows of `A.pvar.parquet`, minus the last chr3 row (`WHERE NOT (CHROM = '3' AND POS = 100000)`, in `__vidx` order), written with `ROW_GROUP_SIZE 500`.
  - **The proof:** a full load fails its row-count check ("variant count mismatch"). The region loader trusts the `.pgen` count (`plink_common.hpp:430`), so a chr1 region reads cleanly.
- `{TEST_DIR}/s1`, `{TEST_DIR}/s2`: `large_example` split by `CHROM IN ('1','2')` and `CHROM = '3'`, written with `FORMAT pfile` (struct genotypes), for a v1 list read.
- `{TEST_DIR}/mx`: the 3 × 2 matrix from `read_pfile_catalog_iid_prune.test:17-19`.

```sql
-- the fixture works: a full read fails, region := passes (both measured on the base, 6f78387)
statement error
SELECT count(*) FROM read_pfile('{TEST_DIR}/A', pvar := '{TEST_DIR}/drop3.pvar.parquet');
----
variant count mismatch

query I
SELECT count(*) FROM read_pfile('{TEST_DIR}/A', pvar := '{TEST_DIR}/drop3.pvar.parquet', region := '1:10000-20000');
----
101

-- RED: the same rows through a WHERE (today: "variant count mismatch")
query I
SELECT count(*) FROM read_pfile('{TEST_DIR}/A', pvar := '{TEST_DIR}/drop3.pvar.parquet')
WHERE CHROM = '1' AND POS BETWEEN 10000 AND 20000;
----
101
```

Also add these RED cases in the same style: genotype orient (808 rows), `POS >= 10000 AND POS <= 20000`, and `CHROM = '1'` alone (1000 rows).

Guards: they pass before and after. Each is an `EXCEPT ALL` in both directions against `region :=` (or against the unfiltered `WHERE` on the source, where no region applies), with `foreach threads 1 4`:
- `WHERE CHROM='2' AND POS BETWEEN 5000 AND 7000` vs `region := '2:5000-7000'` on: `test/data/large_example` (text `.pvar`), `{TEST_DIR}/A` (parquet), `['{TEST_DIR}/s1', '{TEST_DIR}/s2']` (v1 list), and `read_pfile(catalog := '{TEST_DIR}/mx.pfile.json')`. Use variant orient (`ID, genotypes`) and genotype orient (`ID, IID, genotype`).
- These push nothing, and their rows equal the plain `WHERE` on `test/data/large_example`:
  - `CHROM IN ('1', '2') AND POS < 500`;
  - `POS BETWEEN 100 AND 300` (no `CHROM`);
  - `CHROM = '1' OR POS = 100`;
  - `CHROM = '1' AND POS + 0 BETWEEN 100 AND 300`.
- `region := '1:100-5000'` plus `WHERE CHROM = '2'` gives 0 rows; plus `WHERE CHROM = '1' AND POS BETWEEN 3000 AND 9000` gives 21 rows (POS 3000..5000).
- `WHERE CHROM = 'chr1'` gives 0 rows.
- `PREPARE q AS SELECT count(*) FROM read_pfile('{TEST_DIR}/A', pvar := '{TEST_DIR}/drop3.pvar.parquet') WHERE CHROM = '1' AND POS BETWEEN 10000 AND 20000;`, then `EXECUTE q` twice. Each gives 101. This one is RED too.

Compute every expected count with the `region :=` query on the base before writing it in. Do not trust the numbers above without that check: they are 100-step positions (inferred).

- [ ] **Step 2: Run it; confirm the RED cases fail with "variant count mismatch" and every guard passes**

Expected: the first RED query fails, and the runner stops there. Confirm the guards in a temporary copy without the RED cases (`test/sql/zz_tmp_*.test`, deleted after), as in PR #134's Task 4.

- [ ] **Step 3: Commit (tests only)**

`test(read_pfile): a WHERE on CHROM/POS narrows the variant load like region := (RED)`

---

### Task 2: Postpone the `.pvar` load for every variant/genotype-orient read

**Files:**
- Modify: `src/pfile_reader.cpp`, the `bind_data->lazy_variants =` assignment (about line 2214), and any bind-time use of `src.variants` that the suite shows still runs for a non-catalog lazy read.
- Test: `test/sql/read_pfile_lazy_pvar_single.test` (new; Review Focus 1). Its own RED first:

```sql
statement ok
COPY (SELECT * FROM read_pfile('test/data/large_example', genotypes := 'struct')) TO '{TEST_DIR}/g' (FORMAT pfile);

statement ok
COPY (SELECT 'garbage') TO '{TEST_DIR}/g.pvar' (FORMAT csv, header false);

statement ok
DESCRIBE SELECT * FROM read_pfile('{TEST_DIR}/g');                       -- RED: errors today

statement ok
DESCRIBE SELECT * FROM read_pfile('{TEST_DIR}/g', orient := 'genotype'); -- RED: errors today

statement error
SELECT count(*) FROM read_pfile('{TEST_DIR}/g');
----
<the error text the base gives for this file at DESCRIBE; copy it from the RED run>

statement error
DESCRIBE SELECT * FROM read_pfile('{TEST_DIR}/g', orient := 'sample');   -- sample orient still parses at bind
----
<same text>

statement error
DESCRIBE SELECT * FROM read_pfile('{TEST_DIR}/g', variants := ['1:100']); -- variants := still parses at bind
----
<same text>
```

**Interfaces:**
- Produces: `lazy_variants` = `!has_variant_filter && orient ∈ {VARIANT, GENOTYPE}`, for every source kind. `source_ranges` stays catalog-only, and `PfilePushdownComplexFilter`'s source drop keeps its `source_ranges.size() == sources.size()` guard.

- [ ] **Step 1: Write `read_pfile_lazy_pvar_single.test`, run it, and confirm it fails at the first `DESCRIBE`; commit (tests only)**
- [ ] **Step 2: Remove `from_catalog &&` from the condition, build, run the new file and the full suite**

  Expected: the new file passes. Any other failure is a bind-time consumer of `src.variants` on a non-catalog path. Move that consumer to after `EnsureLazyVariants`, or exclude that read from `lazy_variants`. Record a ledger ruling for each one.

- [ ] **Step 3: Format, run the suite green, commit**

`feat(read_pfile): every variant/genotype-orient read parses its .pvar at InitGlobal`

---

### Task 3: GREEN — push a single-contig CHROM/POS constraint into the region

**Files:**
- Modify: `src/pfile_reader.cpp`, `PfilePushdownComplexFilter` (about line 6323 on 6f78387).
- Modify: `docs/functions/read_pfile.md`: the "When each file is read" paragraph (about line 287) and the `region` parameter row. Also `CLAUDE.md`, the read_pfile pushdown bullets.

**Interfaces:**
- Consumes: Task 2's `lazy_variants` for single and list reads.
- Produces: `static bool PushedRegion(const CellPruneConstraints &c, const RegionFilter &explicit_region, RegionFilter &out);`.
  - It returns `false` when nothing can be pushed: `!c.has_chrom`, or `c.chroms.size() != 1`.
  - Otherwise it fills `out` with `chrom` = the one value, `start` = `max(1, c.pos_lo)` and `end` = `c.pos_hi`, both clamped to the `RegionFilter` range.
  - With an active explicit region, it intersects: a different contig, or an empty interval, returns `false`, because the explicit region and the filter then give zero rows anyway.

- [ ] **Step 1: Implement**

In `PfilePushdownComplexFilter`:
- Add `push_region = bind_data.lazy_variants && !bind_data.lazy_variants_loaded && orient ∈ {VARIANT, GENOTYPE}`, and include it in the early-return test.
- After `CollectCellPruneConstraints`, if `push_region && PushedRegion(c, bind_data.region, r)`, set `bind_data.region = r`.
- `EnsureLazyVariants` already loads with `bind_data.region`.
- The pushed filter is not removed from `filters`, so G2 holds.

- [ ] **Step 2: Build; run `read_pfile_where_region_pushdown.test`, `read_pfile_where_prune.test`, `read_pfile_catalog_where_prune.test` and `read_pfile_parquet_pvar_pushdown.test`**

Expected: `All tests passed` for all four.

- [ ] **Step 3: Update the docs; run the full suite and the format check; commit**

`feat(read_pfile): a WHERE on CHROM (+ POS) loads only that region's variants`

---

### Task 4: Measure, for the PR body

- [ ] Run `/tmp/claude-1000/-srv-logical-workspace-home-teague-Projects-plinking-duck-main/ee19a3ea-8fb2-47a1-a90a-16a2487a8457/scratchpad/rg/loc.sh` and `…/rg/warm.sh` on the 10M-variant set (`rg/rg20480`).
- [ ] Also run them on a text-only copy (no `.pvar.parquet`).
- [ ] Report, against the design's table: local, remote cold and remote warm, plus peak memory.
- [ ] Expected: the `WHERE` rows match `region :=` (inferred until measured).

## After the last task

- Push with `mcp__jetsam__git` `["push", "-u", "origin", "feat/where-region-pushdown:refs/heads/feat/where-region-pushdown"]`.
- Open the PR with base `feat/freq-genotype-orient` (stacked on #134). Do not merge.
- Then write the short Part 2 (ATTACH cache) design from the measured warm cost, for Teague's approval.
