# freq := true in genotype orient (PR 1 of 4) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `read_pfile(…, orient := 'genotype', freq := true)` adds the six `plink_freq` columns to every variant × sample row, with each variant's values over the full sample set.

**Architecture:** The variant-orient scan already computes the six values per variant (`src/pfile_reader.cpp`, `PfileDefaultScan`, "freq := true: plink_freq(counts := true)'s values"). Task 1 lifts that block into two shared functions without changing behaviour. Task 3 calls them from the genotype-orient scan (`PfileTidyScan`), once per variant, and repeats the values on each sample row. A `WHERE IID` on a sample-sharded read narrows the sample set (`PruneSamplesByIid`) and would make the counts cover only the kept samples. Task 5 therefore defers that narrowing to InitGlobal, and applies it only when no freq column is projected or filtered.

**Tech Stack:** C++17 DuckDB extension (DuckDB v2.0), pgenlib, sqllogictest.

**Spec:** `docs/superpowers/specs/2026-11-19-attach-genotype-views-design.md`, §2 ("Freq columns in genotype orient"). This plan is PR 1 of the four in its "Delivery" section. PRs 2–4 get their own plans after this one merges.

## Global Constraints

- Spec §2: "Change: allow it, and repeat the variant's six values on each of its sample rows. The values are over the cohort's full sample set, not over the rows a `WHERE IID` keeps."
- `samples := …` still narrows the freq values, as it does in variant orient today (`bind_data.has_sample_subset` → `BuildAlignedSex` and `count_filter_subset`). Only a `WHERE IID` must not.
- Sample orient keeps no freq columns: `freq := true` with `orient := 'sample'` stays an error.
- `dosages := true` with `freq := true` stays an error (existing message).
- Golden rule (CLAUDE.md): implementation OR tests per commit, never both. RED commit first.
- Commit trailer, every commit:
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_01NSRdqdsBejcBWhYiBrD8Ft`.
- The `third_party/plink-ng` submodule shows dirty on purpose. Never stage or commit it.
- Worktree: `/srv/logical/workspace/home/teague/Projects/plinking_duck/trees/fix/121-block0-double-fetch`, branch `feat/freq-genotype-orient`.
- **Tools, per step kind:**
  - Build: `srun -n 1 -c 8 --mem=20G -t 90 -w longbottom make release -j8`, run through Bash. Do not run `make` without `srun`.
  - One test file: `srun -n 1 -c 4 --mem=8G -w longbottom build/release/test/unittest "test/sql/<file>.test"`. Pass one pattern only; several give "No tests ran".
  - Full suite: `srun -n 1 -c 8 --mem=20G -t 60 -w longbottom build/release/test/unittest "test/*"`. Expect `All tests passed`.
  - Format: `clang-format --style=file -i src/pfile_reader.cpp` in place, then from `duckdb/`: `python3 scripts/format.py --all --check --workdir <worktree> --directories src test`. Expect `Passed format-check`.
  - Git: `mcp__jetsam__git` with `cwd=<worktree>` (`add -- <paths>`, `commit -m …`), not raw git through Bash. Load it with ToolSearch `select:mcp__jetsam__git`.
  - Read code: the Read tool, not `cat`/`sed`.

## Review Focus

1. **A `WHERE IID` with a freq column selected, on a sample-sharded catalog.** The freq values must equal the full-cohort values, so no sample block may be dropped. Task 4 pins it with the garbage-`.pgen` matrix.
2. **`SELECT IID, genotype FROM <view over SELECT *> WHERE IID = …`, where freq is not referenced.** The sample-block drop must still happen (the ATTACH view case). At filter-pushdown time `LogicalGet` still lists every column: `UNUSED_COLUMNS` runs after `FILTER_PUSHDOWN` (`duckdb/src/optimizer/optimizer.cpp:277` vs `:392`). Task 4 pins it.
3. **The same prepared statement executed twice.** The deferred narrowing must apply once and give the same rows both times. Task 4 pins it.
4. **chrX with no SEX in the `.psam`.** The counts are NULL and `OBS_CT` is 0, as in `plink_freq`. Task 2 pins it with `psam := 'test/data/sexchr_nosex.psam'`.
5. **`include_genotypes` / `phased := true` beside freq.** The freq values are true counts: a filtered row or a phased output must not change them. Task 2 pins both.

---

### Task 1: Lift the freq computation out of the variant scan (refactor, no behaviour change)

**Files:**
- Modify: `src/pfile_reader.cpp`. The block in `PfileDefaultScan` that starts at the comment "freq := true: plink_freq(counts := true)'s values for this variant" (about line 4260 on 4db7d3b), and the emit block "freq := true columns (before the switch…" (about line 4338).

**Interfaces:**
- Produces, in `src/pfile_reader.cpp` (file-static, above `PfileDefaultScan`):
  - `struct FreqValues { double alt_freq = 0.0; int32_t counts[PfileBindData::FREQ_COL_CT] = {}; bool alt_null = false; bool counts_null = false; };`
  - `static void ComputeFreqValues(const PfileBindData &bind_data, PfileLocalState &lstate, const PfileSource &source, uint32_t vidx, uint32_t output_sample_ct, bool hardcalls_in_buffer, FreqValues &out);`
    `hardcalls_in_buffer` = `lstate.genotype_bytes` already holds this variant's hardcalls for every output sample (today's `genotypes_read`). When it is false and the variant is not autosomal, the function decodes them, as the current code does.
  - `static void WriteFreqColumn(const PfileBindData &bind_data, const FreqValues &f, idx_t file_col, Vector &vec, idx_t row);`
    `file_col` must satisfy `bind_data.IsFreqCol(file_col)`. It applies today's NULL rules: `alt_null` NULLs `ALT_FREQ`; `counts_null` NULLs every count except `OBS_CT`.

- [ ] **Step 1: Run the freq tests to record the baseline**

Run: `srun -n 1 -c 4 --mem=8G -w longbottom build/release/test/unittest "test/sql/read_pfile_freq.test"`
Expected: `All tests passed`. Record the assertion count.

- [ ] **Step 2: Move the two blocks into `ComputeFreqValues` and `WriteFreqColumn`, and call them from `PfileDefaultScan`**

Moving code only. The bodies are the current blocks, with `freq_vals`/`freq_alt`/`freq_alt_null`/`freq_counts_null` replaced by the `FreqValues` fields.

- [ ] **Step 3: Build, then run `read_pfile_freq.test` and the full suite**

Expected: the same assertion count as Step 1, `All tests passed`; the full suite `All tests passed`.

- [ ] **Step 4: Format, then commit**

`refactor(read_pfile): lift the freq computation into ComputeFreqValues / WriteFreqColumn`

---

### Task 2: RED — freq columns in genotype orient

**Files:**
- Create: `test/sql/read_pfile_freq_genotype.test` (header `# group: [sql]`, `require plinking_duck`).
- Modify: `test/sql/read_pfile_freq.test:302-305`. The genotype-orient case becomes `statement ok`. The sample-orient case (`:307-310`) keeps `statement error`, with the expected text `freq := true needs orient := 'variant' or 'genotype'`.

**Interfaces:**
- Consumes: nothing from Task 1. It tests behaviour only.

- [ ] **Step 1: Write the tests**

The oracle is variant orient's freq columns, which `read_pfile_freq.test` already ties to `plink_freq(counts := true)`. Let `X(args)` stand for `read_pfile('test/data/sexchr_example', args)`. Each check runs under `SET threads = 1` and again under `SET threads = 4`.

```sql
-- schema: the six columns follow `genotype`, with the same types as in variant orient
query TT
SELECT column_name, column_type FROM (DESCRIBE SELECT * FROM read_pfile('test/data/sexchr_example', orient := 'genotype', freq := true))
WHERE column_name IN ('genotype','ALT_FREQ','OBS_CT','HOM_REF_CT','HET_CT','HOM_ALT_CT','MISSING_CT');
----
genotype	TINYINT
ALT_FREQ	DOUBLE
OBS_CT	INTEGER
HOM_REF_CT	INTEGER
HET_CT	INTEGER
HOM_ALT_CT	INTEGER
MISSING_CT	INTEGER

-- values: per variant, the tall rows carry exactly variant orient's values (both directions, NULLs included)
-- tall_freq(args)  := SELECT DISTINCT ID, ALT_FREQ, OBS_CT, HOM_REF_CT, HET_CT, HOM_ALT_CT, MISSING_CT FROM X(orient := 'genotype', freq := true, args)
-- var_freq(args)   := SELECT ID, ALT_FREQ, OBS_CT, HOM_REF_CT, HET_CT, HOM_ALT_CT, MISSING_CT FROM X(freq := true, args)
-- assert: count of (tall EXCEPT var) UNION ALL (var EXCEPT tall) = 0, for each args below
-- assert: SELECT count(*) FROM (SELECT ID FROM tall_freq GROUP BY ID HAVING count(*) > 1) = 0   (one value set per variant)
-- assert: SELECT count(*) FROM X(orient := 'genotype', freq := true) = variants x samples of the source
```

Run the value checks with these `args`, so that each fixture and each path is covered:
- (none): autosomes, PAR and chrX non-PAR with SEX.
- `samples := ['<first 3 IIDs of sexchr_example.psam>']`: freq over the subset, the same as variant orient with the same `samples`.
- `psam := 'test/data/sexchr_nosex.psam'`: no SEX, so the chrX counts are NULL and `OBS_CT` is 0 (Review Focus 4).
- `phased := true`: the same values as without `phased` (Review Focus 5).
- `include_genotypes := ['hom_alt']`: the DISTINCT freq rows equal variant orient's for the variants that keep a row (Review Focus 5). Also assert one known variant's `HOM_REF_CT` is the unfiltered value.

Also add `read_pfile('test/data/large_example', …)` with no args: a multi-batch autosomal file, so `PgrGetCounts` is the path.

Projection and filter:

```sql
-- only a freq column projected (no genotype column): the counts still come. round(): the sums run
-- in a thread-dependent order, so compare at 1e-6.
query I
SELECT round((SELECT sum(ALT_FREQ) FROM read_pfile('test/data/large_example', orient := 'genotype', freq := true)), 6)
     = round(8 * (SELECT sum(ALT_FREQ) FROM read_pfile('test/data/large_example', freq := true)), 6);
----
true

-- WHERE on a freq column filters whole variants' rows
query I
SELECT count(*) = 8 * (SELECT count(*) FROM read_pfile('test/data/large_example', freq := true) WHERE ALT_FREQ < 0.2)
FROM read_pfile('test/data/large_example', orient := 'genotype', freq := true) WHERE ALT_FREQ < 0.2;
----
true

-- still errors
statement error
SELECT * FROM read_pfile('test/data/sexchr_example', orient := 'genotype', freq := true, dosages := true);
----
incompatible with dosages := true
```

- [ ] **Step 2: Run the new file and `read_pfile_freq.test`, and confirm they fail for the right reason**

Expected: both fail at their first genotype-orient freq call, with `freq := true needs orient := 'variant' (the default)`. Any other failure (a typo, a wrong count) is a test bug: fix the test.

- [ ] **Step 3: Commit (tests only)**

`test(read_pfile): freq := true in genotype orient (RED)`

---

### Task 3: GREEN — bind, schema and scan for genotype orient

**Files:**
- Modify: `src/pfile_reader.cpp`:
  - the bind check "freq := true needs orient := 'variant'" (about line 1913): allow GENOTYPE; for SAMPLE, raise `read_pfile: freq := true needs orient := 'variant' or 'genotype': sample orient has no per-variant row`;
  - the variant-orient freq schema block (about line 2944);
  - the genotype-orient schema block, which ends `bind_data->geno_total_cols = names.size();` (about line 2573);
  - `PfileTidyScan`'s `load_variant` and `emit_sample_row` lambdas.
- Modify: `docs/functions/read_pfile.md`, the "Frequency Columns (`freq`)" section (line 66 onwards) and line 85. Also `CLAUDE.md`, the conventions line about freq, if one exists (check with Grep; do not add a new one).

**Interfaces:**
- Consumes: `FreqValues`, `ComputeFreqValues`, `WriteFreqColumn` (Task 1).
- Produces: `static void AppendFreqColumns(PfileBindData &bind_data, vector<Name> &names, vector<LogicalType> &return_types);`. This is the existing variant-orient block (the name-clash check plus the six columns) moved into a function. Use the element type that `names` has in that block today (`CompatMakeName` / `CompatNameStr`). Genotype orient calls it after the `genotype` column and before `geno_total_cols` is set. Variant orient calls it where the block is now.
- Produces: `PfileLocalState::FreqValues tall_freq;` (or a member of the local state that the scan owns), filled once per variant.

- [ ] **Step 1: Implement**

- In `load_variant`, when `gstate.need_freq`, call `ComputeFreqValues(…, vidx, output_sample_ct, hardcalls_in_buffer, lstate.tall_freq)`. Here `hardcalls_in_buffer` is true only on the paths that ran `GenoarrToBytesMinus9` into `lstate.genotype_bytes` for this variant: the plain and the phased paths, not the dosage path, which bind refuses with freq.
- Compute after the count-filter early return, so a skipped variant costs nothing.
- In `emit_sample_row`, before the column `switch`, call `WriteFreqColumn` for each `IsFreqCol` column, as `PfileDefaultScan` does.
- `need_freq` already reads `input.column_ids`, so a WHERE-only freq column is covered.

- [ ] **Step 2: Build; run `read_pfile_freq_genotype.test` and `read_pfile_freq.test`**

Expected: `All tests passed` for both.

- [ ] **Step 3: Update the docs**

Replace "Variant orient only" wording with:
- "In genotype orient, each row repeats its variant's six values. The values are over the read's sample set (`samples :=`), not over the rows that a `WHERE IID` keeps."
- Line 85 becomes: "`freq := true` is an error with `orient := 'sample'` or `dosages := true`, and …".

- [ ] **Step 4: Full suite and format, then commit (implementation + docs)**

Expected: the full suite `All tests passed`, `Passed format-check`.
`feat(read_pfile): freq := true in genotype orient — each row repeats its variant's values`

---

### Task 4: RED — a WHERE IID must not narrow the freq values (sample-sharded)

**Files:**
- Create: `test/sql/read_pfile_freq_genotype_iid.test`. Copy the setup of `test/sql/read_pfile_catalog_iid_prune.test:13-35`: the `mx` and `badhi` 3 × 2 matrices from `large_example`, where `badhi` has a garbage `.pgen` in every `hi` cell.

**Interfaces:**
- Consumes: Task 3's behaviour.

- [ ] **Step 1: Write the tests**

```sql
-- (a) freq selected + WHERE IID: every block is needed, so the garbage 'hi' cells are opened
statement error
SELECT IID, ALT_FREQ FROM read_pfile(catalog := '{TEST_DIR}/badhi.pfile.json', orient := 'genotype', freq := true) WHERE IID = 'SAMP2';
----
data.pgen

-- (b) and the values are full-cohort: equal to plink_freq over the source, for every variant
query I
SELECT count(*) FROM (
  (SELECT ID, ALT_FREQ, OBS_CT, HOM_REF_CT, HET_CT, HOM_ALT_CT, MISSING_CT
     FROM read_pfile(catalog := '{TEST_DIR}/mx.pfile.json', orient := 'genotype', freq := true) WHERE IID = 'SAMP2'
   EXCEPT ALL
   SELECT ID, ALT_FREQ, OBS_CT, HOM_REF_CT, HET_CT, HOM_ALT_CT, MISSING_CT FROM plink_freq('test/data/large_example', counts := true))
  UNION ALL
  (SELECT ID, ALT_FREQ, OBS_CT, HOM_REF_CT, HET_CT, HOM_ALT_CT, MISSING_CT FROM plink_freq('test/data/large_example', counts := true)
   EXCEPT ALL
   SELECT ID, ALT_FREQ, OBS_CT, HOM_REF_CT, HET_CT, HOM_ALT_CT, MISSING_CT
     FROM read_pfile(catalog := '{TEST_DIR}/mx.pfile.json', orient := 'genotype', freq := true) WHERE IID = 'SAMP2'));
----
0
```

Adjust the `plink_freq` column list to its real `counts := true` names. Read `docs/functions/plink_freq.md`; if a name differs, alias it in the query, and do not change the reader.

Guards. These pass before Task 5 and must still pass after it:

```sql
-- (c) freq NOT referenced, through a SELECT * view: the 'hi' blocks are still dropped (Review Focus 2)
statement ok
CREATE VIEW vg AS SELECT * FROM read_pfile(catalog := '{TEST_DIR}/badhi.pfile.json', orient := 'genotype', freq := true);

query I
SELECT count(*) FROM (SELECT IID, genotype FROM vg WHERE IID = 'SAMP2');
----
3000

-- (d) the same without the view
query I
SELECT count(genotype) >= 0 FROM read_pfile(catalog := '{TEST_DIR}/badhi.pfile.json', orient := 'genotype', freq := true) WHERE IID = 'SAMP2';
----
true

-- (e) prepared and executed twice: the same rows both times (Review Focus 3)
statement ok
PREPARE p AS SELECT count(*) FROM vg WHERE IID = 'SAMP2';

query I
EXECUTE p;
----
3000

query I
EXECUTE p;
----
3000
```

Also run (a) and (b) under `SET threads = 1` and `SET threads = 4`.

- [ ] **Step 2: Run the file, and confirm the failure**

Expected: (a) FAILS, because the statement succeeds: today's prune drops `hi`. (b) FAILS with a non-zero count: the values cover SAMP1–SAMP4 only. (c), (d) and (e) PASS. If (c) fails, the view already defeats the prune on main. Then stop and report it, because it is a finding about main and not part of this task.

- [ ] **Step 3: Commit (tests only)**

`test(read_pfile): a WHERE IID must not narrow freq in genotype orient (RED)`

---

### Task 5: GREEN — defer the IID sample narrowing to InitGlobal when freq is on

**Files:**
- Modify: `src/pfile_reader.cpp`:
  - `PruneSamplesByIid` (about line 6233);
  - `PfileBindData` (fields, next to the `lazy_variants*` members, about line 559);
  - `PfileInitGlobal` (about line 3719), before `EnsureLazyVariants`.
- Modify: `docs/functions/read_pfile.md`, the paragraph on the IID sample-block drop in genotype orient (find it with Grep for "IID" near "sample block"). Add one sentence: "With `freq := true`, a query that selects or filters a freq column reads every sample block, because the values cover the full sample set."

**Interfaces:**
- Produces, in `PfileBindData`:
  - `vector<uint32_t> pending_sample_keep;`
  - `bool has_pending_sample_keep = false;`
  - `bool pending_sample_keep_applied = false;`
- Produces: `static void ApplySampleKeep(PfileBindData &bind_data, vector<uint32_t> keep);`. This is the tail of today's `PruneSamplesByIid`, from `if (keep.empty())` onwards, moved into a function: it sets `sample_indices`, `has_sample_subset` and `subset_sample_ct`, and rebuilds `count_filter_subset` when `genotype_filter.active`.

- [ ] **Step 1: Implement**

- `PruneSamplesByIid` computes `keep` as now.
  - If `bind_data.include_freq` is false, call `ApplySampleKeep` at once, as today, so a read without freq does not change.
  - Else store `keep` in `pending_sample_keep` and set `has_pending_sample_keep = true`.
- In `PfileInitGlobal`, before `EnsureLazyVariants`:
  - Compute `freq_needed`: any `input.column_ids` entry with `bind_data.IsFreqCol(id)`.
  - Then take `bind_data.lazy_variants_lock`, and `const_cast` as `EnsureLazyVariants` does.
  - If `has_pending_sample_keep && !pending_sample_keep_applied && !freq_needed`, call `ApplySampleKeep(…, pending_sample_keep)`.
  - Set `pending_sample_keep_applied = true` in both cases.
- **Ruling, for the executor:** the first InitGlobal decides. A re-execution with a different projection cannot happen for one bind, because a prepared statement keeps its plan.
- **Do not move `BuildAlignedSex`.** When the narrowing applies, no freq column is read, so `aligned_sex` is not used.

- [ ] **Step 2: Build; run `read_pfile_freq_genotype_iid.test`, `read_pfile_catalog_iid_prune.test` and `read_pfile_freq_genotype.test`**

Expected: `All tests passed` for all three.

- [ ] **Step 3: Measure, for the PR body**

Measure on the 20 × 2 matrix at `/tmp/claude-1000/-srv-logical-workspace-home-teague-Projects-plinking-duck-main/ee19a3ea-8fb2-47a1-a90a-16a2487a8457/scratchpad/att2/mx`: median of 5, `srun -c 8`, `.timer on`. Time three queries:
- `SELECT count(*) FROM read_pfile('mx', orient := 'genotype', freq := true) WHERE POS BETWEEN 199000 AND 199999 AND IID IN ('per1','per150')`, without freq selected;
- the same query with `ALT_FREQ` selected;
- the same query with no `freq := true`.

Report all three times. The first and the third should match (inferred).

- [ ] **Step 4: Update the docs; run the full suite and the format check; commit (implementation + docs)**

Expected: the full suite `All tests passed`, `Passed format-check`.
`feat(read_pfile): a WHERE IID keeps every sample block when a freq column is read`

---

## After the last task

Push with `mcp__jetsam__git` `["push", "-u", "origin", "feat/freq-genotype-orient:refs/heads/feat/freq-genotype-orient"]`. Open the PR with `gh pr create`; its body states:
- the RED and GREEN commits;
- the three timings;
- the ruling in Task 5.

Do not merge. Teague merges in the GitHub UI.
