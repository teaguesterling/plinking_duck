# WHERE genotype pushdown (views PR 2a) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** In genotype orient, a `WHERE` on `genotype` becomes the scan's `GenotypeRangeFilter`, the filter behind `include_genotypes`. The scan then skips a variant with no matching call before decoding it, and never emits a non-matching row.

**Measured gap** (1M samples): `WHERE genotype > 0` takes 2.17 s, against 1.05 s with `include_genotypes := ['het','hom_alt']`. `WHERE genotype = 2` over 3 IDs takes 0.80 s, against 0.45 s.

**Architecture:** `PfilePushdownComplexFilter` gains a third collector, `CollectGenotypeAllowed`. It walks the AND-conjunction for predicates on the `genotype` column (`geno_genotype_col`) and turns each one it recognizes into an allowed set over {0, 1, 2, missing}, intersecting them. The result is intersected into `bind_data.genotype_filter`. The pushed filter is never removed, so a predicate that is not recognized only narrows less.

**Spec:** `docs/superpowers/specs/2026-11-19-attach-genotype-views-design.md`, §2b. Approved by Teague 2026-11-19. It is split from §2c (CPRA, PR 2b; the builder's call).

## Global Constraints
- §2b: a pushed filter whose allowed set is a subset of {0, 1, 2, missing} compiles to `GenotypeRangeFilter`.
- §2b: the pushed filter still runs on the output.
- §2b: a filter that cannot be compiled stays a plain residual filter.
- §2b: with `include_genotypes :=` as well, the two intersect.
- §2b: freq and `af_range` / `ac_range` keep counting every call.
- The scope is genotype orient, not phased, not dosages: there the `genotype` column is a TINYINT, and missing is NULL (`docs/genotype-encoding.md:120`). Sample and variant orient have no scalar `genotype` column.
- **Mapping** (a missing call is NULL):

  | Predicate | Allowed set |
  |---|---|
  | `= c` | {c} |
  | `<> c` | {0,1,2} \ {c} |
  | `<`, `<=`, `>`, `>=` | the subset of {0,1,2} |
  | `BETWEEN` | the intersection of its two bounds |
  | `IN (…)` | the listed values that are in {0,1,2} |
  | `IS NULL` | {missing} |
  | `IS NOT NULL` | {0,1,2} |

  A constant that is not an integer, or any other form (OR, NOT, a function, a cast on the column), gives no constraint.
- **The hook builds `count_filter_subset`** when `has_sample_subset` holds and the subset is not built yet, as bind does at `src/pfile_reader.cpp:2714`. It sets the filter before `PruneSamplesByIid`, which reads `genotype_filter.active`.
- **Test-only counter**, in the `plinking_debug_*` pattern: `plinking_debug_genotype_where_pushed()` counts the reads whose `WHERE` compiled into the filter. `plinking_debug_genotype_where_pushed_reset()` resets it.
- Golden rule, trailers, the submodule, hooks and the tools: as in `docs/superpowers/plans/2026-11-19-attach-cache.md`.
- This branch is stacked on `feat/attach-cache` (#136).

## Review Focus
1. **A predicate that is too narrow loses rows.** Every recognized form is checked against a table oracle, including the missing (NULL) semantics: `<> 0` excludes missing, and so does `IN`.
2. **A cast on the column** (`genotype > 0.5`, `genotype::INT = 1`) must give the same rows as the oracle, pushed or not.
3. **freq values and `af_range` with a pushed filter** stay the unfiltered counts.
4. **A pushed filter together with `samples :=` and `WHERE IID`** on a sample-sharded matrix: `count_filter_subset` must exist and be rebuilt.
5. **`OR` and `NOT`** push nothing, and the counter shows it.

### Task 1: RED — the counter and the row oracle
**Files:** create `test/sql/read_pfile_genotype_where_pushdown.test`.
- **Oracle:** `CREATE TABLE t AS SELECT * FROM read_pfile('test/data/large_example', orient := 'genotype')`, plus the 3 × 2 matrix `mx` (as in `read_pfile_catalog_iid_prune.test:17-19`).
- For each predicate P below: `EXCEPT ALL` in both directions between `read_pfile(…) WHERE P` and `t WHERE P`, at 1 and 4 threads.
- **Predicates:**
  - `genotype = 0`, `= 1`, `= 2`
  - `> 0`, `>= 1`, `< 2`, `<= 1`, `<> 0`
  - `IN (1, 2)`, `IN (0, 2)`
  - `BETWEEN 1 AND 2`
  - `IS NULL`, `IS NOT NULL`
  - `> 0 AND < 2`
  - `= 2 AND include_genotypes := ['het']` (empty)
  - `> 0` with `samples := ['SAMP1','SAMP2','SAMP5']`
  - `= 2` with `WHERE IID IN ('SAMP2','SAMP6')` on `mx`
  - `= 0 OR = 2` (not pushed)
  - `NOT (genotype = 0)` (not pushed)
  - `genotype > 0.5`
  - `genotype + 0 > 0`
- Each predicate also gets the counter: reset, run the query, then read `plinking_debug_genotype_where_pushed()`. It is 1 for the pushed forms and 0 for OR / NOT. Both cast forms are measured first, and their counts written in.
- **freq:** `SELECT DISTINCT ID, ALT_FREQ, HET_CT … WHERE genotype = 2` equals variant orient's freq for those IDs.
- [ ] Step 1: write the file. Run it, and confirm it fails at the first counter assertion (the function does not exist). Commit (tests only).

### Task 2: GREEN — compile the filter in the hook
**Files:** `src/pfile_reader.cpp`, `PfilePushdownComplexFilter` and a new `static bool CollectGenotypeAllowed(ClientContext &, const LogicalGet &, const Expression &, idx_t geno_col, bool allowed[4]);`. Index 3 is missing. It returns true when it narrowed anything.
- Gate: `push_genotype = orient == GENOTYPE && !include_phased && !include_dosages && !lazy_variants_loaded`.
- Apply: if anything was narrowed, intersect into `genotype_filter` (or set it, if it is inactive), set `active`, and build `count_filter_subset` if `has_sample_subset`. Increment the counter. All of this comes before `PruneSamplesByIid`.
- Register the two debug functions next to the `plinking_debug_sparse_*` ones.
- Docs: `docs/functions/read_pfile.md`, "Filter Pushdown", one paragraph. `CLAUDE.md`, one bullet.
- [ ] Step 1: implement, build, and run the file, the `read_pfile_include_genotypes*` tests and `read_pfile_freq_genotype.test`.
- [ ] Step 2: run the suite and the format check, then commit.

### Task 3: Measure, for the PR body
The 1M-sample `…/m1/m` set, median of 5:
- [ ] `WHERE POS BETWEEN 1 AND 50 AND genotype > 0`, against `include_genotypes := ['het','hom_alt']` (before: 2.17 s against 1.05 s);
- [ ] `WHERE ID IN (3 IDs) AND genotype = 2`, against `include_genotypes := ['hom_alt']` (before: 0.80 s against 0.45 s).
