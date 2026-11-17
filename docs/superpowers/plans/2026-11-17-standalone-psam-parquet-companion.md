# Standalone FORMAT psam parquet companion — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the standalone `FORMAT psam` writer parity with the standalone `FORMAT pvar` writer — a `psam_parquet` option that emits the `.psam.parquet` companion, and a drop of a stale companion on a rewrite without it, placed so a failed query can never leave a fresh `.psam` beside an old companion.

**Architecture:** Mirror the pvar writer's `target_path`/db-at-bind pattern and reuse the existing `WritePsamParquetFromIids` helper. Two deliverables: (1) emit the companion on success, (2) drop a pre-existing companion at `InitGlobal`, before the new text `.psam` is written. Writer-only; reader, bundle, catalog, and the pvar writer are untouched.

**Tech Stack:** C++ DuckDB extension, `FileSystem` (`MoveFile`/`RemoveFile`/`FileExists`), sqllogictest `.test` files, DuckDB v2.0.0 only.

**Spec:** `docs/superpowers/specs/2026-11-17-standalone-psam-parquet-companion-design.md` (read it — it carries the intent, §5 approaches, §6 the mechanism, §6a the owner decision, §7 builder's calls, §8 tests). Branch `feat/psam-parquet-standalone` off `fb7a134`.

## Global Constraints

- DuckDB **v2.0.0 only**.
- **Never commit `third_party/plink-ng`** — the CMake configure patches it (0-line diff); it shows dirty, pointer unchanged. Exclude it from every commit.
- `clang-format` clean (`make format-fix`) or Format CI reddens.
- Match surrounding conventions: capture the `DatabaseInstance` exactly as `PvarWriteBindData` does (db-at-bind); reuse `WritePsamParquetFromIids` (`src/pfile_writer.cpp:1247`); mirror `PvarCopyFinalize`'s `tmp_`-move branch keyed on `g.path != target_path`.
- Work only in `trees/feat-psam-parquet-standalone`; do not touch the shared `main` worktree.
- Build on Slurm: `srun -n 1 -c 12 --mem=24G -p ci -t 30 bash -lc 'VCPKG_TOOLCHAIN_PATH=$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake GEN=ninja make release'` (the build is pre-warmed, so incremental fits `ci`). A **new** `.test` file is discovered only at configure/build time, so build after Task 1 adds it.

## Review Focus

- **ROW_STREAM companion order must equal the text `.psam` write order.** An unordered accumulation (e.g. reusing the dup-detection `set`) silently misaligns IIDs against the text file. — Task 1 asserts the companion's row order.
- **`psam` path override.** With `psam '/p/x.psam'`, the companion is `/p/x.psam.parquet` and the stale-drop must key on the override, not the `COPY` target. — Task 1 (emit to override) + Task 2 (drop at override).
- **`force` rewrite WITH `psam_parquet` replaces the old companion**, leaving exactly one, not two files or an appended one. — Task 1.
- **Error path: an upstream query that throws after `InitGlobal` leaves no stale companion** (the §6a case; a Finalize-only drop fails this). — Task 2.
- **Zero-row input.** An empty IID set writes a valid single-column empty companion; a reader returns zero rows, not a crash. — Task 1.

---

## Task 1: `psam_parquet` emit on standalone `FORMAT psam`

**Files:**
- Modify: `src/pfile_writer.cpp` — `PsamWriteBindData` (`:1406`), `PsamCopyBind` (`:1446`), `PsamWriteGlobalState` (`:1426`), `PsamCopySink` (`:1520`), `PsamCopyFinalize` (`:1556`).
- Create: `test/sql/write_psam_parquet_standalone.test`.

**Interfaces:**
- Consumes: `WritePsamParquetFromIids(DatabaseInstance &db, const string &psam_path, const vector<string> &iids)` (`:1247`); the `pvar_parquet` option-parse shape (`PvarCopyBind:299–304`); `PvarCopyFinalize`'s `tmp_`-move branch (`:586–592`).
- Produces: `PsamWriteBindData::target_path` (string — the real final `.psam` target: `path_override` if set, else the `COPY` file path) and `PsamWriteBindData::psam_parquet` (bool), both read by Task 2. `PsamWriteBindData::db` (captured `DatabaseInstance*`). `PsamWriteGlobalState::streamed_iids` (`vector<string>`, ROW_STREAM write-order IIDs) and `::pq_on` (bool, gates companion cleanup).

- [ ] **Step 1: Write the failing test — OUT_OF_BAND emit.** In the new `.test`, write a standalone psam with the option, then read the companion directly:
```sql
statement ok
COPY (SELECT 1 AS x) TO '__TEST_DIR__/oob.psam' (FORMAT psam, samples ['s1', 's2', 's3'], psam_parquet);

query I
SELECT IID FROM read_parquet('__TEST_DIR__/oob.psam.parquet');
----
s1
s2
s3
```
(Use the repo's test-temp-dir convention as `write_psam_parquet_companion.test` does.)

- [ ] **Step 2: Run it to verify it fails.** Build on Slurm (the new `.test` needs discovery), then run the suite filtered to `write_psam_parquet_standalone`. Expected: FAIL — `psam_parquet` is an unknown option, or no companion is written so `read_parquet` errors.

- [ ] **Step 3: Implement the option + bind capture.** In `PsamCopyBind`, parse `psam_parquet` with the exact bare/true/false shape of `pvar_parquet`; set `bind->target_path = path_override.empty() ? input.info.file_path : path_override`; capture `bind->db` as `PvarWriteBindData` captures its database. Add the three fields to `PsamWriteBindData` (and its `Copy`/`Equals`).

- [ ] **Step 4: Implement the companion write in Finalize.** Add `streamed_iids` and `pq_on` to `PsamWriteGlobalState`; in `PsamCopySink` (ROW_STREAM) append each validated IID to `streamed_iids` in write order (alongside the existing dup-set insert). In `PsamCopyFinalize`, when `bind.psam_parquet`, inside a try-guard set `pq_on` then call `WritePsamParquetFromIids(*bind.db, g.path, iids)` with `iids` = `bind.iids` (OUT_OF_BAND) or `g.streamed_iids` (ROW_STREAM); if `g.path != bind.target_path`, move `<g.path>.parquet` to `<target_path>.parquet` replacing any existing (mirror `PvarCopyFinalize:586–592`). Extend `PsamWriteGlobalState::CleanupPartial` to remove `<path>.parquet` when `pq_on`.

- [ ] **Step 5: Run the OOB test to verify it passes.** Incremental build on Slurm, run the filtered suite. Expected: PASS.

- [ ] **Step 6: Add the remaining Task-1 assertions.** In the same `.test`:
  - **ROW_STREAM order** — `COPY (SELECT iid FROM (VALUES ('b'),('a'),('c')) t(iid)) TO '__TEST_DIR__/rs.psam' (FORMAT psam, psam_parquet)`; assert `SELECT IID FROM read_parquet('__TEST_DIR__/rs.psam.parquet')` is `b,a,c` (write order, not sorted).
  - **Schema guard** — `SELECT column_name, column_type FROM (DESCRIBE SELECT * FROM read_parquet('__TEST_DIR__/oob.psam.parquet'))` is exactly one row, `IID` / `VARCHAR`.
  - **force replaces** — rewrite `oob.psam` with `force, psam_parquet` and a different `samples` set; assert `SELECT count(*) FROM glob('__TEST_DIR__/oob.psam.parquet')` = 1 and the companion holds the new IIDs.
  - **Zero-row** — `COPY (SELECT … WHERE false) TO '__TEST_DIR__/empty.psam' (FORMAT psam, samples [], psam_parquet)` (or a row-stream with no rows); assert `SELECT count(*) FROM read_parquet('__TEST_DIR__/empty.psam.parquet')` = 0, no error.
  - **`psam` path override** — `COPY (…) TO '__TEST_DIR__/ignored.psam' (FORMAT psam, psam '__TEST_DIR__/ovr.psam', samples ['s1'], psam_parquet)`; assert the companion is at `__TEST_DIR__/ovr.psam.parquet` (keyed on the override, not the `COPY` target) and holds `s1`, and that no `__TEST_DIR__/ignored.psam.parquet` exists.

- [ ] **Step 7: Run the full suite to verify Task-1 assertions pass and nothing regressed.** Expected: the new file green; no other suite reddened.

- [ ] **Step 8: `make format-fix`, then commit** (exclude `third_party/plink-ng`): the writer changes + the new `.test`.

## Task 2: InitGlobal stale-drop (§6a) + error path

**Files:**
- Modify: `src/pfile_writer.cpp` — `PsamCopyInitGlobal` (`:1496`).
- Extend: `test/sql/write_psam_parquet_standalone.test`.

**Interfaces:**
- Consumes: `PsamWriteBindData::target_path` from Task 1.
- Produces: nothing for later tasks.

- [ ] **Step 1: Write the failing test — stale-drop on a clean rewrite.** Create `x.psam` + `x.psam.parquet` (standalone, `psam_parquet`), then rewrite `x.psam` via `FORMAT psam, force` **without** `psam_parquet` and a different IID set. Assert `SELECT count(*) FROM glob('__TEST_DIR__/x.psam.parquet')` = 0 (the stale companion is gone). Add an override variant: create the companion at an override path, rewrite with the same `psam '<override>'` and no `psam_parquet`, and assert `<override>.parquet` is gone (the drop keys on `target_path` = the override).

- [ ] **Step 2: Write the failing test — error-path invariant (§6a).** Create `y.psam` + `y.psam.parquet`, then run an **out-of-band** `FORMAT psam, force` rewrite whose feeding query throws after InitGlobal — e.g. a projection that casts a late row to a failing type — so Finalize is skipped. Assert `SELECT count(*) FROM glob('__TEST_DIR__/y.psam.parquet')` = 0. (Use `statement error` for the COPY, then the glob query.)

- [ ] **Step 3: Run both to verify they fail.** Run the filtered suite. Expected: FAIL — the stale companion survives (no InitGlobal drop yet).

- [ ] **Step 4: Implement the InitGlobal drop.** In `PsamCopyInitGlobal`, before the OUT_OF_BAND text write (`:1506`) and before opening the ROW_STREAM handle (`:1509`): if `bind.target_path` is non-empty and `<target_path>.parquet` exists, `fs.RemoveFile` it. Unconditional of `psam_parquet` — when on, the fresh companion is written at Finalize; when off, none is written. This runs for both modes, so the drop never depends on Finalize running.

- [ ] **Step 5: Run both tests to verify they pass.** Incremental build, run the filtered suite. Expected: PASS — both globs return 0.

- [ ] **Step 6: Run the full suite.** Expected: the new file green (all Task-1 + Task-2 assertions); no regression across the other suites. Read the tail; record counts.

- [ ] **Step 7: `make format-fix`, then commit** (exclude `third_party/plink-ng`).

## Final: whole-branch review

- [ ] **Whole-branch review** on the most capable model, against §5 of the spec (both design choices) and this plan's Review Focus. Confirm: the ROW_STREAM companion order equals text order; the `psam` override routes both the emit and the drop to `<override>.parquet`; a `force`+`psam_parquet` rewrite leaves exactly one companion; the error-path leaves none; zero-row is valid. Fix Critical/Important in one pass, each reproduced by a failing test first.
- [ ] **STOP. Do not push or open a PR.** Report the review outcome and the suite counts, and ask the owner before `ship`. (This mirrors #122: build through green + review, then the push/PR is the owner's separate call.)
