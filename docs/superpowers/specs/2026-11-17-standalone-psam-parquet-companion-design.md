# Standalone `FORMAT psam` parquet companion — design

**Repo:** `teaguesterling/plinking_duck` (DuckDB extension reading PLINK-2 files in SQL).
**Branch:** `feat/psam-parquet-standalone` off `main` (`fb7a134`).
**Status:** approved (owner chose approach B, "writer for both pvar and psam", 2026-11-17;
and the InitGlobal stale-drop, §6a, after the builder check). Read-side is unchanged;
this is a writer-only build.

## 1. Intent and why

The reader prefers a `.psam.parquet` companion over the text `.psam` when
`plinking_use_parquet_companions=true` (projection, lazy IIDs, sharded order). Every
writer that emits a `.psam` must keep that companion either correct or absent — a
stale companion is read as current, and returns the wrong sample identity.

Three of the four psam-writing paths already do this:

- the `FORMAT pfile` bundle, and
- every `pfile_catalog` layout (plain / `shard_by` / `shard_samples` / matrix)

write the companion under the `psam_parquet` option, and drop a stale one on a
rewrite without the option, through `PvarMoveWithCompanion(..., drop_stale=true)`.

The **standalone `FORMAT psam` writer does neither**
(`src/pfile_writer.cpp`, `PsamCopy*`, ~1402–1567). It has no `psam_parquet` option, so
it cannot emit a companion; and on a `force` rewrite it leaves any pre-existing
`<x>.psam.parquet` untouched. A reader with companions on then returns stale samples.

The standalone **`FORMAT pvar`** writer already solved the symmetric problem for
`.pvar.parquet`: the `pvar_parquet` option to refresh the companion
(`PvarCopyBind:299`), and a finalize-time drop of a stale companion when the option is
off (`PvarCopyFinalize:575–593`, keyed on a bind-captured `target_path` that is immune
to DuckDB's `tmp_` overwrite-rename). This build brings standalone psam to that same
parity. **Standalone pvar is not changed** — it is already complete and tested
(`write_pvar_parquet_companion.test`, `parquet_companions.test`).

## 2. Goals

1. `COPY (…) TO 'x.psam' (FORMAT psam, psam_parquet)` writes `x.psam.parquet` beside the
   text `.psam` — same rows, same order — readable by the reader's companion path.
2. `COPY (…) TO 'x.psam' (FORMAT psam, force)` **without** `psam_parquet` removes any
   pre-existing `x.psam.parquet` (it would describe other samples).
3. Both psam identity modes are covered: `OUT_OF_BAND` (`samples [...]` or a struct
   genotypes column) and `ROW_STREAM` (a per-row IID column).

## 3. Non-goals and anti-goals

- No change to the standalone pvar writer, the pfile bundle, or any catalog layout —
  all already correct. Touching them is out of scope.
- No new companion format. Reuse `WritePsamParquetFromIids`
  (`src/pfile_writer.cpp:1247`) exactly as the bundle and catalog use it: one `IID`
  VARCHAR column, zstd, one default row group.
- No reader change.
- **Not a partial fix.** The option alone (without the stale-drop) still lets a reader
  read stale samples; the stale-drop alone leaves psam asymmetric with pvar. Both
  halves ship together.

## 4. Constraints

- DuckDB **v2.0.0 only**.
- **Never commit `third_party/plink-ng`** — the CMake configure patches it (0-line
  diff); it shows dirty, pointer unchanged. Exclude from every commit.
- `clang-format` clean (`make format-fix`), or Format CI reddens.
- Match surrounding conventions: the pvar writer's `target_path` + db-at-bind pattern,
  `WritePsamParquetFromIids`, and `FileSystem` `MoveFile` / `RemoveFile` / `FileExists`.
- Work only in `trees/feat-psam-parquet-standalone`; do not touch the shared `main`
  worktree.
- Build on Slurm: cold build on `batch` (`-t 90`), incremental on `ci` (`-t 30`).
- A **new** `.test` file is discovered only at configure/build time.

## 5. Approaches

- **A — stale-drop only (owner rejected).** Add just the finalize-time drop of a
  pre-existing `<x>.psam.parquet` to standalone psam; no `psam_parquet` option. ~10
  lines, closes the correctness gap, but standalone psam still cannot emit a companion,
  asymmetric with standalone pvar.
- **B — full parity (chosen).** Add the `psam_parquet` option (refresh) **and** the
  stale-drop (safety-when-off), making standalone psam symmetric with standalone pvar.
  More code — `ROW_STREAM` IID retention and db-at-bind — but it is the complete
  companion-aware behavior and matches the pvar template exactly. Owner chose B,
  "writer for both pvar and psam"; pvar already has it, so this brings psam up to the
  same level.

## 6. The decided shape

Mirror the standalone pvar writer. Concretely:

- **`PsamWriteBindData`** (`:1406`) gains: `bool psam_parquet = false`,
  `string target_path`, and a captured `DatabaseInstance *db` (the companion writer
  needs a `Connection`).
- **`PsamCopyBind`** (`:1446`): parse `psam_parquet` with the same bare/true/false
  shape as `pvar_parquet` (`:299–304`); set
  `target_path = path_override.empty() ? input.info.file_path : path_override`; capture
  the database from `context`.
- **`PsamWriteGlobalState`** (`:1426`) gains `vector<string> streamed_iids` — the
  `ROW_STREAM` ordered IID accumulation (`OUT_OF_BAND` already has the full list in
  `bind.iids`). Its `CleanupPartial` also removes `<path>.parquet`.
- **`PsamCopySink`** (`:1520`, `ROW_STREAM`): append each validated IID to
  `streamed_iids`, in addition to the existing dup-detection `set` insert. Order is
  write order, which is row order — what the text `.psam` and the reader expect.
- **`PsamCopyInitGlobal`** (`:1496`): **drop the stale companion here, before the new
  text `.psam` is written** (owner decision, §6a). If `target_path` is non-empty and
  `<target_path>.parquet` exists, remove it. This runs for both modes and before the
  `OUT_OF_BAND` text write at `:1506`, so a later upstream-query error can never leave a
  new text `.psam` beside an old companion.
- **`PsamCopyFinalize`** (`:1556`): after the text `.psam` is synced and durable, and
  only when `psam_parquet` is **on**, write the fresh companion, inside the try-guard so
  a failure is cleaned: `WritePsamParquetFromIids(db, g.path, iids)`, with `iids` =
  `bind.iids` (`OUT_OF_BAND`) or `g.streamed_iids` (`ROW_STREAM`). If
  `g.path != target_path` (a `tmp_` overwrite), move it to `<target_path>.parquet`,
  replacing any existing. The stale-drop is NOT repeated here — InitGlobal already
  cleared it, and with the option off there is nothing to write.

## 6a. Error-path invariant (owner decision, 2026-11-17)

§1's invariant is **companion correct or absent**, and it must hold on the error path,
not only on clean completion. The builder check found that `OUT_OF_BAND` writes the
complete, synced new `.psam` at `InitGlobal` (`:1506`), while a Finalize-only stale-drop
is skipped when the feeding `SELECT` errors mid-stream — leaving a new, valid `.psam`
beside the OLD `.psam.parquet`, which a reader returns silently. Standalone pvar shares
the window but leaves a PARTIAL (loud) text file; psam `OUT_OF_BAND` leaves a COMPLETE
(silent) one, the worse case, and a `psam` path override escapes DuckDB's own `tmp_`
cleanup. So the stale-drop moves to **InitGlobal, before the text write** — the drop no
longer depends on Finalize running. This diverges from "mirror pvar exactly" (§5/§6) by
design: the invariant is the authority; mirroring pvar was the means. A fully atomic
staging-and-swap (the bundle's pattern) was considered and set aside as a larger change
than the invariant requires.

## 7. What is left to the builder

- **db capture at bind.** The bundle uses `*bind.pvar_bind->db`; standalone psam has no
  `pvar_bind`, so capture the `DatabaseInstance` directly from the bind `ClientContext`
  (`DatabaseInstance::GetDatabase(context)` or equivalent). Confirm the handle the
  helper needs.
- **`OUT_OF_BAND` companion timing.** The text `.psam` is written at `InitGlobal` for
  this mode; the companion can go there or at Finalize. Prefer **Finalize**, the single
  reconciliation point that serves both modes.
- **The `tmp_` rename question.** The pvar writer handles DuckDB handing `InitGlobal` a
  `tmp_` path on overwrite. This now affects ONLY the fresh-companion move at Finalize
  (option on), not the stale-drop (moved to InitGlobal, keyed on `target_path`).
  Determine empirically (plan task 1) whether standalone psam sees a `tmp_` path. If it
  does, the Finalize move branch is required; if not, it is dead but harmless. Keep the
  branch either way, keyed on `g.path != target_path`.
- Exact fixtures and assertions (§8 sets the shape).

### Builder's calls (settled by the builder check, recorded here)

- `CleanupPartial` companion removal is gated by a `pq_on`-style flag set just before
  `WritePsamParquetFromIids`, with removal keyed on `g.path` — mirroring pvar — rather
  than the unconditional wording earlier in §6.
- The nested COPY from Finalize (the `WritePsamParquetFromIids` helper runs a COPY on a
  fresh `Connection`, as the bundle already does) is de-risked by a cheap early smoke
  step in task 1, not treated as an open risk.
- A zero-row `ROW_STREAM` companion (an empty IID-only parquet) is builder-testable, not
  a planning risk.

## 8. Tests (TDD, RED → GREEN)

New `test/sql/write_psam_parquet_standalone.test`:

1. **`OUT_OF_BAND` emit.** `COPY (…) TO 'x.psam' (FORMAT psam, samples ['a','b','c'],
   psam_parquet)`. Assert `x.psam.parquet` exists; a read with
   `plinking_use_parquet_companions=true` returns `a,b,c` in order.
2. **`ROW_STREAM` emit.** `COPY (SELECT iid FROM … ORDER BY …) TO 'y.psam'
   (FORMAT psam, psam_parquet)`. Same assertions; companion order equals row order.
3. **Stale-drop.** Write `x.psam` + `x.psam.parquet` (with `psam_parquet`), then rewrite
   `x.psam` via `FORMAT psam, force` **without** `psam_parquet`, using a **different**
   IID set. Assert `x.psam.parquet` no longer exists, and a read with companions on
   returns the **new** IIDs — not the stale ones. (Sum/row-count alone cannot catch
   this; compare the IID set.)
4. **Companion-schema guard.** Assert the companion is a single `IID` VARCHAR column, so
   a regenerated fixture that stops writing the parquet, or writes a different schema,
   fails loudly.
5. **Error-path invariant (§6a).** Write `x.psam` + `x.psam.parquet` (with
   `psam_parquet`), then run a `FORMAT psam, force` **out-of-band** rewrite whose feeding
   `SELECT` errors mid-stream (e.g. a cast that throws on a late row), so Finalize is
   skipped. Assert `x.psam.parquet` is gone (the InitGlobal drop ran) — not left stale
   beside a new text `.psam`. This is the case the owner decision exists to close; a
   Finalize-only drop fails it.

No perf target: this is a correctness and parity build, not a performance one.

---
*Anchors (`fb7a134`): standalone psam `PsamCopy*` `src/pfile_writer.cpp:1402–1567`;
`WritePsamParquetFromIids` `:1247`; pvar template `PvarCopyBind:299`,
`PvarCopyFinalize:575–593`; `pvar_parquet` field `:112`. Re-grep before editing.*
