# The five ATTACH genotype views (views PR 3) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `ATTACH … (TYPE pfile_catalog)` exposes the five genotype views of the design.

| View | One row per | Body |
|---|---|---|
| `c.genotypes` | variant × sample (tall) | `IID`, `genotype`, `CPRA`, freq |
| `c.genotypes_by_variant` | variant | one column per IID, `CPRA`, freq |
| `c.genotypes_by_sample` | sample | one column per variant ID |
| `c.genotypes_by_variant_list` | variant | today's `c.genotypes` |
| `c.genotypes_by_sample_list` | sample | today's `c.genotypes_by_sample` |

This is a breaking change to `c.genotypes`.

**Spec:** `docs/superpowers/specs/2026-11-19-attach-genotype-views-design.md` §1, §1b, §3. Approved 2026-11-19.

## Constraints (spec, verbatim decisions)
- §1: the view SQL.
  - `c.genotypes` = `read_pfile(catalog := …, orient := 'genotype', freq := true)`.
  - `c.genotypes_by_variant` = `read_pfile(catalog := …, genotypes := 'columns', freq := true)`.
  - `c.genotypes_by_sample` = `read_pfile(catalog := …, orient := 'sample', genotypes := 'columns', <naming option>)`.
  - The `_list` views are today's SQL.
  - Each variant-row view also names `CPRA`, which PR 2b made a virtual column.
  - The non-catalog fallbacks (a bare file set) get the same shapes through `GenotypesSql`.
- §1b: on a matrix with holes, the four wide and list views raise "not supported" at bind, and `DESCRIBE` raises too. The tall view is a `UNION ALL` of per-cell tall reads; a hole gives no rows.
- §3: the naming option is view-only. With it, a `.`, empty or duplicate ID column is named by its `CPRA`, then `_<n>` in variant order. Without it, `read_pfile` keeps today's behaviour, including the duplicate-ID error.
- §1: the list views follow the order of `c.samples` / `c.variants`. A test checks element *i* against row *i* at 1 and 4 threads; if that order can drift, add an `ORDER BY`.
- §1: the docs carry the breaking-change note.
- **Rulings (Claude):**
  - **The option** is `unique_column_names := true`, a BOOLEAN, default false.
  - **The holes error:** the wide and list views reuse `read_pfile`'s own "a matrix with holes is not supported yet" error. A message that names `c.genotypes` would need a new parameter, so the docs name it instead.
  - **Freq on the holed tall view** counts each row's own sample block, as today's per-block fan-out does.
  - **Migrating the existing tests:** `<alias>.genotypes` becomes `<alias>.genotypes_by_variant_list`, and `<alias>.genotypes_by_sample` becomes `<alias>.genotypes_by_sample_list`. The holed-matrix cases change meaning on purpose: the list view errors, and the counts move to the tall view.
- Golden rule, trailers, the submodule, hooks and tools: as before. Stacked on `feat/cpra-virtual-column` (#138).

## Review Focus
1. The tall view's rows equal `read_pfile(source, orient := 'genotype', freq := true)` for every kind: a bare file set, a v1 variant-sharded catalog, a v1 sample-sharded catalog, and a matrix.
2. The wide by-variant view equals `genotypes := 'columns'` on the source, with the IID column names.
3. The by-sample view's column names, with `.` and duplicate IDs.
4. List order at 4 threads.
5. The acceptance queries Q1–Q3 on the tall view equal the source oracle.

### Task 1: RED
- New `test/sql/attach_pfile_catalog_views.test` (Review Focus 1–5). Also cover `read_pfile` with `unique_column_names`: an ID of `.`, a duplicate ID, and the option absent (still the duplicate-ID error).
- Migrate the existing tests as ruled above.
- Confirm the new file and the migrated files fail because the views or the option are missing; commit (tests only).

### Task 2: GREEN
- `src/pfile_catalog_storage.cpp`: the five views per kind.
- `src/pfile_reader.cpp`: `unique_column_names`.
- Docs: `docs/formats/pfile_catalog.md` (the view table, the breaking note, holes, the list order), `docs/functions/read_pfile.md` (the option), and `CLAUDE.md`.
- Run the suite and the format check; commit.
