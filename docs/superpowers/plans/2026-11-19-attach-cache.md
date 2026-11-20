# The ATTACH cache (Part 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A query on an attached catalog reuses the catalog grid, the companion paths it discovered, and the parsed `.psam` from earlier queries on the same ATTACH, so a warm remote query makes no metadata round trips for them. `ATTACH … (CACHE false)` keeps today's behaviour.

**Architecture:**
- **The cache object:** `PfileAttachCache`, an `ObjectCacheEntry` with an internal mutex, one per ATTACH. It lives in DuckDB's `ObjectCache` under a key unique to that ATTACH.
- **Lifetime:** `PfileCatalogAttach` creates the entry. The view generator, which the attached catalog owns, deletes it in its destructor, so `DETACH` frees it.
- **The views:** every `read_pfile` call in the view SQL passes `attach_cache := '<key>'`.
- **The reader:** `read_pfile` resolves the key at bind (`ObjectCache::Get`). A missing entry, for example one evicted under memory pressure, means no cache. Three load points check it first and fill it on a miss:
  1. `LoadCatalogGrid`;
  2. `LoadPfileSource`'s companion discovery;
  3. the `.psam` loaders.

**Tech Stack:** C++17 DuckDB extension (DuckDB v2.0; `duckdb/storage/object_cache.hpp`), sqllogictest.

**Spec:** `docs/superpowers/specs/2026-11-19-where-region-pushdown-and-attach-cache-design.md`, Part 2. Approved by Teague 2026-11-19, with "Build our cache, document settings" after the `enable_http_metadata_cache` measurement.

## Global Constraints

- Spec: the ATTACH is a snapshot. The cache is never re-validated against the files.
- Spec: `CACHE false` passes no key. A direct `read_pfile` call (no `attach_cache`) is never cached.
- Spec: variant metadata is not cached.
- Teague (2026-11-19): document `enable_http_metadata_cache` and `parquet_metadata_cache` as the remote settings. **Do not set them.** Measured, warm remote query, 40 ms latency:
  - catalog: 0.426 s → 0.183 s with `enable_http_metadata_cache`;
  - bare file set: 0.518 s → 0.310 s.
- The `CACHE` option must be removed from `options.options` before `PfileCatalogAttach` returns. The backing DuckCatalog's storage manager rejects unknown options (`duckdb/src/storage/storage_manager.cpp:147`).
- **Thread safety:** two queries on one ATTACH can bind at once, so every cache read and write takes the entry's mutex. Readers get copies, not references.
- **Stacked:** the branch `feat/attach-cache` is stacked on `feat/where-region-pushdown` (PR #135).
- **Golden rule, trailers, never staging `third_party/plink-ng`, never bypassing hooks, and the tools:** as in `docs/superpowers/plans/2026-11-19-where-region-pushdown.md`.

## Review Focus

1. **A file that changes after its first read.** With the cache, the old content is served; without it (`CACHE false`, or a fresh ATTACH), the new content. Task 2 pins every cached item this way.
2. **Two ATTACHes of one catalog, and `DETACH` then re-ATTACH under the same name.** They must not share an entry. Task 1 makes the key unique; Task 2 pins it.
3. **A cache entry evicted by DuckDB's LRU between queries.** `Get` returns nothing, and the query reads from disk. Task 1: the code treats a miss as "no cache".
4. **Concurrent queries on one ATTACH.** No data race; each bind takes a copy. Task 2 runs the guard queries at 4 threads. Concurrent binds are inferred safe from the mutex, not tested.
5. **The copy cost of a large `.psam`.** A 1M-sample text `.psam` cached as rows is copied on every hit. Task 3 measures the hit against the parse.

---

### Task 1: The cache object, the ATTACH option, and the view key

**Files:**
- Create: `src/include/pfile_attach_cache.hpp`, with `class PfileAttachCache : public ObjectCacheEntry`:
  - `static string ObjectType()` returning `"plinking_attach_cache"`;
  - `GetObjectType`, and `GetEstimatedCacheMemory` (a running byte estimate);
  - typed `Get…` / `Put…` pairs for the three items (Interfaces below).
- Modify: `src/pfile_catalog_storage.cpp`:
  - parse `CACHE` from `options.options` (BOOLEAN; default true; a non-boolean value is an error), then erase it;
  - build the key `"plinking_attach:" + <UUID>` (`UUID::ToString(UUID::GenerateRandomUUID())`);
  - `Put` a new `PfileAttachCache` into `ObjectCache::GetObjectCache(context)`;
  - add `, attach_cache := '<key>'` to every `read_pfile(` call that `GenotypesSql`, `MatrixGenotypesSql` and the `catalog :=` views emit;
  - give `PfileCatalogViewGenerator` the key and `DatabaseInstance &`, so its destructor can call `ObjectCache::GetObjectCache(db).Delete(key)`.
- Modify: `src/pfile_reader.cpp`: register the named parameter `attach_cache` (VARCHAR). Bind resolves it to `shared_ptr<PfileAttachCache> PfileBindData::attach_cache` (null when absent or evicted).
- Test: `test/sql/attach_pfile_catalog_cache.test` (new), the first part:
  - `ATTACH 'test/data/large_example' AS c (TYPE pfile_catalog, CACHE false)` works, and so does `CACHE true`;
  - `CACHE 'x'` is an error naming `CACHE`;
  - `c.genotypes` gives the same rows in both cases;
  - `SELECT sql FROM duckdb_views() WHERE database_name = 'c' AND view_name = 'genotypes'` contains `attach_cache :=` only for the cached ATTACH.

**Interfaces:**
- Produces, in `pfile_attach_cache.hpp`:
  - `bool GetGrid(const string &catalog_path, CatalogGrid &out) const;` / `void PutGrid(const string &catalog_path, const CatalogGrid &g);`
  - `struct CachedCompanions { string pgen, pvar, psam, eff_prefix; };` with `bool GetCompanions(const string &prefix, CachedCompanions &out) const;` / `void PutCompanions(const string &prefix, const CachedCompanions &c);`
  - `bool GetSamples(const string &kind, const string &path, SampleInfo &info, PfileSampleMetadata *meta) const;` / `void PutSamples(const string &kind, const string &path, const SampleInfo &info, const PfileSampleMetadata *meta);`
    - `kind` is the loader: `"full"` (LoadPfileSampleMetadata), `"source"` (…FromSource), `"info"` (LoadSampleMetadata) or `"count"` (LoadSampleCount).
    - `PfileSampleMetadata` moves to a header if it is file-local today; the builder states where it went.
- Produces: `PfileBindData::attach_cache` (`shared_ptr<PfileAttachCache>`).

- [ ] Step 1: write the test's first part; run it; confirm it fails at the first `CACHE false` ATTACH ("Unrecognized option for attach"); commit (tests only).
- [ ] Step 2: implement; build; run the new test, `attach_pfile_catalog.test` and `attach_pfile_catalog_kinds.test`. Expected: `All tests passed`.
- [ ] Step 3: run the suite and the format check, then commit: `feat(attach): CACHE option and a per-ATTACH cache entry that the views name`.

### Task 2: Use the cache at the three load points

**Files:**
- Modify: `src/pfile_reader.cpp`:
  - the `LoadCatalogGrid` call (about line 1956): `GetGrid`, else load and `PutGrid`;
  - `LoadPfileSource`'s discovery block (about lines 1064–1130), only when all three overrides are empty: `GetCompanions(prefix)`, else discover and `PutCompanions`. This needs `PfileAttachCache *` as a new trailing parameter, default `nullptr`;
  - the `.psam` loads in `PfileBind` and `LoadConcatSampleMetadata`: wrap each loader call in a `GetSamples`/`PutSamples` pair keyed by kind and path. `CodedSexesFromMetadata` runs on the result, so a cached `meta` gives the same SEX.
- Test: `test/sql/attach_pfile_catalog_cache.test`, the second part. Each case:
  1. ATTACH;
  2. run one query, so the cache fills;
  3. overwrite the file with garbage through `COPY (SELECT 'garbage') TO … (FORMAT csv, header false)`;
  4. then:

  | Item | Fixture | `CACHE true`, the same ATTACH | `CACHE false` / a new ATTACH / after `DETACH` + `ATTACH` |
  |---|---|---|---|
  | grid | a 3 × 2 matrix (`read_pfile_catalog_iid_prune.test:17-19`); overwrite `mx.pfile.json` | the same rows | an error |
  | companions | a text-only copy of `large_example` (writer `FORMAT pfile`); add a garbage `x.pvar.parquet` | the same rows (the text `.pvar` stays chosen) | an error (discovery prefers the parquet) |
  | `.psam` | the matrix; overwrite `GROUP=hi/samples.psam` (the writer's path; check it with `read_pfile_catalog_layout`) | the same rows | an error |

  Run the `CACHE true` queries at `SET threads = 4` too. The oracle for "the same rows" is a `CREATE TABLE … AS` of the first query.
  Use a separate fixture directory per row, so that a garbage write cannot leak into another case.

- [ ] Step 1: write the second part; run it; confirm the `CACHE true` cases fail (they error today); commit (tests only).
- [ ] Step 2: implement; build; run the file and the full suite. Expected: `All tests passed`.
- [ ] Step 3: update the docs:
  - `docs/formats/pfile_catalog.md`, the "Attaching a catalog" section: the `CACHE` option, the snapshot rule, and the two DuckDB settings with their measured effect;
  - `CLAUDE.md`, one bullet.

  Then run the format check and commit: `feat(attach): queries on one ATTACH reuse its catalog grid, companion paths and .psam`.

### Task 3: Measure, for the PR body

Run these with the branch build, warm and cold, with and without `CACHE false`, with and without `enable_http_metadata_cache`:
- [ ] `…/att2/warmcat2.sh`, the 20 × 2 catalog;
- [ ] `…/rg/warm2.sh`, the 10M bare file set.

Also measure locally:
- [ ] a cache hit on a 1M-sample text `.psam` against its parse (Review Focus 5). The fixture is `…/pqcmp/v2/m` with `plinking_use_parquet_companions = false`, `ATTACH` with and without `CACHE false`, and the second query timed.

## After the last task

Push. Open the PR with base `feat/where-region-pushdown` (stacked on #135), or `main` if #135 has merged. Do not merge.
