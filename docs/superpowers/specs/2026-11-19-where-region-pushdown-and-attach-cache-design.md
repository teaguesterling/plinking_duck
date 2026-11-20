# A WHERE on CHROM/POS narrows the variant load; then a capped ATTACH cache

Date: 2026-11-19. Status: design, for Teague's approval.
Decisions by Teague are marked **(Teague)**. Calls made by Claude are marked **(Claude's call)**.

## Intent

A `WHERE CHROM = … AND POS BETWEEN …` must cost what the same `region :=` costs. Today it costs
34× more locally, and it repeats on every query. After that, a cache that `ATTACH` holds keeps
the parsed per-cell data that is small enough to keep. `ATTACH` can turn the cache off.
**(Teague: "Both, but allow caching to be disabled in attach".)**

## What exists today (measured 2026-11-19; 10M variants × 8 samples, V2 `.pvar.parquet` in 20,480-row groups)

| Query, local | Time | Peak memory |
|---|---|---|
| `read_pfile(p, region := '3:30000000-30100000')` | 0.10 s | 116 MB |
| `read_pfile(p) WHERE CHROM = '3' AND POS BETWEEN 30000000 AND 30100000` (same 349 rows) | 3.42 s | 2,165 MB |
| `c.genotypes WHERE …` (the same region) | 3.48 s | |
| `c.genotypes WHERE ID = 'rs…'` | 3.35 s | |
| the same `c.genotypes` query twice in one session | 6.65 s | |

- **Remote, cold, 40 ms latency:** `region :=` read 5 MB of the `.pvar.parquet` in 58 requests. `c.genotypes WHERE …` read 59.8 MB in 433 requests and took 5.8 s. A warm repeat took 3.5 s, almost all CPU.
- **The `.pgen` header is not the cost.** Its index is 15 MB here (1.5 bytes per variant). A cold query reads it whole (`PgfiInitPhase2`). Later queries in the same session read none of it, because the external file cache from #126 keeps it.
- **Why the `WHERE` is slow.** `lazy_variants` (`src/pfile_reader.cpp:2104`) is on only for a catalog read. A single file set, and a v1 list read, parse every `.pvar` at bind. A catalog cell parses its `.pvar` at InitGlobal, but `EnsureLazyVariants` narrows it only by the explicit `region :=` (`LoadSourceVariants(context, src, bind_data.region)`). The `WHERE` bounds are used only to drop whole cells (`PfilePushdownComplexFilter`), and to narrow the `.pgen` read (#120).

## Goals

- **G1.** A pushed `CHROM`/`POS` filter narrows the `.pvar` / `.pvar.parquet` load to its bounds, as `region :=` does. This holds for a single file set, a v1 list read, and every cell of a catalog read, in variant and genotype orient.
- **G2.** The rows do not change: the pushed filter still runs on the output rows. The oracle is the same query without the narrowing.
- **G3.** A capped cache, held by the attached catalog, is on by default. `ATTACH … (TYPE pfile_catalog, CACHE false)` turns it off. **(Teague: it can be turned off. Claude's call: the option name `CACHE`, a BOOLEAN, default `true`.)**

## Non-goals

- `WHERE ID =` speed. The ID index of the views design (its PR 4) covers it.
- Sample orient. Its schema needs the variant count at bind.
- A partial read of the `.pgen` index. pgenlib reads all of it at open, and the external file cache already keeps it for the session.

## Part 1: WHERE → region narrowing (one PR, first)

1. **Defer the variant load.** Turn `lazy_variants` on for every variant-orient and genotype-orient read with no `variants :=`, not only for catalog reads.
   - The schema must still need nothing but the sample axis and the `.pvar` header (the extra column names and types).
   - A read whose schema needs more keeps the bind-time load: for example, if extra `.pvar` columns need a full parse to type them.
   - The builder finds every bind-time use of `src.variants` on those paths, and lists it in the PR. **(Claude's call.)**
2. **Collect the bounds.** At `PfilePushdownComplexFilter`, collect the `CHROM` / `POS` constraints that `CollectCellPruneConstraints` already finds, and store them in the bind data as the "pushed region": a contig set, and `POS` low and high.
3. **Apply them.** At InitGlobal, `EnsureLazyVariants` loads each surviving source with the effective region: explicit `region :=` ∩ the pushed region. This applies to `LoadSourceVariants` and `BuildEffectiveVariantList`.
   - A contig set of more than one contig loads each of them **(Claude's call)**.
   - A constraint the region form cannot express (an `OR`, a function of `POS`) pushes no region, as for cell drops today.
4. **The `.pgen` read** keeps #120's narrowing. Nothing changes there.
5. **The same once-per-bind rule as `lazy_variants`.** A prepared statement reuses the region that its first InitGlobal applied.

**Tests (RED first):**
- **Same rows:** `WHERE CHROM/POS` with the narrowing equals `region :=`. Do this for a file set with `.pvar.parquet` and for one with a text `.pvar`, for a v1 list, and for a matrix catalog, at 1 and 4 threads.
- **Unrecognized forms:** an `OR`, a `POS + 0 BETWEEN …`, and a contig set.
- **A load counter shows the narrowing.** It is a test-only, or `EXPLAIN`-visible, count of variants loaded. Timing-free. The builder picks the method and states it. **(Claude's call.)**

**Measure, in the PR body:** the table above again, local and remote, cold and warm, with a text `.pvar` and with a `.pvar.parquet`. Expected: the `WHERE` matches `region :=` (inferred).

## Part 2: the ATTACH cache (second PR)

**What Part 1 left, measured 2026-11-19** (remote, 40 ms latency, one run each):

| Warm query (second run in a session) | Time | Requests per query |
|---|---|---|
| bare file set, `c.genotypes WHERE CHROM … POS …` (10M variants) | 0.51 s | 8 HEAD (companion discovery: `.pgen`, `.psam`, `.psam.parquet`, `.fam.parquet`, 4 × `.pvar.parquet`), 1 `.psam` GET, 2-5 `.pvar.parquet` GETs |
| 20 x 2 matrix catalog, `c.genotypes WHERE POS …` | 0.44 s | `.pfile.json` HEAD + GET, 2 x `.psam` HEAD + GET, `.pvar` HEAD + GET (560 KB, again), 2 x `.pgen` HEAD |
| the bare file set with `SET parquet_metadata_cache = true` | 0.51 s | `.pvar.parquet` GETs 5 -> 2; no time saved |

So a warm query spends its time on about 10 serial metadata round trips. The `.pgen` index is
already kept by the external file cache (#126). DuckDB's footer cache does not help.

**The cache.** The attached catalog owns one cache object, held in DuckDB's `ObjectCache` under a key
unique to that ATTACH and dropped at DETACH. The views pass that key to `read_pfile` (a new
internal named parameter, `attach_cache := '<key>'`). `read_pfile` then takes from the cache, and
fills it on a miss:

1. the parsed catalog grid (`LoadCatalogGrid`): no `.pfile.json` HEAD or GET;
2. the companion paths that discovery resolved, for each cell or file set: no HEAD requests;
3. each `.psam`'s parsed sample metadata (`SampleInfo` and `PfileSampleMetadata`): no `.psam` GET
   or parse.

**Not cached in this PR:** the variant metadata. A full load is about 2 GB per 10M variants
(inferred from a 2,165 MB peak), and with a `.pvar.parquet` the region load is already small. A text
`.pvar` stays a per-query parse (the 560 KB GET above). A capped variant cache is a possible
follow-up.

**Validity.** The ATTACH is a snapshot: the cache is never re-validated against the files. A file
that changes after ATTACH is seen after `DETACH` and `ATTACH`, or with `CACHE false`. **(Claude's
call; validating by mtime would cost one HEAD per file, which is the cost being removed.)**

**`ATTACH … (TYPE pfile_catalog, CACHE false)`** passes no key, so every query reads and parses
its own metadata, as today. **(Teague.)** A direct `read_pfile` call is not cached.

**Expected:** a warm query does only its data reads (the `.pgen` region records), so it takes about
0.05-0.1 s instead of 0.44-0.51 s (inferred). Memory: the grid and the paths are small. A `.psam` of
1M samples is tens of MB (inferred, to be measured).

**Tests:** the same rows with and without `CACHE false`. A request-count check that a warm query
on the cached ATTACH opens no `.pfile.json` or `.psam`: through the range server, or a test-only
counter, the builder picks. `DETACH` frees the entry. Two ATTACHes of one catalog do not share a
key.

## Order with the views work

Part 1 comes before PR 1 of `2026-11-19-attach-genotype-views-design.md`. Both edit
`PfilePushdownComplexFilter` and InitGlobal in `src/pfile_reader.cpp`. The PR 1 plan
(`docs/superpowers/plans/2026-11-19-freq-genotype-orient.md`) is re-checked against Part 1's merged
code before it runs.
