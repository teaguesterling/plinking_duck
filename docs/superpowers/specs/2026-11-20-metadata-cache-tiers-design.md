# Metadata cache tiers for an 8M-sample × 170M-variant cohort — design

Status: **approved** (Teague, 2026-11-20; §6 answered). Author: Claude Opus, from Teague's request:
"caching all metadata may be impossible on smaller memory machines … smart caching or caching tiers
that don't take many gigabytes", and "is it worth creating a single global (sync'd) pgen header cache?"

## 1. The target, and the layout it implies

- **Cohort:** 8 million samples × 170 million variants (Teague).
- **Layout:** the variants are split into several chunks per chromosome, and the samples into groups by
  ancestry (Teague). That is a v2 matrix: one `.pgen` per (chunk, group) cell.
- **Consequence:** each cell's `.pgen` header covers one chunk, not 170M variants. Cells of the same chunk
  in different groups have **different** headers, because record lengths depend on the samples. So a
  header is cached per cell file.

## 2. What costs memory today (measured 2026-11-20, longbottom)

| Item | Size | Held by | Evidence |
|---|---|---|---|
| `.pgen` header, decoded (`vrtypes` + `var_fpos`) | 9 B per variant | **each scan thread**, per cell it opens | probe: 90.0 MB and 52 ms per open of a 10M-variant `.pgen`; +88 MB RSS per open |
| ATTACH sample metadata (`.psam`) | ~103 B per sample (IID only) | the ATTACH cache, **copied out on every query** | 103 MB at 1M samples |
| Variant metadata (`.pvar`) | ~200 B per variant | per query, never cached | 2.1 GB peak for a 10M-variant full scan |
| Catalog grid, companion paths, parquet footers | KB–MB | ATTACH cache / DuckDB | — |

**Per-thread header copies, measured on `read_pfile`** (10M-variant file set, `count(genotypes)`):

| Query | 1 thread | 16 threads |
|---|---|---|
| chr1 region | 306 MB | 1479 MB |
| full scan with `freq` | 2168 MB | 3519 MB |

**One shared header, measured with a pgenlib probe** (`scratchpad/hdr/share.cc`: N threads read every
variant with `PgrGet`, checksums equal):

| Threads | private header per thread (today) | one shared header |
|---|---|---|
| 1 | 86 MB | 86 MB |
| 4 | 344 MB | 86 MB |
| 16 | 1377 MB | 86 MB |

Time was the same in both modes (447 / 1720 / 1055 ms against 444 / 1676 / 1046 ms).

**The cohort-shaped fixture** (10 chunks of 1M variants × 5 groups = 50 cells, from the 10M-variant
file set; written by `COPY … (FORMAT pfile, shard_by, shard_samples, index)` in 95 s), today's code:

| Query (`c.genotypes`, ATTACH) | 1 thread | 16 threads |
|---|---|---|
| one ID (5 cells) | 49 ms warm, 71 MB | 49–74 ms, 80 MB |
| 100 kb region | 87 ms, 97 MB | 71 ms, 120 MB |
| chr8, `genotype = 2` | 3.0 s, 288 MB | 0.9–1.0 s, 1108 MB |
| full scan with `freq` | 8.0 s, 1521 MB | 3.4 s, 2283 MB |

The chr8 query's +820 MB at 16 threads fits 16 threads × 5 groups × ~9 MB per 1M-variant cell header
(inferred; not separately attributed).

**At the target (arithmetic, not measured):** one header per cell is 9 B × chunk variants. All cells
together are 9 B × 170M × groups: 7.7 GB for 5 groups. Today each of 16 threads can hold its own copy of
each cell it touches.

## 3. pgenlib already supports one shared header

- We open `.pgen` files in pgenlib's per-variant mode (mode 2). There, "multiple PgenReaders can be based
  off of" one `PgenFileInfo`; each later reader opens its own `FILE*` from `fname`
  (`pgenlib_read.h:478-505`).
- After `PgfiInitPhase2`, the header arrays are only read. So threads and queries can share them.
- The probe closes the header's `shared_ff` after phase 2, so every reader opens its own handle. Each
  reader keeps its own `PgenReader` and buffers (small: 0.5 MB at 1M samples).
- **Memory mapping does not help the header.** `var_fpos` is computed (a prefix sum) from 1–4-byte record
  lengths on disk, so the on-disk bytes are not the in-memory array. pgenlib removed its mmap mode on
  2022-03-14 (`pgenlib_read.h:432-437`; it points to commit `c470317` for "a webserver backend that can
  address multiple queries in parallel"). Remote files cannot be mapped.

## 4. The tiers

| Tier | What | Where | Bound |
|---|---|---|---|
| 0 | Catalog grid, companion paths, parquet footers | ATTACH cache (exists); DuckDB's parquet metadata cache | small; no limit needed |
| 1 | **Decoded `.pgen` headers, one per cell file, shared by all threads and queries** | new, process-wide `ObjectCache` entries | byte limit, least-recently-used out |
| 2 | Sample metadata per sample block (`.psam`) | ATTACH cache, changed: shared, not copied; compact IIDs | byte limit; over the limit → not cached |
| 3 | Variant metadata | never cached | served by the ID index and region loads (#140, #141) |
| 4 (later) | Headers in their compact on-disk form, decoded per 64k-variant block | pgenlib patch | only if tier 1's limit proves too small |

### Tier 1: the shared header cache (the first PR)

- **Key:** the `.pgen` path plus its identity: size and modification time (local), or ETag /
  Last-Modified (remote, through the VFS). **A stale header gives wrong records, not lost rows.** It is
  worse than a stale ID index, so the identity check is required.
- **Value:** a `shared_ptr` to the phase-2 `PgenFileInfo` and its allocation, with `shared_ff` closed.
  `GetEstimatedCacheMemory` returns the allocation size. DuckDB's `ObjectCache` is an LRU cache that
  counts each entry against the buffer pool (`object_cache.hpp:60-140`). An evicted entry stays alive
  while a reader holds it.
- **Build once:** the first query that needs a cell decodes its header under a per-key lock. Others wait,
  then share it.
- **Readers:** each thread still runs `PgrInit(fname, …)` for its own `PgenReader`, through
  `plinking_pgen_fopen` (local, remote and `localize` work as today).
- **Setting:** `plinking_pgen_header_cache` (BOOLEAN, default true). The byte limit is DuckDB's own
  object-cache limit, which counts against `memory_limit`.
- **Scope of the first PR:** `read_pfile` (`PfileLocalState`) and `ConcatPgenReader`. The other readers
  (`plink_freq`, `plink_hardy`, `plink_missing`, `plink_ld`, `plink_score`, `plink_glm`, `plink_pca`,
  `read_pgen`) follow in a second PR; they have the same 10 call sites.
- **Expected gain** (from the measurements above): header memory no longer grows with threads; the
  decode (52 ms per 10M variants; ~5 ms per 1M-variant cell) and the remote header fetch happen once
  per process, not once per thread per query.

### Tier 2: sample metadata

- Today a `Get` copies the whole entry out, at ~103 B per sample. 8M samples, IID only, would be
  ~820 MB per copy (inferred, linear).
- **Change:** hold the entry as `shared_ptr<const …>`, so a query reads it in place; store IIDs in one
  character buffer plus offsets (~14 B per sample for a 10-character IID; inferred).
- **Per sample block:** with ancestry groups, each group's `.psam` is its own entry, and a `WHERE IID` or
  an ancestry filter skips whole groups.
- **Limit:** an entry over `plinking_attach_cache_max_samples_bytes` (a proposed setting) is not cached;
  that query reads the `.psam.parquet` with only the projected columns (v0.6.1 behaviour).

## 5. Rejected

- **Cache all metadata on ATTACH:** ~2 GB per 10M variants for the variant axis alone (measured).
- **mmap the headers:** see §3.
- **One process-wide header per variant chunk:** wrong, since each ancestry group's header differs.

## 6. Decisions (Teague, 2026-11-20)

1. `ATTACH … (CACHE false)` does **not** turn off the tier-1 header cache. Only
   `plinking_pgen_header_cache` does; the header cache is process-wide.
2. The bound is DuckDB's object-cache limit (8 GiB by default, counted against `memory_limit`).
   No separate byte limit.
3. Chunks vary, but are likely **10–25M variants** each, with the samples split by ancestry. One
   cell's header is then 90–225 MB (9 B per variant, arithmetic); today 16 threads would hold
   1.4–3.6 GB for one cell.
4. The cost of caching is accepted: a single-threaded full scan now keeps every cell's header
   (measured on the 50-cell fixture: 1521 → 1949 MB peak).

## 7. Delivery

1. Tier 1 in `read_pfile` and `ConcatPgenReader`, with RED tests: a counter of header decodes (as
   `plinking_debug_genotype_where_pushed`) proves one decode per cell across threads and across queries;
   a rewritten `.pgen` (new mtime/size) is decoded again.
2. Tier 1 in the other readers.
3. Tier 2.
4. Tier 4 only if measurements at the real chunk size call for it.
