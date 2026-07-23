# File Handling

How `read_pfile` (and `read_pgen` / the `plink_*` analysis functions) find, resolve,
and open PLINK 2 filesets — from a single local prefix to a variant-sharded biobank
on S3.

## The fileset

A PLINK 2 fileset is three files sharing a prefix:

| file | contents | reader |
|---|---|---|
| `.pgen` | compressed genotypes (binary) | pgenlib |
| `.pvar` (or `.bim`) | variant metadata | text / parquet |
| `.psam` (or `.fam`) | sample metadata | text / parquet |

`read_pfile('data/cohort')` opens `data/cohort.pgen` + `.pvar` + `.psam`. Legacy
PLINK 1 `.bim`/`.fam` are auto-detected. Each path can be overridden explicitly:
`read_pfile('data/cohort', pvar := 'data/cohort.annotated.pvar')`.

## Companion files

For each companion, discovery prefers a **parquet** sidecar when present
(`plinking_use_parquet_companions`, default `true`): `cohort.pvar.parquet` /
`cohort.psam.parquet` are read via DuckDB's native parquet reader before falling
back to text (`.pvar`/`.bim`, `.psam`/`.fam`, optionally `.zst`). Parquet companions
are a major speedup at sample scale — a wide `.psam` reads only the projected
columns. You can also point a companion at any table/view/query result (see
[read_pvar](../functions/read_pvar.md) / [read_psam](../functions/read_psam.md)).

## Path resolution

Reader inputs are resolved through DuckDB's virtual filesystem before opening, so
the same conveniences as `read_csv` apply:

- **`file_search_path`** — relative prefixes resolve under a search directory:
  `SET file_search_path = '/data/cohort'; SELECT * FROM read_pfile('chr22');`
- **Globs** — `read_pfile('data/chr*.pgen')` expands to every matching shard.
- **Protocols** — a registered protocol (e.g. the scalarfs `pathmacro:` catalog)
  resolves to concrete paths: `read_pfile('pathmacro:cohort?gene=BRCA1')`.

A single glob or protocol URL can therefore fan out to a whole **sharded fileset**.

## Multi-file / sharded reads

The first argument may be a `LIST(VARCHAR)` of prefixes (or a glob/protocol that
expands to one), read as one logical table:

```sql
SELECT * FROM read_pfile(['data/chr1', 'data/chr2', 'data/chr3']);
```

This targets the biobank layout where a fileset is **sharded by variant** — many
`.pgen`/`.pvar` files that share one sample set. Variants row-concatenate in list
order; in `orient := 'sample'` each subject's genotype array spans every shard's
variants.

- **`combine_samples`** — how the shards' sample sets combine: `'implicit'`
  (default; trust that all shards have the same IIDs in the same order — no
  per-shard check, no extra I/O), `'identical'` (verify IIDs match across shards;
  a mismatch is a hard error). `'union'`/`'intersect'`/`'concatenate'` are reserved.
- **Shared `psam`** — `read_pfile([...], psam := 'cohort.psam')` uses one shared
  `.psam` for every shard (the per-shard-`.pgen`/`.pvar`, one-`.psam`-elsewhere
  layout). `pgen`/`pvar` overrides remain single-file only.
- `region :=` filters each shard to the range and returns the union.

See [read_pfile](../functions/read_pfile.md#multi-file-input) for the full contract.

## Remote / cloud reads (`s3://`, `https://`, …)

With the `httpfs` extension loaded, every reader and analysis function can read a
`.pgen` **directly from any DuckDB filesystem** — the `.pvar`/`.psam` companions too:

```sql
LOAD httpfs;
SELECT IID FROM read_pfile('s3://bucket/cohort', orient := 'sample',
    genotypes := 'counts', region := '17:43000000-43125000',
    include_genotypes := ['het','hom_alt']);   -- BRCA1 carriers, over S3
```

Credentials/secrets and settings come from the query context automatically (via
DuckDB Secrets). This is genuinely ahead of the `plink2` CLI, which is local-only.

### How `.pgen` bytes are read — `plinking_pgen_io`

`.pgen` I/O is the one part of pgenlib that historically bypassed the VFS (raw
`FILE*`). `plinking_pgen_io` governs it:

| value | behavior |
|---|---|
| `'auto'` *(default)* | remote/VFS paths read through the VFS; local paths use native `fopen` (zero overhead) |
| `'native'` | always native `fopen` — errors on a remote path |
| `'vfs'` | always through the VFS (even local) — for testing |
| `'localize'` | *reserved* — download remote to a temp then read locally (not yet implemented) |

Under `'auto'`, **local reads are byte-for-byte the classic path** (no overhead);
only remote/VFS paths take the VFS route.

### What gets fetched

A **targeted query fetches only the bytes it needs** — the variant offset index
plus the region's/variants' records — via HTTP range reads, not the whole file. A
read-ahead block cache collapses per-read over-fetch to ~1×, so carrier/region
queries over a large remote `.pgen` are efficient.

**Full scans** of a remote `.pgen` read more (each parallel thread reads the index
and its variant range independently). For full-scan-over-remote workloads, load the
community **`cache_httpfs`** extension (a block cache over `httpfs`) or reduce
threads. See [Optimizations](optimizations.md#remote-cloud-pgen-reads).

### Limitations

- Split-index filesets (separate `.pgi`) are not yet supported — embedded-index
  `.pgen` only (the `plink2 --make-pgen` default).
- Remote **writes** are not supported (there is no writer yet).
