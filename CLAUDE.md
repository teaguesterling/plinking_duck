# PlinkingDuck

DuckDB extension for reading PLINK 2 genomics file formats in SQL.

## Project Structure

- `src/` — C++ extension source
  - `plinking_duck_extension.cpp` — entry point, registers all table functions
  - `pvar_reader.cpp` / `.hpp` — `read_pvar()`: .pvar and .bim files
  - `psam_reader.cpp` / `.hpp` — `read_psam()`: .psam and .fam files
  - `pgen_reader.cpp` / `.hpp` — `read_pgen()`: .pgen binary genotype files (uses pgenlib), ARRAY(TINYINT, N) output
  - `pfile_reader.cpp` / `.hpp` — `read_pfile()`: unified reader with orient parameter (variant/genotype/sample modes)
  - `plink_common.cpp` / `.hpp` — shared P2 infrastructure: RAII wrappers, file utilities, sample subsetting, region filtering, OrientMode/GenotypeMode enums
  - `plink_freq.cpp` / `.hpp` — `plink_freq()`: per-variant allele frequencies via PgrGetCounts
  - `plink_hardy.cpp` / `.hpp` — `plink_hardy()`: per-variant HWE exact test p-values via PgrGetCounts
  - `plink_missing.cpp` / `.hpp` — `plink_missing()`: per-variant and per-sample missingness via PgrGetMissingness
  - `plink_ld.cpp` / `.hpp` — `plink_ld()`: pairwise linkage disequilibrium (r², D, D') via PgrGet/PgrGetP
  - `plink_score.cpp` / `.hpp` — `plink_score()`: polygenic risk scoring via PgrGetD
  - `plink_glm.cpp` / `.hpp` — `plink_glm()`: per-variant GWAS regression (linear, logistic, Firth) via plink2's regression engine
  - `plink2_glm_logistic_math.cpp` / `.hpp` — extracted plink2 logistic/Firth regression math functions
- `test/sql/` — DuckDB sqllogictest files (positive + negative per reader)
- `test/data/` — test fixtures (small VCF-derived PLINK files)
- `docs/planning/` — design docs and implementation plans
- `third_party/plink-ng/` — pgenlib submodule (read path only)

## Build & Test

- `make -j 4` to build (blq command: `build`)
- `make test` to run all sqllogictests (blq command: `test`)
- Build outputs: `./build/release/duckdb`, `./build/release/test/unittest`
- pgenlib_write.cc is intentionally excluded from the build

## Architecture Patterns

- **Table function lifecycle**: Bind → InitGlobal → InitLocal (per-thread) → Scan
- **Projection pushdown**: all readers support it; skip expensive work (e.g., PgrGet) when columns aren't referenced
- **VFS integration**: use DuckDB FileSystem for all file I/O (not raw stdio)
- **Bulk file reads**: read entire file into memory then parse (see ReadFileLines in plink_common)
- **Parallel scan**: atomic batch claiming on variant index; per-thread local state
- **Shared P2 infrastructure** (`plink_common.hpp/cpp`): AlignedBuffer RAII wrapper, VariantMetadata, SampleSubset (sample_include + interleaved_vec + cumulative_popcounts), VariantRange region filtering, companion file discovery. All P2 utility functions build on these.

## pgenlib Notes

- PgenFileInfo must outlive PgenReader (pgr holds reference to pgfi)
- Per-thread: each thread needs its own PgenFileInfo + PgenReader (FILE* not thread-safe)
- Buffer sizing: use `NypCtToAlignedWordCt` (not `DivUp`) for SIMD-safe genovec allocation
- Cleanup order: CleanupPgr before CleanupPgfi
- `GenoarrToBytesMinus9`: unpacks 2-bit genovec → int8 array (-9 = missing)
- `PgrGetCounts()`: fast genotype counting without decompression; requires `sample_include_interleaved_vec` (built via `FillInterleavedMaskVec`)
- `FillInterleavedMaskVec`: builds interleaved bit vector from sample_include mask; buffer must be `BitCtToAlignedWordCt`-aligned (not `DivUp`)
- `PgrSampleSubsetIndex`: init from sample_include + cumulative_popcounts for subsetting with PgrGetCounts
- `PgrGetMissingness()`: fast missingness extraction without decompression; returns bitarray (1 = missing), needs `sample_include` + `pssi` (not `interleaved_vec`), requires caller-provided `genovec_buf` scratch buffer
- `PopcountWords()`: count set bits in bitarray; `BitCtToWordCt(sample_ct)` for word_ct arg
- `ctzw()`: count trailing zeros — used with `word &= word - 1` pattern to iterate set bits in missingness bitarray

## Test Data

- `pgen_example.*`: 4 variants × 4 samples (main pgen test dataset)
- `pgen_example.bim`: .bim companion for pgen_example (4 variants)
- `pgen_orphan.*`: .pgen + .pvar only, no .psam (index-only mode)
- `all_missing.*`: 2 variants × 2 samples, all genotypes missing
- `large_example.*`: 3000 variants × 8 samples, 3 chroms × 1000 each (multi-batch/parallel tests)
- `streaming_example.*`: 50000 variants × 8 samples, 3 chroms (20K/15K/15K), biallelic A/G (streaming/threading tests)
- `wes_chr10.*`: 30000 variants × 10000 samples, chr10 WES-scale dataset (~11MB .pgen, parallel correctness tests)
- `example.*`: pvar/psam/bim/fam test data for P1-001/P1-002 — do not overwrite
- `generate_test_data.sh`: regenerates pgen fixtures (needs plink2)
- `plink2`: `/srv/physical/longbottom-nvme2/teague/Dev/spack/opt/spack/linux-haswell/plink2-2.0.0-a.6.9-7p3r57u7ynsbctzz7zqspoavti4vllrp/bin/plink2` (spack; the old `/mnt/aux-data/...` path is gone). Or `spack location -i plink2`. NOTE: plink2's CLI only writes **embedded**-index `.pgen`; split-index (`.pgen.pgi`) fixtures need plink-ng's `pgen_compress -i` (build via `make pgen_compress` in the submodule, or compile `pgen_compress.cc`+`pgenlib_write.cc` against the extension's prebuilt pgenlib objects).

## CI / Platform Support

- **Linux amd64 + arm64**: fully supported, build + all tests pass
- **macOS amd64 + arm64**: fully supported since v0.1.1
  - Requires `DUCKDB_CPP_EXTENSION_ENTRY` macro for loadable extension linking (macOS linker uses `-exported_symbol` for the `_duckdb_cpp_init` entry point)
  - Hidden visibility on pgenlib/libdeflate/zstd static libs is fine on macOS (absorbed into final .dylib)
- **Windows amd64 (MinGW/RTools)**: disabled in CI pending further investigation
  - MinGW/RTools use GCC which handles `__attribute__` syntax
  - MSVC (`windows_amd64`) also unsupported (can't compile libdeflate)
- **Excluded platforms** (pre-existing plink-ng portability issues):
  - `linux_amd64_musl`: compiles with rawmemchr shim (`src/musl_compat.h`), but pgenlib SIMD segfaults at runtime
  - `windows_amd64` (MSVC): can't compile plink-ng libdeflate (`__attribute__` syntax)
  - `wasm_*`: plink-ng requires POSIX APIs (pthreads, file I/O)
- Exclusions configured in `.github/workflows/MainDistributionPipeline.yml` via `exclude_archs`
- `src/musl_compat.h`: inline rawmemchr for musl; force-included via CMake `check_symbol_exists`

## Conventions

- sqllogictest format specifiers: `T` = varchar, `I` = integer, `R` = real
- DuckDB function overloads: can't register two with same positional args; use `LogicalType::ANY` for type-dispatched named params
- Named parameters go in `read_pgen.named_parameters["key"]`
- `read_pfile`/`read_pvar` first positional arg accepts a single prefix `VARCHAR` **or** a `LIST(VARCHAR)` of prefixes/paths (mirrors `read_csv`/`read_json`); a list row-concatenates variants in file order. All files in a list must share the **same sample set (same IIDs, same order)** — the variant-sharded, identical-`.psam` layout; alignment is **positional** (trust-the-caller), but a mismatched `sample_ct` across files is a hard error. All three orients row-concat variants across files: `variant` + `genotype` at scan time, `sample` by pre-reading each shard into a matrix concatenated along the variant axis (genotypes array dimension = total variants across all shards). Named `pgen`/`pvar` overrides are single-file only (per-shard ambiguous — error if combined with a multi-element list); a `psam` override **is** allowed multi-file and applies as the one shared sample metadata for every shard ("use this psam for all pgens" — the per-shard-pgen+pvar, shared-psam biobank layout; sample_ct must still match every shard's `.pgen`).
- `combine_samples := '<mode>'` governs how sources' sample sets combine (multi-file): `implicit` (default — trust positional alignment, no per-shard IID check, no extra I/O), `identical` (verify every shard's IIDs are equal in the same order — one `.psam` parse per shard; a mismatch is a hard error), `concatenate` (P13-008, SAMPLE-sharded: same variants, disjoint samples; output samples = shards' samples in list order). `union`/`intersect` (real sample-set joins) are reserved and error as "not yet implemented".
- `concatenate` design: `sources` holds ONE PfileSource (shared variant axis — shard 0's `.pvar` or `pvar :=`, which is allowed only in this mode; shards 1..K's `.pvar` never read), `sample_shards` holds per-shard `.pgen`/sample_ct/offset, and `raw_sample_ct` = Σ (combined space). Every per-variant read site dispatches to `ConcatPgenReader` (`src/pgen_concat_reader.*`: per-shard pgr + split subset, bit-appends genovec/dosage/phase at the shard offset, SUMS genocounts) — `lstate.concat` / `pr_cat` / `cf_cat`. Emit code is shared, so all orients/modes/filters work. Rejects: `psam :=`/`pgen :=`, variant-count mismatch, differing `.psam` columns (genotype/sample orient), mixed text+parquet `.psam` in genotype/sample orient. Parquet companions (P13-009): `PfileSampleMetadata.source_paths` + `PsamSelectSql` = UNION ALL of per-shard `read_parquet(…, file_row_number)` tagged with the shard index, outer ORDER BY (shard, row) — used by `BuildProjectedPsamCdc` and `EnsureSourceIids`. Sparse streaming (difflist) is off in concat. `'concatenate'` trusts IID uniqueness; `combine_samples := 'disjoint'` (P13-009) = concatenate + a cross-shard IID uniqueness check (forces the IID load, errors naming both shards). Fixtures: `*_cat_*` = contiguous sample splits (oracle = unsplit fileset), `test/sql/read_pfile_sample_concat.test`. ATTACH lifts the sample-sharded reject via this mode. `read_pgen` multi-file input is planned (not yet available). Pairs with the scalarfs `pathmacro:` protocol (a catalog macro returns a `VARCHAR[]` of shard prefixes that feeds `read_pfile([...])`).
- `COPY … (FORMAT pfile, shard_samples {label: ['iid', …], …})` (P13-007, `src/pfile_writer.cpp`) is the write-side mate: one shared `out.pvar`, `out/<label>/data.{pgen,psam}`, `out.pfile.json` (`concat_dim:"samples"`, `shared.pvar`, per-file `psam`, sample-ordinal partitions). Groups must tile the input IID order as contiguous runs IN ORDER (strict). N `.pgen` writers open at once; files open lazily on the first row; cleanup removes all staging + created dirs (non-recursive SINGLE mode). The `.pgen` row path is split into `PgenRowGenotypes` (row guards, once per row) + `PgenAppendCallRange` (pack calls[lo,hi)). Test: `test/sql/write_pfile_sample_shard.test`.
- `read_pfile(catalog := '<x>.pfile.json')` reads a **v1 or v2 catalog** as one 2-D cohort (no positional arg; a zero-positional overload). v1 maps onto the grid (`LoadV1CatalogGrid`): concat_dim variants = N×1, samples = 1×K, single = 1×1; only the concat axis has declared extents (`CatalogGrid::has_*_extent`). **Discovery:** a positional VARCHAR prefix `p` reads `p.pfile.json` when it exists — the catalog WINS over a `p.pgen` beside it (Teague's rule, = ATTACH order); a positional `*.pfile.json` is a catalog; `catalog := false` disables, `catalog := true` requires (param type ANY). Local paths only (no `://`, `pathmacro:`, glob); a list never discovers. pgen/pvar/psam overrides with any catalog error and name `catalog := false`. Cells sorted row-major (variant block, then sample block); bind maps cell i to variant block `i / K`, sample block `i % K`. Concat mode generalizes to a grid: `sources` = one `PfileSource` per **variant block**, each with its own `shard_pgen_paths` (parallel to the bind-global `sample_shards` split); every `ConcatPgenReader::Open` uses `ShardsFor(src_idx)` + that source's `raw_variant_ct`. K = 1 is the plain variant-sharded path with per-cell template `.pvar`/`.psam`. Holes are a bind error (2-D meaning undefined). Each cell's `.pgen` counts must equal its blocks' declared extents. K > 1 defaults `combine_samples` to `'concatenate'`; `'implicit'`/`'identical'` error. ATTACH's matrix `c.genotypes_by_sample` = `read_pfile(catalog := …, orient := 'sample')` (replaced the v0.16.0 `range(error())` guard). ATTACH's `c.genotypes` = `read_pfile(catalog := <abs>)` for EVERY catalog (v1 and matrix; one row per variant, so a WHERE prunes blocks — Teague approved 2026-11-09, breaking for matrices), EXCEPT a holed matrix, which keeps the per-cell fan-out (`MatrixGenotypesSql`); a bare fileset keeps `GenotypesSql`.
- **Catalog pruning (P13-006 §10.4 stage 2):** `PruneCatalogGrid` drops cells right after `LoadCatalogGrid`, from the EXPLICIT params only (a `WHERE` does not prune yet). `region :=` drops a variant block whose declared range (`CatalogGrid::variant_ranges`: v1 file `range`, or a matrix block's CHROM enum value + POS-interval bounds) proves no overlap via `ShardEntry::RangeOverlaps` (verbatim contig, conservative keep); if none survive, block 0 is kept so the result equals a full read. INTEGER `samples :=` drops sample blocks holding none of the cohort-wide indices and REWRITES `sel.samples` into the surviving blocks; IIDs, invalid indices, and `combine_samples := 'disjoint'` prune nothing. Shape guards (variants/variant_range across blocks, the combine_samples contradiction) use `declared_*_block_ct`, never the pruned counts. Do not reuse `PruneCells` here: it assumes a CHROM key (read_pfile_catalog rejects `region` without one), but read_pfile must read such a matrix and keeps every block. Proof test: `read_pfile_catalog_prune.test` corrupts a cell that pruning must skip.
- **Lazy `.pvar` (P13-006 §10.4 stage 3):** `PfileBindData::lazy_variants` = catalog input AND orient variant/genotype AND no `variants :=`. Bind then loads each source with `need_pvar=false` (pgen header only), stores the catalog's `.pvar` path, and skips region/variant_range. `EnsureLazyVariants` (first call in `PfileInitGlobal`, once per bind data under `lazy_variants_lock`, via `const_cast`: InitGlobal's bind data is const but was created mutable) runs `LoadSourceVariants` + `BuildEffectiveVariantList` + `ApplyVariantRange` per source and rebuilds `variant_offsets` with `assign`. Anything that reads `sources[*].variants` / effective lists / `variant_offsets` must run after it. Sample orient stays eager (schema = variant count). The `.psam` files are still read at bind.
- **WHERE cell pruning (P13-006 §10.4 stage 4):** `fn.pushdown_complex_filter = PfilePushdownComplexFilter`. Only for a `lazy_variants` read not yet loaded: `CollectCellPruneConstraints` walks the AND-conjunction of bound Expressions (v2.0: comparisons and BETWEEN are `BoundFunctionExpression`s, use the `BoundComparisonExpression` / `BoundBetweenExpression` helpers; IN is a `BoundOperatorExpression`, children via `GetChildren()`), maps column refs via `get.GetColumnIds()[binding.column_index].GetPrimaryIndex()` (CHROM_COL 0, POS_COL 1), folds constants with `TryEvaluateScalar`, and drops `sources` (+ the parallel `source_ranges`) whose declared range cannot match. It NEVER removes a filter (rows = full read), keeps source 0 if all drop, and ignores OR / casts on the column / anything unrecognized. Dropping block 0 is safe: `ShardsFor` takes each source's own `.pgen` paths.
- **Lazy `.pgen` header (P13-006 §10.5, stage 5):** `PfileBindData::lazy_pgen` = `lazy_variants` AND `grid.has_variant_extent` AND NOT `plinking_pgen_io := 'localize'` (`PgenIoLocalizes`; localize copies every `.pgen` at bind). A cell defers when its sample count is known without the header: `grid.has_sample_extent` (v2 matrix) or `!concat && si > 0` (v1 variant concat; cell 0 opens at bind and gives `raw_sample_ct`). Bind calls `LoadPfileSource(..., read_header=false, resolve_overrides=false)`: a catalog cell's `.pgen` = prefix + `.pgen` and its `.pvar`/`.psam` are taken as written, with NO existence check (§10.8; a v1 entry that omits `pvar` still discovers it), sets `raw_variant_ct` / the sample count from the cell's declared extents, and marks `PfileSource::pgen_header_deferred`. `EnsureLazyVariants` runs `CheckDeferredPgenHeaders` per surviving source BEFORE its `.pvar` load: `ReadPgenHeaderCounts` on each shard (`ShardsFor`) and the bind-time count messages. Measured on a 10 × 4 grid: `WHERE CHROM = '3'` opens 4 `.pgen` (was 40), `DESCRIBE` 0; local time unchanged (interleaved A/B). Test: `read_pfile_catalog_lazy_pgen.test`.
- **WHERE IID sample prune (P13-006 §10.7):** the same hook, for genotype orient over sample shards (`IsConcat()`, >1 shard): `CellPruneConstraints::iid_col` (= `geno_sample_col_start` + the psam IID column) collects `IID = c` / `IID IN (…)`; `PruneSamplesByIid` narrows `sample_indices` (sorted, every duplicate IID kept, ∩ `samples :=`, one sample kept if none match) and rebuilds `count_filter_subset` under `include_genotypes`. Skipped with `af_range`/`ac_range` (frequency over the sample set). `CheckDeferredPgenHeaders` skips a shard with no selected sample. Sample orient (§10.9): only counts/stats with no filter param (`smp_headers_at_init`, decided from `named_parameters` before the source loop): bind defers each catalog cell's `.pgen` header (`lazy_pgen`) but still loads the `.pvar` (`LoadSourceVariants` after the declared counts are set); `EnsureLazyVariants` checks the headers at InitGlobal. `CellPruneConstraints::chrom_col`/`pos_col` are `INVALID_INDEX` in sample orient (column 0 is the psam's IID). Per-element sample orient is NOT pruned (matrix pre-read at bind, before the hook). Test: `read_pfile_catalog_iid_prune.test`.
- **Writers never overwrite without `force`** (Teague, 2026-11-10): every COPY (`pvar`/`pgen`/`psam`/`pfile`, plain and sharded) errors at BIND when its target exists (`RefuseExistingOutput` in `pfile_writer.cpp`); `force` / `force true` allows it. Plain `pfile` checks its 3 effective files + `<p>.pfile.json`; sharded checks `<p>.{pfile.json,pgen,pvar,psam}` + the `<p>/` tree. With `force`, a plain `pfile` write removes `<p>.pfile.json` just before promotion (else discovery would read the stale catalog over the old shards with the new `.psam`). Tests that overwrite must pass `force true`.
- `read_pfile` orient parameter: `'variant'` (default, one row per variant), `'genotype'` (one row per variant×sample), `'sample'` (one row per sample, genotypes across variants)
- `orient := 'sample'` pre-reads all genotypes into a matrix at bind time; genotypes array dimension = effective variant count (not sample count)
- `read_pgen` accepts `orient` parameter but only `'variant'` is valid (no sample metadata for other modes)
- `phased := true` changes genotype element type from `TINYINT` to `ARRAY(TINYINT, 2)` representing `[allele1, allele2]` haplotype pairs (0=REF, 1=ALT); missing genotypes → NULL
- Phase direction: `[0, 1]` = REF|ALT (also the unphased het convention), `[1, 0]` = ALT|REF; hom ref = `[0, 0]`, hom alt = `[1, 1]`
- `phased` works across all orient modes and genotypes modes (array, list); incompatible with `dosages := true`
- `genotypes := 'columns'` outputs one scalar TINYINT column per sample (variant orient) or per variant (sample orient); column names from IIDs or variant IDs
- `genotypes := 'columns'` has no column count limit (DuckDB has no hard column limit)
- `genotypes := 'columns'` + `orient := 'genotype'` is an error (genotype mode already produces scalar output)
- `af_range := {min: 0.01, max: 0.5}` filters variants by allele frequency; uses PgrGetCounts (no decompression) for cheap pre-filtering
- `ac_range := {min: 2, max: 10}` filters variants by allele count; can combine with af_range
- Count filters work in all orient modes and in both `read_pfile` and `read_pgen`; variant/genotype orient filter at scan time, sample orient filters at bind time (before schema)
- AF denominator is `2 * non_missing_samples` (matches plink_freq), not `2 * total_samples`
- `STD_ARRAY_DECL(uint32_t, 4, genocounts)` for PgrGetCounts (pgenlib uses std::array, not C arrays)
- `include_genotypes := ['het', 'hom_alt']` filters by hardcall category — canonical genotype filter. Categories: `hom_ref` (0), `het` (1), `hom_alt` (2), `missing` (-9). Case-insensitive, whitespace-trimmed; unknown labels error. Arbitrary (incl. non-contiguous, e.g. `['hom_ref','hom_alt']`) subsets. `missing` is a first-class category — selecting it retains missing genotypes/rows.
- `genotype_range := {min: 1, max: 2}` is the numeric alias of `include_genotypes` (contiguous ranges only over {0,1,2}); optional `include_missing := true` field retains missing. Specifying both `include_genotypes` and `genotype_range` errors.
- Both compile to `GenotypeRangeFilter` (`allowed[3]` mask + `include_missing`); all engine sites use `AllowsCall(dosage)` / `CheckGenotypeRange`. `include_genotypes`/`genotype_range` incompatible with `dosages := true`.
- Filter semantics by orient (row grain differs): **variant** orient — keep variant if any sample matches (`any_pass`), per-element NULL out non-matching sample values; **genotype** orient — skip rows whose genotype doesn't match; **sample** orient — skip sample rows with no matching genotype (built as a `sample_keep` list at bind, scan iterates it so only surviving rows materialize), per-element NULL out for array/list/columns, and **counts/stats aggregate modes report TRUE counts** (out-of-range hardcalls are NOT nulled to -9, so they are never miscounted as `missing`).
- `GenotypeRangeResult` from `CheckGenotypeRange`: `any_pass` = skip variant entirely when false (true when a selected category or, with `include_missing`, missing is present); `all_pass` = skip per-element check when true.

## Extension Config Options

- `plinking_max_matrix_elements` (BIGINT, default 16G): max genotype matrix elements for `orient := 'sample'` pre-read. Controls memory guard for variants × samples product. Set via `SET plinking_max_matrix_elements = <value>`.
- `plinking_max_threads` (BIGINT, default 0): max threads for parallel scan operations. 0 = default (hardcoded cap of 16), >0 = cap at this value. Set via `SET plinking_max_threads = <value>`. Applied in all parallel table functions (`read_pfile`, `read_pgen`, `plink_freq`, `plink_hardy`, `plink_glm`, `plink_missing`, `plink_ld`, `plink_score`, `plink_pca`).
- `plinking_use_parquet_companions` (BOOLEAN, default true): auto-discover `.pvar.parquet` and `.psam.parquet` companion files. When true, parquet companions are preferred over text formats in the companion file search order. Set via `SET plinking_use_parquet_companions = false` to disable. Parquet companions are loaded via a separate Connection invoking DuckDB's built-in parquet reader.
- `plinking_pgen_io` (VARCHAR, default `'auto'`): how `.pgen` bytes are read. `auto` = remote/VFS paths (`s3://`, `http`, any other `scheme://` incl. `file://`, `pathmacro:`, glob, `file_search_path`) go through DuckDB's VFS (Path V — a `fopencookie`/`funopen` `FILE*` over a `FileHandle`, opened extension-side with the query context so httpfs secrets apply), local paths use native `fopen` (Path L, zero overhead). `native` = always `fopen` (errors on remote). `vfs` = always VFS (even local; for testing). `localize` = materialize a local temp copy of the `.pgen` at bind (via the context FS, so httpfs secrets apply), rewrite the path, then read natively (Path L) — best for remote **full scans** (downloads once, `1×`, vs. `auto`'s per-thread range re-fetch). `localize` ALWAYS copies, even a local source (so the copy path is testable without a network — see `test/sql/read_pfile_localize.test`); no size guard yet. Temp dir = `plinking_localize_dir` (else DuckDB `temporary_directory`, created if missing); per-query, unlinked when bind data dies (`PgenLocalizeGuard`, move-only, held per-source in each reader's bind data). Wiring: `LocalizePgenIfRequested(context, pgen_path, guard)` before `PgenIoUseVfs` in each reader's bind; `pfile_reader` localizes per shard in `LoadPfileSource`. The `.pgen` open is the only VFS-bypassing part of pgenlib (raw `FILE*`); a minimal re-appliable patch (`patches/plink-ng/`) redirects its 3 `fopen` sites to `plinking_pgen_fopen` (`src/plinking_pgen_vfs.*`), which the extension wires via `PgenVfsScope` (`src/pgen_vfs_opener.*`) around every reader/analysis open. A read-ahead block cache in the cookie collapses httpfs per-read over-fetch to ~1×; targeted queries fetch only the index + region records. Full scans over remote re-read more (per-thread index) — use `localize`, load `cache_httpfs`, or reduce threads. **Split-index (`.pgi`) is supported transparently:** pgenlib auto-derives `<pgen>.pgi` for external-index `.pgen` (header `0x2x`) and opens it via `plinking_pgen_fopen` (so local/remote/`vfs` all work); `localize` also copies the `.pgi` companion. Non-standard `.pgi` names (not `<pgen>.pgi`) aren't discovered. Fixture: `pgen_split.*` via `pgen_compress -i` (test/sql/read_pfile_pgi.test).
- `plinking_localize_dir` (VARCHAR, default `''`): directory for `plinking_pgen_io := 'localize'` temp copies. Empty = DuckDB's `temporary_directory`, else `.`. Created if absent (race-tolerant). Temps are per-query.
- `.pgen` VFS shim lives in the `pgenlib` static lib (`src/plinking_pgen_vfs.cpp`, no DuckDB dep — symbol resolves for every target linking pgenlib); the DuckDB-side opener/scope in `src/pgen_vfs_opener.cpp`. Vendored pgenlib patches under `patches/plink-ng/*.patch` are applied idempotently at CMake configure (`git apply --reverse --check` guard); the submodule stays pinned.
- Config options registered via `DBConfig::GetConfig(db).AddExtensionOption()` in `Load()`, read at bind time via `context.TryGetCurrentSetting()`

<!-- blq:agent-instructions -->
## blq - Build Log Query

Run builds and tests via blq MCP tools, not via Bash directly:
- `mcp__blq_mcp__commands` - list available commands
- `mcp__blq_mcp__run` - run a registered command (e.g., `run(command="test")`)
- `mcp__blq_mcp__register_command` - register new commands
- `mcp__blq_mcp__status` - check current build/test status
- `mcp__blq_mcp__errors` - view errors from runs
- `mcp__blq_mcp__info` - detailed run info (supports relative refs like `+1`, `latest`)
- `mcp__blq_mcp__output` - search/filter captured logs (grep, tail, head, lines)

Do NOT use shell pipes or redirects in commands (e.g., `pytest | tail -20`).
Instead: run the command, then use `output(run_id=N, tail=20)` to filter.
<!-- /blq:agent-instructions -->
