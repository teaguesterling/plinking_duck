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
- `plink2`: Binary is in `/mnt/aux-data/teague/Dev/spack/opt/spack/linux-zen3/plink2-2.0.0-a.6.9-e2l3wx22vy6tkmdwzhaexlml4rs3okx3`

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
- `combine_samples := '<mode>'` governs how sources' sample sets combine (multi-file): `implicit` (default — trust positional alignment, no per-shard IID check, no extra I/O), `identical` (verify every shard's IIDs are equal in the same order — one `.psam` parse per shard; a mismatch is a hard error). `union`/`intersect`/`concatenate` (real sample-set joins) are reserved and error as "not yet implemented". `read_pgen` multi-file input is planned (not yet available). Pairs with the scalarfs `pathmacro:` protocol (a catalog macro returns a `VARCHAR[]` of shard prefixes that feeds `read_pfile([...])`).
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
- `plinking_pgen_io` (VARCHAR, default `'auto'`): how `.pgen` bytes are read. `auto` = remote/VFS paths (`s3://`, `http`, `pathmacro:`, glob, `file_search_path`) go through DuckDB's VFS (Path V — a `fopencookie`/`funopen` `FILE*` over a `FileHandle`, opened extension-side with the query context so httpfs secrets apply), local paths use native `fopen` (Path L, zero overhead). `native` = always `fopen` (errors on remote). `vfs` = always VFS (even local; for testing). `localize` = reserved (download-to-temp, not yet implemented). The `.pgen` open is the only VFS-bypassing part of pgenlib (raw `FILE*`); a minimal re-appliable patch (`patches/plink-ng/`) redirects its 3 `fopen` sites to `plinking_pgen_fopen` (`src/plinking_pgen_vfs.*`), which the extension wires via `PgenVfsScope` (`src/pgen_vfs_opener.*`) around every reader/analysis open. A read-ahead block cache in the cookie collapses httpfs per-read over-fetch to ~1×; targeted queries fetch only the index + region records. Full scans over remote re-read more (per-thread index) — load `cache_httpfs` or reduce threads. `.pgi` split-index not yet supported.
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
