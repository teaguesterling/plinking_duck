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
- `example.*`: pvar/psam/bim/fam test data for P1-001/P1-002 — do not overwrite
- `generate_test_data.sh`: regenerates pgen fixtures (needs plink2)
- `plink2`: Binary is in `/mnt/aux-data/teague/Dev/spack/opt/spack/linux-zen3/plink2-2.0.0-a.6.9-e2l3wx22vy6tkmdwzhaexlml4rs3okx3`

## CI / Platform Support

- **Linux amd64 + arm64**: fully supported, build + all tests pass
- **macOS amd64 + arm64**: fully supported since v0.1.1
  - Requires `DUCKDB_CPP_EXTENSION_ENTRY` macro for loadable extension linking (macOS linker uses `-exported_symbol` for the `_duckdb_cpp_init` entry point)
  - Hidden visibility on pgenlib/libdeflate/zstd static libs is fine on macOS (absorbed into final .dylib)
- **Windows amd64 (MinGW/RTools)**: fully supported since v0.1.2
  - No code changes needed — MinGW/RTools use GCC which handles `__attribute__` syntax
  - Only MSVC (`windows_amd64`) remains unsupported
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
- `read_pfile` orient parameter: `'variant'` (default, one row per variant), `'genotype'` (one row per variant×sample, replaces `tidy := true`), `'sample'` (one row per sample, genotypes across variants)
- `orient := 'sample'` pre-reads all genotypes into a matrix at bind time; genotypes array dimension = effective variant count (not sample count)
- `tidy := true` still works as an alias for `orient := 'genotype'`; specifying both `orient` and `tidy` is an error
- `read_pgen` accepts `orient` parameter but only `'variant'` is valid (no sample metadata for other modes)
- `phased := true` changes genotype element type from `TINYINT` to `ARRAY(TINYINT, 2)` representing `[allele1, allele2]` haplotype pairs (0=REF, 1=ALT); missing genotypes → NULL
- Phase direction: `[0, 1]` = REF|ALT (also the unphased het convention), `[1, 0]` = ALT|REF; hom ref = `[0, 0]`, hom alt = `[1, 1]`
- `phased` works across all orient modes and genotypes modes (array, list); incompatible with `dosages := true`
- `genotypes := 'columns'` outputs one scalar TINYINT column per sample (variant orient) or per variant (sample orient); column names from IIDs or variant IDs
- `genotypes := 'columns'` has a 256-column limit without explicit `samples :=` / `variants :=` / `region :=` filter
- `genotypes := 'columns'` + `orient := 'genotype'` is an error (genotype mode already produces scalar output)

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
