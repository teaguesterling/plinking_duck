# Development

## Building from Source

### Prerequisites

- C++17 compiler (GCC or Clang)
- CMake 3.12+
- Make
- Git (with submodule support)
- **Eigen3** — required (header-only); `plink_pca` needs it, and it is supplied
  by vcpkg. See [Build dependency: Eigen3](#build-dependency-eigen3).

### Build

```sh
git clone --recurse-submodules https://github.com/teaguesterling/plinking_duck.git
cd plinking_duck
make -j4
```

Build outputs:

| Path | Description |
|------|-------------|
| `build/release/duckdb` | DuckDB shell with extension auto-loaded |
| `build/release/test/unittest` | Test runner |
| `build/release/extension/plinking_duck/plinking_duck.duckdb_extension` | Loadable extension binary |

### Build dependency: Eigen3

`plink_pca` is implemented with **Eigen3** (header-only), which is a **required**
build dependency. It is listed in `vcpkg.json`, so the released binaries and any
build through the vcpkg toolchain always include `plink_pca`.

```sh
sudo apt install libeigen3-dev   # Debian/Ubuntu
brew install eigen               # macOS
vcpkg install eigen3             # vcpkg
```

A plain local `make` that does **not** use the vcpkg toolchain, on a machine
without a system Eigen3, fails at **configure time** with a clear message
(issue #53) rather than silently omitting `plink_pca`:

```
CMake Error at /path/to/plinking_duck/CMakeLists.txt:223 (message):
  Eigen3 not found. It is a required build dependency of plinking_duck
  (needed by plink_pca) and is declared in vcpkg.json ...
```

Install Eigen3 (above) or build through vcpkg to resolve it. CMake looks for a
vcpkg/CONFIG package first, then falls back to a system install under
`/usr/include/eigen3` or `/usr/local/include/eigen3`.

### Clean Build

```sh
make clean && make -j4
```

## Running Tests

```sh
make test
```

Tests use DuckDB's [sqllogictest](https://duckdb.org/docs/dev/sqllogictest/intro.html) framework. Test files are in `test/sql/` and test data in `test/data/`.

!!! note "`plink_pca` tests always run"

    Because Eigen3 is a hard build dependency (issue #53), a build either
    configures with `plink_pca` present or fails at configure time — there is no
    successfully-built suite missing `plink_pca`. So `test/sql/plink_pca.test`
    and `test/sql/plink_pca_negative.test` always run and are **not** guarded
    with `require`: a skip would turn the suite green precisely when `plink_pca`
    stopped being built, hiding a real regression.

### Test Data Files

| File | Description |
|------|-------------|
| `example.*` | .pvar/.psam/.bim/.fam test data for text file readers |
| `pgen_example.*` | 4 variants x 4 samples (main pgen test dataset) |
| `pgen_example.bim` | .bim companion for pgen_example (4 variants) |
| `pgen_orphan.*` | .pgen + .pvar only, no .psam (index-only mode testing) |
| `all_missing.*` | 2 variants x 2 samples, all genotypes missing |
| `large_example.*` | 3000 variants x 8 samples, 3 chroms x 1000 each (multi-batch/parallel tests) |

Test data can be regenerated with `test/data/generate_test_data.sh` (requires `plink2` binary).

## Project Structure

```
src/
├── plinking_duck_extension.cpp     # Entry point, registers all functions
├── pvar_reader.cpp / .hpp          # read_pvar()
├── psam_reader.cpp / .hpp          # read_psam()
├── pgen_reader.cpp / .hpp          # read_pgen()
├── pfile_reader.cpp / .hpp         # read_pfile()
├── plink_common.cpp / .hpp         # Shared infrastructure
├── plink_freq.cpp / .hpp           # plink_freq()
├── plink_hardy.cpp / .hpp          # plink_hardy()
├── plink_missing.cpp / .hpp        # plink_missing()
├── plink_ld.cpp / .hpp             # plink_ld()
└── plink_score.cpp / .hpp          # plink_score()
test/
├── sql/                            # sqllogictest files
└── data/                           # Test fixtures
third_party/
└── plink-ng/                       # pgenlib submodule (read path only)
```

## Architecture

### Table Function Lifecycle

DuckDB table functions follow a four-phase lifecycle:

1. **Bind** -- parse parameters, load metadata, define output schema
2. **InitGlobal** -- set up shared state (variant range, projection flags)
3. **InitLocal** -- per-thread state (pgenlib readers, buffers)
4. **Scan** -- emit rows in batches

### Shared Infrastructure (`plink_common`)

The `plink_common` module provides reusable components:

- **AlignedBuffer** -- RAII wrapper for cache-aligned allocations
- **VariantMetadata** -- pre-loaded variant info (CHROM, POS, ID, REF, ALT)
- **SampleSubset** -- sample bitmask, interleaved vec, cumulative popcounts
- **VariantRange** -- parsed region filter
- **LoadVariantMetadata** -- single-read `.pvar`/`.bim` parser
- **ResolveSampleIndices** -- `LIST(INTEGER)` / `LIST(VARCHAR)` sample parameter dispatch
- **FindCompanionFile** -- `.pgen` → `.pvar`/`.psam` discovery

### pgenlib Integration

PlinkingDuck links against the pgenlib library from the [plink-ng](https://github.com/chrchang/plink-ng) project (read path only, `pgenlib_write.cc` is excluded).

Key pgenlib patterns:

- **Two-phase init**: `PgfiInitPhase1` → allocate → `PgfiInitPhase2` → `PgrInit`
- **Per-thread readers**: each DuckDB thread gets its own `PgenFileInfo` + `PgenReader` (FILE* is not thread-safe)
- **Cleanup order**: `CleanupPgr` before `CleanupPgfi`
- **Buffer sizing**: use `NypCtToAlignedWordCt` (not `DivUp`) for SIMD-safe genovec allocation

### Dependencies

| Dependency | Source | Notes |
|------------|--------|-------|
| zstd | pgenlib vendored copy | Hidden visibility (DuckDB bundles its own) |
| libdeflate | plink-ng submodule | DuckDB doesn't include it |
| zlib | System | |
| simde | plink-ng submodule | Header-only SIMD emulation |
| Eigen3 | vcpkg, or system `libeigen3-dev` | **Required** (header-only) by `plink_pca`; a build that cannot find it fails at configure time (see [above](#build-dependency-eigen3)) |

## CI / Platform Support

CI runs via GitHub Actions using DuckDB's extension distribution pipeline.

### Supported Platforms

- **Linux amd64**: fully supported
- **Linux arm64**: fully supported

### Excluded Platforms

| Platform | Reason |
|----------|--------|
| `linux_amd64_musl` | pgenlib SIMD segfaults at runtime (compiles with `rawmemchr` shim) |
| `osx_amd64`, `osx_arm64` | Linker symbol visibility issue |
| `windows_amd64`, `windows_amd64_mingw` | MSVC can't compile plink-ng libdeflate (`__attribute__` syntax) |
| `wasm_*` | plink-ng requires POSIX APIs (pthreads, file I/O) |

Platform exclusions are configured in `.github/workflows/MainDistributionPipeline.yml`.

## Conventions

- **sqllogictest format specifiers**: `T` = varchar, `I` = integer, `R` = real
- **Named parameters**: use DuckDB's `:=` syntax (e.g., `pvar := 'path.pvar'`)
- **Type-dispatched parameters**: `samples` and `variants` use `LogicalType::ANY` and dispatch on actual value type in the bind function
