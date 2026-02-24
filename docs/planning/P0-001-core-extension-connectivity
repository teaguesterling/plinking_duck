# plinking_duck Implementation Plan

## Overview

A DuckDB C++ extension that reads PLINK2 files (.pgen/.pvar/.psam) and exposes
plink utility commands as SQL functions.

**Table functions**: `read_pvar`, `read_psam`, `read_pgen`, `read_pfile`
**Utility functions** (future): `plink_freq`, `plink_hardy`, `plink_ld`, etc.

## Architecture Decisions

### C++ Extension API (not C)
- Template is C++, pgenlib is C++, DuckDB C++ API is more ergonomic
- Use `TableFunction` class with typed bind/init/scan callbacks
- Subclass `TableFunctionData`, `GlobalTableFunctionState`, `LocalTableFunctionState`
- Register via `ExtensionLoader::RegisterFunction()`

### Threading: DuckDB only, no pgenlib pthreads
- pgenlib's `plink2_thread.h/.cc` provides its own threading primitives (ThreadGroup, etc.)
- We must NOT use these. DuckDB controls all threading via its parallel scan model.
- pgenlib's `PgenReader` is designed to be per-thread (share one `PgenFileInfo`,
  create one `PgenReader` per thread). This maps perfectly to DuckDB's
  `GlobalTableFunctionState` (shared) / `LocalTableFunctionState` (per-thread) pattern.
- `plink2_thread.cc` must still compile (other pgenlib code references symbols from it),
  but we never call `SpawnThreads()`, `JoinThreads()`, etc.

### Dependency strategy
- pgenlib is in submodule at `third_party/plink-ng/` (pinned at commit `4ce97faa`)
- **zstd**: DuckDB bundles zstd at `duckdb/third_party/zstd/`. Attempt to use DuckDB's
  zstd headers and link. If symbol conflicts arise, compile pgenlib's vendored zstd
  (`third_party/plink-ng/2.0/zstd/`) as a static lib with `-fvisibility=hidden`.
- **libdeflate**: pgenlib vendors this at `third_party/plink-ng/2.0/libdeflate/`.
  DuckDB does NOT bundle libdeflate. Compile it from the submodule.
- **zlib**: System library (`find_package(ZLIB)`). DuckDB has `miniz` but that's internal.
- **simde**: Header-only SIMD portability layer at `third_party/plink-ng/2.0/simde/`.
  Just add to include path. No linking.
- **pthreads**: System library. Needed for plink2_thread.cc compilation even though
  we don't use pgenlib's threading.

### Genotype output design
- **Default (variant-centric)**: One row per variant. Genotypes as `LIST(TINYINT)` where
  values are alt allele counts: 0=hom_ref, 1=het, 2=hom_alt, NULL=missing.
- **Tidy mode** (`tidy := true`, future): One row per (variant, sample). Joins variant
  metadata + sample metadata + genotype. Massive expansion but SQL-friendly.
- **Dosage** (`dosages := true`): Adds `dosages LIST(FLOAT)` column (0.0-2.0 scale).
- Projection pushdown is critical: skip genotype decoding if genotype columns aren't selected.

---

## Phase 0: Build Infrastructure

**Goal**: Get pgenlib compiling as part of the extension build.

### 0.1 Install cmake (prerequisite)
System doesn't have cmake installed. Need `cmake >= 3.5`.
```
sudo apt-get install cmake
# or: pip install cmake
# or: conda install cmake
```

### 0.2 Remove OpenSSL dependency
The template ships with OpenSSL as an example dependency. Remove it:
- `CMakeLists.txt`: Remove `find_package(OpenSSL REQUIRED)` and OpenSSL link lines
- `vcpkg.json`: Remove `"openssl"` from dependencies
- `src/plinking_duck_extension.cpp`: Remove `#include <openssl/opensslv.h>` and
  the openssl version function

### 0.3 CMakeLists.txt: Add pgenlib as a static library

```cmake
# Paths
set(PLINK_NG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/plink-ng/2.0)
set(PGENLIB_DIR ${PLINK_NG_DIR}/include)
set(LIBDEFLATE_DIR ${PLINK_NG_DIR}/libdeflate)
set(ZSTD_DIR ${PLINK_NG_DIR}/zstd)
set(SIMDE_DIR ${PLINK_NG_DIR}/simde)

# --- libdeflate (C library) ---
add_library(plink_libdeflate STATIC
    ${LIBDEFLATE_DIR}/lib/adler32.c
    ${LIBDEFLATE_DIR}/lib/crc32.c
    ${LIBDEFLATE_DIR}/lib/deflate_compress.c
    ${LIBDEFLATE_DIR}/lib/deflate_decompress.c
    ${LIBDEFLATE_DIR}/lib/gzip_compress.c
    ${LIBDEFLATE_DIR}/lib/gzip_decompress.c
    ${LIBDEFLATE_DIR}/lib/utils.c
    ${LIBDEFLATE_DIR}/lib/zlib_compress.c
    ${LIBDEFLATE_DIR}/lib/zlib_decompress.c
    ${LIBDEFLATE_DIR}/lib/arm/arm_cpu_features.c
    ${LIBDEFLATE_DIR}/lib/x86/x86_cpu_features.c
)
target_include_directories(plink_libdeflate PUBLIC ${LIBDEFLATE_DIR})
target_compile_definitions(plink_libdeflate PRIVATE LIBDEFLATE_STATIC)
set_target_properties(plink_libdeflate PROPERTIES
    C_VISIBILITY_PRESET hidden
    POSITION_INDEPENDENT_CODE ON
)

# --- zstd (C library, pgenlib's vendored copy) ---
# Compile ONLY decompression + common (we only read, never write zstd)
# Use hidden visibility to avoid symbol conflicts with DuckDB's bundled zstd
file(GLOB PLINK_ZSTD_COMMON_SRCS ${ZSTD_DIR}/lib/common/*.c)
file(GLOB PLINK_ZSTD_DECOMPRESS_SRCS ${ZSTD_DIR}/lib/decompress/*.c)
file(GLOB PLINK_ZSTD_COMPRESS_SRCS ${ZSTD_DIR}/lib/compress/*.c)
add_library(plink_zstd STATIC
    ${PLINK_ZSTD_COMMON_SRCS}
    ${PLINK_ZSTD_DECOMPRESS_SRCS}
    ${PLINK_ZSTD_COMPRESS_SRCS}
)
target_include_directories(plink_zstd PUBLIC
    ${ZSTD_DIR}/lib
    ${ZSTD_DIR}/lib/common
)
target_compile_definitions(plink_zstd PRIVATE ZSTD_DISABLE_ASM)
set_target_properties(plink_zstd PROPERTIES
    C_VISIBILITY_PRESET hidden
    POSITION_INDEPENDENT_CODE ON
)

# --- pgenlib (C++ library) ---
add_library(pgenlib STATIC
    ${PGENLIB_DIR}/plink2_base.cc
    ${PGENLIB_DIR}/plink2_bits.cc
    ${PGENLIB_DIR}/plink2_bgzf.cc
    ${PGENLIB_DIR}/plink2_htable.cc
    ${PGENLIB_DIR}/plink2_memory.cc
    ${PGENLIB_DIR}/plink2_string.cc
    ${PGENLIB_DIR}/plink2_text.cc
    ${PGENLIB_DIR}/plink2_thread.cc
    ${PGENLIB_DIR}/plink2_zstfile.cc
    ${PGENLIB_DIR}/pgenlib_misc.cc
    ${PGENLIB_DIR}/pgenlib_read.cc
    ${PGENLIB_DIR}/pgenlib_ffi_support.cc
    ${PGENLIB_DIR}/pvar_ffi_support.cc
)
target_include_directories(pgenlib PUBLIC
    ${PGENLIB_DIR}
    ${LIBDEFLATE_DIR}
    ${SIMDE_DIR}
)
target_link_libraries(pgenlib PRIVATE plink_libdeflate plink_zstd)
target_compile_definitions(pgenlib PRIVATE
    ZSTD_DISABLE_ASM
    LIBDEFLATE_STATIC
)
# Note: do NOT add -DNO_CPP11_TYPE_ENFORCEMENT; that's only for Cython
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
target_link_libraries(pgenlib PRIVATE ZLIB::ZLIB Threads::Threads)
set_target_properties(pgenlib PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    C_VISIBILITY_PRESET hidden
    POSITION_INDEPENDENT_CODE ON
    CXX_STANDARD 17
)

# --- Link to extension targets ---
target_link_libraries(${EXTENSION_NAME} pgenlib plink_libdeflate plink_zstd ZLIB::ZLIB Threads::Threads)
target_link_libraries(${LOADABLE_EXTENSION_NAME} pgenlib plink_libdeflate plink_zstd ZLIB::ZLIB Threads::Threads)
target_include_directories(${EXTENSION_NAME} PRIVATE ${PGENLIB_DIR})
target_include_directories(${LOADABLE_EXTENSION_NAME} PRIVATE ${PGENLIB_DIR})
```

### 0.4 Smoke test
Replace the template scalar functions with a minimal table function that just
verifies pgenlib links correctly:
```cpp
// In plinking_duck_extension.cpp LoadInternal():
#include "pgenlib_read.h"
// Verify we can call pgenlib functions at link time
// (actual table functions come in later phases)
```

### 0.5 Validation
```bash
make clean && make
# Should produce build/release/extension/plinking_duck/plinking_duck.duckdb_extension
```

---

## Phase 1: `read_pvar` Table Function

**Goal**: Parse .pvar (and .bim) files into a DuckDB table.

### File: `src/pvar_reader.cpp` (+ `src/include/pvar_reader.hpp`)

### Schema
The .pvar format has an optional `#CHROM` header line. If present, columns are named
from the header. If absent, assume .bim format.

#### .pvar (with header):
| Column | Type | Notes |
|--------|------|-------|
| CHROM | VARCHAR | Chromosome |
| POS | INTEGER | Base-pair position |
| ID | VARCHAR | Variant identifier |
| REF | VARCHAR | Reference allele |
| ALT | VARCHAR | Alternate allele(s); could be LIST(VARCHAR) if multi-allelic |
| QUAL | FLOAT | Optional |
| FILTER | VARCHAR | Optional |
| INFO | VARCHAR | Optional, raw semicolon-separated string |
| CM | DOUBLE | Optional, centimorgan position |

#### .bim (no header, 6 columns):
`CHROM  ID  CM  POS  ALT  REF` (note: different column order than .pvar)

### Implementation pattern

```cpp
// --- Bind Data ---
struct PvarBindData : public TableFunctionData {
    string file_path;
    vector<string> column_names;  // discovered from header
    vector<LogicalType> column_types;
    bool is_bim_format;           // no header → .bim column order
};

// --- Global State ---
struct PvarGlobalState : public GlobalTableFunctionState {
    mutex lock;
    idx_t current_line;           // for sequential scan
    bool finished;
};

// --- Local State ---
struct PvarLocalState : public LocalTableFunctionState {
    // File handle (each thread gets its own for parallel reads)
    // For Phase 1, single-threaded is fine.
    unique_ptr<FileHandle> file_handle; // or std::ifstream
    string line_buffer;
    idx_t current_line;
    bool finished;
};

// --- Bind ---
// 1. Open file, read first lines to detect format
// 2. If line starts with "##" → skip (pvar comment)
// 3. If line starts with "#CHROM" → parse header columns
// 4. If neither → .bim format, use default 6-column schema
// 5. Populate return_types and names from discovered schema
// 6. Store file_path and schema in PvarBindData

// --- Scan ---
// 1. Read lines from file
// 2. Split on tabs
// 3. Map fields to columns based on schema
// 4. Handle missing values (. → NULL)
// 5. Fill DataChunk vectors, set cardinality
```

### Named parameters
- None initially. Could add `columns` filter later.

### Registration
```cpp
TableFunction read_pvar("read_pvar", {LogicalType::VARCHAR}, PvarScan, PvarBind, PvarInitGlobal, PvarInitLocal);
read_pvar.projection_pushdown = true;
loader.RegisterFunction(read_pvar);
```

### Test: `test/sql/read_pvar.test`
Need test .pvar and .bim files. Create small synthetic ones in `test/data/`.

---

## Phase 2: `read_psam` Table Function

**Goal**: Parse .psam (and .fam) files into a DuckDB table.

### File: `src/psam_reader.cpp` (+ `src/include/psam_reader.hpp`)

### Schema

#### .psam (with header):
| Column | Type | Notes |
|--------|------|-------|
| FID | VARCHAR | Family ID (optional, may not exist) |
| IID | VARCHAR | Individual/sample ID (always present) |
| SID | VARCHAR | Source ID (optional) |
| PAT | VARCHAR | Paternal ID (optional) |
| MAT | VARCHAR | Maternal ID (optional) |
| SEX | INTEGER | 1=male, 2=female, 0/NA=unknown |
| *phenotypes* | varies | Any remaining columns are phenotypes/covariates |

#### .fam (no header, 6 columns):
`FID  IID  PAT  MAT  SEX  PHENO1`

### Implementation
Nearly identical pattern to `read_pvar`. Text parsing, header detection, dynamic schema.

### Test: `test/sql/read_psam.test`

---

## Phase 3: `read_pgen` Table Function

**Goal**: Read .pgen binary genotype files with variant metadata from .pvar and
sample info from .psam.

### File: `src/pgen_reader.cpp` (+ `src/include/pgen_reader.hpp`)

### pgenlib initialization (bind phase)

```cpp
// In PgenBindData:
// Phase 1: Open file, read header
plink2::PgenFileInfo pgfi;
plink2::PreinitPgfi(&pgfi);
plink2::PgenHeaderCtrl header_ctrl;
uintptr_t pgfi_alloc_cacheline_ct;
char errstr_buf[plink2::kPglErrstrBufBlen];

PglErr err = plink2::PgfiInitPhase1(
    pgen_path.c_str(),
    nullptr,           // .pgi path (auto-detect)
    UINT32_MAX,        // raw_variant_ct (read from file)
    UINT32_MAX,        // raw_sample_ct (read from file)
    &header_ctrl,
    &pgfi,
    &pgfi_alloc_cacheline_ct,
    errstr_buf
);
// pgfi.raw_variant_ct and pgfi.raw_sample_ct now available

// Phase 2: Load variant record types
unsigned char* pgfi_alloc;
plink2::cachealigned_malloc(pgfi_alloc_cacheline_ct * plink2::kCacheline, &pgfi_alloc);
uint32_t max_vrec_width;
uintptr_t pgr_alloc_cacheline_ct;
plink2::PgfiInitPhase2(header_ctrl, 0, 0, 0, 0,
    pgfi.raw_variant_ct, &max_vrec_width, &pgfi,
    pgfi_alloc, &pgr_alloc_cacheline_ct, errstr_buf);
```

### Schema
| Column | Type | Source | Notes |
|--------|------|--------|-------|
| CHROM | VARCHAR | .pvar | Chromosome |
| POS | INTEGER | .pvar | Position |
| ID | VARCHAR | .pvar | Variant ID |
| REF | VARCHAR | .pvar | Reference allele |
| ALT | VARCHAR | .pvar | Alternate allele(s) |
| genotypes | LIST(TINYINT) | .pgen | Alt allele counts [0,1,2,NULL] per sample |
| dosages | LIST(FLOAT) | .pgen | Optional, if dosage data present & requested |
| phased | LIST(BOOLEAN) | .pgen | Optional, if phase data present & requested |

### Bind phase
1. Accept .pgen path as first positional argument
2. Auto-discover .pvar and .psam by replacing the extension
   (or accept explicit paths as named parameters: `pvar := '...'`, `psam := '...'`)
3. Initialize PgenFileInfo (two-phase init as above)
4. Parse .pvar header for variant metadata columns
5. Parse .psam to get sample_ct and sample IDs
6. Store PgenFileInfo, variant metadata, sample metadata in bind data
7. Register columns: variant metadata + genotypes + optional dosages/phased

### Global state
```cpp
struct PgenGlobalState : public GlobalTableFunctionState {
    atomic<uint32_t> next_variant_idx{0};  // atomic counter for work-stealing
    uint32_t total_variants;

    idx_t MaxThreads() const override {
        // Partition variant ranges across threads
        // Each thread processes a contiguous block of variants
        return min(total_variants / 1000 + 1, (uint32_t)16);
    }
};
```

### Local state (per-thread)
```cpp
struct PgenLocalState : public LocalTableFunctionState {
    plink2::PgenReader pgr;              // per-thread reader
    unsigned char* pgr_alloc;            // reader's working memory
    uintptr_t* genovec;                  // 2-bit packed genotype buffer
    int32_t* genotype_buf;               // expanded int32 buffer for output
    // Optional:
    uintptr_t* dosage_present;
    uint16_t* dosage_main;
    // Variant range assigned to this thread
    uint32_t variant_start;
    uint32_t variant_end;
    uint32_t current_variant;
    bool finished;
};
```

### Scan function
```cpp
void PgenScan(ClientContext &context, TableFunctionInput &data, DataChunk &output) {
    auto &bind = data.bind_data->Cast<PgenBindData>();
    auto &global = data.global_state->Cast<PgenGlobalState>();
    auto &local = data.local_state->Cast<PgenLocalState>();

    idx_t row_count = 0;
    while (row_count < STANDARD_VECTOR_SIZE) {
        // Claim next variant (atomic)
        uint32_t vidx = global.next_variant_idx.fetch_add(1);
        if (vidx >= global.total_variants) {
            break;
        }

        // Fill variant metadata columns (from pre-loaded .pvar data)
        // Use projection pushdown: only fill columns DuckDB actually needs
        for (auto col_id : projected_columns) {
            if (col_id == COL_CHROM) { /* set chrom */ }
            else if (col_id == COL_POS) { /* set pos */ }
            // ...
            else if (col_id == COL_GENOTYPES) {
                // Read genotypes from .pgen
                PgrGet(nullptr, pssi, sample_ct, vidx, &local.pgr, local.genovec);
                // Convert 2-bit packed → int32 array (missing = -9)
                GenoarrToInt32sMinus9(local.genovec, sample_ct, local.genotype_buf);
                // Fill LIST(TINYINT) vector
                auto list_vec = output.data[col_output_idx];
                // ... populate list entries, converting -9 to NULL ...
            }
        }
        row_count++;
    }
    output.SetCardinality(row_count);
}
```

### Named parameters
| Parameter | Type | Default | Purpose |
|-----------|------|---------|---------|
| `pvar` | VARCHAR | auto | Explicit .pvar path |
| `psam` | VARCHAR | auto | Explicit .psam path |
| `dosages` | BOOLEAN | false | Include dosage column |
| `phased` | BOOLEAN | false | Include phase column |

### Projection pushdown (critical for performance)
If the query is `SELECT chrom, pos, id FROM read_pgen('file.pgen')`, we must NOT
decode genotypes at all. The `column_ids` from `TableFunctionInitInput` tells us
which columns DuckDB actually needs. Only call `PgrGet()` if the genotypes column
is in the projection list.

### Parallel scan strategy
- Variants are independent in .pgen (random access by variant index)
- Each thread claims batches of variants via atomic counter
- Each thread has its own PgenReader (required by pgenlib)
- Shared PgenFileInfo is immutable after init

### Test: `test/sql/read_pgen.test`
Need synthetic .pgen/.pvar/.psam test files. Use plink2 to create them from VCF,
or create programmatically.

---

## Phase 4: `read_pfile` Table Function

**Goal**: Unified reader combining all three files, with tidy mode support.

### File: `src/pfile_reader.cpp` (+ `src/include/pfile_reader.hpp`)

### Interface
```sql
-- Prefix-based (discovers .pgen, .pvar, .psam automatically)
SELECT * FROM read_pfile('data/example');

-- Explicit paths
SELECT * FROM read_pfile(pgen := 'a.pgen', pvar := 'b.pvar', psam := 'c.psam');

-- Tidy mode: one row per (variant, sample)
SELECT * FROM read_pfile('data/example', tidy := true);

-- With sample metadata joined in
SELECT chrom, pos, iid, sex, genotype
FROM read_pfile('data/example', tidy := true)
WHERE sex = 1 AND chrom = '22';
```

### Schema (default, variant-centric)
Same as `read_pgen` but with additional sample metadata accessible.

### Schema (tidy mode)
| Column | Type | Source |
|--------|------|--------|
| CHROM | VARCHAR | .pvar |
| POS | INTEGER | .pvar |
| ID | VARCHAR | .pvar |
| REF | VARCHAR | .pvar |
| ALT | VARCHAR | .pvar |
| FID | VARCHAR | .psam |
| IID | VARCHAR | .psam |
| SEX | INTEGER | .psam |
| genotype | TINYINT | .pgen |
| dosage | FLOAT | .pgen (optional) |
| *phenotype columns* | varies | .psam |

### Tidy mode scan pattern
Follows the DuckHTS BCF tidy pattern:
```
For each variant:
    Read genotypes for all samples
    For each sample:
        Emit one row: (variant_metadata..., sample_metadata..., genotype)
        If output chunk is full, return (resume at this sample on next call)
```

This requires tracking both `current_variant` and `current_sample` in local state.

### Named parameters
| Parameter | Type | Default | Purpose |
|-----------|------|---------|---------|
| `pgen` | VARCHAR | auto | Explicit .pgen path |
| `pvar` | VARCHAR | auto | Explicit .pvar path |
| `psam` | VARCHAR | auto | Explicit .psam path |
| `tidy` | BOOLEAN | false | Tidy output (one row per variant x sample) |
| `dosages` | BOOLEAN | false | Include dosage data |
| `phased` | BOOLEAN | false | Include phase data |
| `samples` | LIST(VARCHAR) | all | Filter to specific sample IDs |
| `region` | VARCHAR | all | Filter to genomic region (chr:start-end) |

### Sample subsetting
If `samples` parameter is provided:
1. Look up sample indices from .psam
2. Create sample_include bitarray
3. Use `PgrSetSampleSubsetIndex()` for efficient pgenlib subsetting
4. Only decode requested samples (huge performance win for biobank-scale data)

---

## Phase 5: Plink Utility Commands

**Goal**: Expose key plink2 computations as SQL functions.

These can be table-generating functions or aggregate functions. Priority order
based on common usage:

### 5.1 `plink_freq(path)` — Allele frequencies
```sql
SELECT * FROM plink_freq('data/example');
-- Returns: CHROM, POS, ID, REF, ALT, ALT_FREQ, OBS_CT
```
Uses `PgrGetCounts()` which is fast (doesn't fully decompress genotypes).

### 5.2 `plink_hardy(path)` — Hardy-Weinberg equilibrium
```sql
SELECT * FROM plink_hardy('data/example');
-- Returns: CHROM, POS, ID, REF, ALT, HOM_REF_CT, HET_CT, HOM_ALT_CT, P_HWE
```

### 5.3 `plink_missing(path)` — Missingness rates
```sql
SELECT * FROM plink_missing('data/example');
-- Returns: variant missingness and sample missingness tables
```
Uses `PgrGetMissingness()` which is efficient.

### 5.4 `plink_ld(path, variant1, variant2)` — Linkage disequilibrium
```sql
SELECT * FROM plink_ld('data/example', 'rs123', 'rs456');
-- Returns: R2, D_PRIME, etc.
```

### 5.5 Future: `plink_pca`, `plink_glm`, `plink_score`
These are more complex and may benefit from integration with DuckDB's expression
evaluation rather than being standalone table functions.

---

## File Organization

```
src/
  include/
    plinking_duck_extension.hpp    # Extension class declaration
    pvar_reader.hpp                # read_pvar declarations
    psam_reader.hpp                # read_psam declarations
    pgen_reader.hpp                # read_pgen declarations
    pfile_reader.hpp               # read_pfile declarations
    plink_utils.hpp                # Shared utilities (path resolution, etc.)
  plinking_duck_extension.cpp      # Extension entry point, registers all functions
  pvar_reader.cpp                  # Phase 1
  psam_reader.cpp                  # Phase 2
  pgen_reader.cpp                  # Phase 3
  pfile_reader.cpp                 # Phase 4
  plink_freq.cpp                   # Phase 5.1
  plink_hardy.cpp                  # Phase 5.2
test/
  data/
    example.pgen                   # Synthetic test data
    example.pvar
    example.psam
    example.bim                    # .bim format test
    example.fam                    # .fam format test
  sql/
    read_pvar.test
    read_psam.test
    read_pgen.test
    read_pfile.test
third_party/
  plink-ng/                        # Git submodule (already added)
```

---

## Key pgenlib API Reference

### Initialization (two-phase)
```cpp
plink2::PgenFileInfo pgfi;
plink2::PgenReader pgr;
plink2::PreinitPgfi(&pgfi);
plink2::PreinitPgr(&pgr);

// Phase 1: header
PgfiInitPhase1(filename, nullptr, UINT32_MAX, UINT32_MAX,
               &header_ctrl, &pgfi, &pgfi_alloc_cacheline_ct, errstr_buf);

// Allocate
cachealigned_malloc(pgfi_alloc_cacheline_ct * kCacheline, &pgfi_alloc);

// Phase 2: variant records
PgfiInitPhase2(header_ctrl, 0, 0, 0, 0, pgfi.raw_variant_ct,
               &max_vrec_width, &pgfi, pgfi_alloc,
               &pgr_alloc_cacheline_ct, errstr_buf);

// Per-thread reader
cachealigned_malloc(pgr_alloc_cacheline_ct * kCacheline, &pgr_alloc);
PgrInit(filename, max_vrec_width, &pgfi, &pgr, pgr_alloc);
```

### Reading functions
| Function | Purpose |
|----------|---------|
| `PgrGet(sample_include, pssi, sample_ct, vidx, &pgr, genovec)` | Basic 2-bit genotypes |
| `PgrGetP(...)` | Genotypes + phase |
| `PgrGetD(...)` | Genotypes + dosage |
| `PgrGetDp(...)` | Genotypes + phase + dosage (full) |
| `PgrGetCounts(...)` | Just genotype counts (fast, no decompression) |
| `PgrGetMissingness(...)` | Just missingness bitarray |

### FFI conversion
| Function | Purpose |
|----------|---------|
| `GenoarrToInt32sMinus9(genovec, ct, buf)` | 2-bit → int32 array (missing = -9) |
| `GenoarrToBytesMinus9(genovec, ct, buf)` | 2-bit → int8 array (missing = -9) |
| `Dosage16ToDoublesMinus9(...)` | Dosage → double array |

### Genotype encoding (2-bit packed, "nyparray")
- `0b00` = homozygous reference (0 alt alleles)
- `0b01` = heterozygous (1 alt allele)
- `0b10` = homozygous alternate (2 alt alleles)
- `0b11` = missing

### Cleanup
```cpp
CleanupPgr(&pgr, &reterr);
CleanupPgfi(&pgfi, &reterr);
aligned_free(pgr_alloc);
aligned_free(pgfi_alloc);
```

---

## DuckDB C++ Table Function Pattern Reference

```cpp
// === Bind Data ===
struct MyBindData : public TableFunctionData {
    string file_path;
    // ... immutable data discovered at bind time
};

// === Global State ===
struct MyGlobalState : public GlobalTableFunctionState {
    atomic<idx_t> next_row{0};
    idx_t MaxThreads() const override { return N; }
};

// === Local State ===
struct MyLocalState : public LocalTableFunctionState {
    // Per-thread resources (file handles, buffers, etc.)
};

// === Bind ===
unique_ptr<FunctionData> MyBind(ClientContext &context,
                                TableFunctionBindInput &input,
                                vector<LogicalType> &return_types,
                                vector<string> &names) {
    auto result = make_uniq<MyBindData>();
    result->file_path = input.inputs[0].GetValue<string>();
    // ... discover schema ...
    return_types.push_back(LogicalType::VARCHAR);
    names.push_back("col_name");
    return std::move(result);
}

// === Init Global ===
unique_ptr<GlobalTableFunctionState> MyInitGlobal(ClientContext &context,
                                                   TableFunctionInitInput &input) {
    return make_uniq<MyGlobalState>();
}

// === Init Local ===
unique_ptr<LocalTableFunctionState> MyInitLocal(ExecutionContext &context,
                                                 TableFunctionInitInput &input,
                                                 GlobalTableFunctionState *global_state) {
    auto result = make_uniq<MyLocalState>();
    // ... open file handles, allocate buffers ...
    return std::move(result);
}

// === Scan ===
void MyScan(ClientContext &context, TableFunctionInput &data, DataChunk &output) {
    auto &bind = data.bind_data->Cast<MyBindData>();
    auto &global = data.global_state->Cast<MyGlobalState>();
    auto &local = data.local_state->Cast<MyLocalState>();

    idx_t count = 0;
    while (count < STANDARD_VECTOR_SIZE) {
        // ... read data, fill vectors ...
        count++;
    }
    output.SetCardinality(count);
}

// === Registration ===
TableFunction func("my_func", {LogicalType::VARCHAR}, MyScan, MyBind, MyInitGlobal, MyInitLocal);
func.projection_pushdown = true;
func.named_parameters["option"] = LogicalType::BOOLEAN;
loader.RegisterFunction(func);
```

---

## Build & Test Commands

```bash
# Build
make clean && make

# Run DuckDB shell with extension
./build/release/duckdb -unsigned -cmd "LOAD 'build/release/extension/plinking_duck/plinking_duck.duckdb_extension';"

# Run tests
make test

# Debug build
make debug
make test_debug
```

---

## Known Risks & Mitigations

1. **zstd symbol conflicts**: pgenlib vendors zstd, DuckDB bundles zstd. Hidden visibility
   on the pgenlib static lib should prevent conflicts. If not, we may need to rename
   symbols via objcopy or use a linker version script.

2. **plink2_thread.cc compilation**: We include it for symbol resolution but never call
   threading functions. If it pulls in unwanted deps, stub it out with empty implementations.

3. **Large genotype arrays**: For biobank-scale data (500K+ samples), a single
   `LIST(TINYINT)` column entry contains 500K elements. DuckDB handles large lists
   but memory pressure is a concern. The tidy mode is the escape hatch.

4. **pgenlib memory model**: Uses `cachealigned_malloc` (memalign/posix_memalign).
   This should work on all target platforms but is worth testing.

5. **cmake not installed on dev machine**: Need to install before building.
   `sudo apt-get install cmake` or use conda/pip.
