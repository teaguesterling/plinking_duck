# Getting Started

## Building from Source

PlinkingDuck is built as a DuckDB extension using the standard extension build system.

### Prerequisites

- C++17 compiler (GCC or Clang)
- CMake 3.12+
- Make
- Git (with submodule support)

### Build

```sh
git clone --recurse-submodules https://github.com/teaguesterling/plinking_duck.git
cd plinking_duck
make -j4
```

This produces:

| Output | Description |
|--------|-------------|
| `build/release/duckdb` | DuckDB shell with extension auto-loaded |
| `build/release/test/unittest` | Test runner |
| `build/release/extension/plinking_duck/plinking_duck.duckdb_extension` | Loadable extension binary |

### Run Tests

```sh
make test
```

## Loading the Extension

When using the built DuckDB shell (`build/release/duckdb`), the extension is automatically loaded. For other DuckDB installations:

```sql
LOAD 'path/to/plinking_duck.duckdb_extension';
```

## First Query

Given a PLINK 2 fileset (`cohort.pgen`, `cohort.pvar`, `cohort.psam`):

```sql
-- Check what variants are in the dataset
SELECT CHROM, POS, ID, REF, ALT
FROM read_pvar('cohort.pvar')
LIMIT 10;
```

```sql
-- Check the sample list
SELECT IID, SEX
FROM read_psam('cohort.psam')
LIMIT 10;
```

```sql
-- Read genotypes for a specific region
SELECT CHROM, POS, ID, genotypes
FROM read_pgen('cohort.pgen')
WHERE CHROM = '1';
```

```sql
-- Tidy mode: one row per variant-sample pair
SELECT chrom, pos, id, iid, genotype
FROM read_pfile('cohort', tidy := true)
WHERE genotype IS NOT NULL
LIMIT 20;
```

## PLINK File Formats

PlinkingDuck supports both PLINK 2 and legacy PLINK 1 file formats. See the [File Formats](file-formats.md) reference for details on what's supported.
