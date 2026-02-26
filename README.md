# PlinkingDuck

A DuckDB extension for reading [PLINK 2](https://www.cog-genomics.org/plink/2.0/) genomics file formats directly in SQL.

PlinkingDuck brings PLINK genotype, variant, and sample data into DuckDB, letting you query genomics datasets with standard SQL instead of format-specific tools.

## Usage

### read_pvar(path)

Read PLINK2 `.pvar` (variant information) or legacy `.bim` files into SQL-queryable tables.
Format is auto-detected: files with a `#CHROM` header line are parsed as `.pvar`,
files without a header are parsed as legacy `.bim` (6-column, tab-delimited).

```sql
-- Read a .pvar file
SELECT * FROM read_pvar('path/to/file.pvar');

-- Read a legacy .bim file (auto-detected, columns normalized to .pvar order)
SELECT * FROM read_pvar('path/to/file.bim');

-- Filter variants by chromosome
SELECT CHROM, POS, ID, REF, ALT
FROM read_pvar('data/variants.pvar')
WHERE CHROM = '22';
```

**Output columns (.pvar):** CHROM, POS, ID, REF, ALT, and any optional columns
present in the header (QUAL, FILTER, INFO, CM).

**Output columns (.bim):** CHROM, POS, ID, REF, ALT, CM (normalized from
the .bim file order of CHROM, ID, CM, POS, ALT, REF).

Dot values (`.`) in any column are returned as `NULL`. Projection pushdown is
supported — only columns referenced in the query are parsed.

## Building
### Managing dependencies
DuckDB extensions uses VCPKG for dependency management. Enabling VCPKG is very simple: follow the [installation instructions](https://vcpkg.io/en/getting-started) or just run the following:
```shell
git clone https://github.com/Microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh
export VCPKG_TOOLCHAIN_PATH=`pwd`/vcpkg/scripts/buildsystems/vcpkg.cmake
```
Note: VCPKG is only required for extensions that want to rely on it for dependency management. If you want to develop an extension without dependencies, or want to do your own dependency management, just skip this step. Note that the example extension uses VCPKG to build with a dependency for instructive purposes, so when skipping this step the build may not work without removing the dependency.

### Build steps
Now to build the extension, run:
```sh
make
```
The main binaries that will be built are:
```sh
./build/release/duckdb
./build/release/test/unittest
./build/release/extension/plinking_duck/plinking_duck.duckdb_extension
```
- `duckdb` is the binary for the duckdb shell with the extension code automatically loaded.
- `unittest` is the test runner of duckdb. Again, the extension is already linked into the binary.
- `plinking_duck.duckdb_extension` is the loadable binary as it would be distributed.

### `read_pfile(prefix [, pgen, pvar, psam, tidy, samples, variants, region])`

Read a complete PLINK 2 fileset (`.pgen` + `.pvar` + `.psam`) from a common prefix.
Supports tidy output mode, sample subsetting, variant filtering, and region filtering.

```sql
-- Read all variants (one row per variant, genotypes as list)
SELECT * FROM read_pfile('data/example');

-- Tidy mode (one row per variant x sample)
SELECT chrom, pos, iid, genotype
FROM read_pfile('data/example', tidy := true);

-- Sample subset (by name or 0-based index)
SELECT * FROM read_pfile('data/example',
    tidy := true, samples := ['SAMPLE1', 'SAMPLE3']);

-- Region filter
SELECT * FROM read_pfile('data/example', region := '22:1-50000000');

-- Variant filter
SELECT * FROM read_pfile('data/example', variants := ['rs1', 'rs2']);

-- Combined filters
SELECT iid, genotype
FROM read_pfile('data/example', tidy := true,
    region := '1:10000-50000', samples := ['SAMPLE1']);
```

**Default mode output:** CHROM, POS, ID, REF, ALT, genotypes `LIST(TINYINT)` — same as `read_pgen`.

**Tidy mode output:** CHROM, POS, ID, REF, ALT, plus all `.psam` columns (FID, IID, SEX, etc.),
plus scalar `genotype` (TINYINT). One row per variant x sample combination.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `pgen` | VARCHAR | Explicit `.pgen` path (overrides prefix) |
| `pvar` | VARCHAR | Explicit `.pvar`/`.bim` path (overrides prefix) |
| `psam` | VARCHAR | Explicit `.psam`/`.fam` path (overrides prefix) |
| `tidy` | BOOLEAN | Tidy output: one row per (variant, sample). Default: false |
| `samples` | LIST(VARCHAR) or LIST(INTEGER) | Filter to specific samples by IID or 0-based index |
| `variants` | LIST(VARCHAR) or LIST(INTEGER) | Filter to specific variants by ID or 0-based index |
| `region` | VARCHAR | Filter to genomic region (`chr:start-end` or `chr`) |

**Projection pushdown:** Genotype decoding is skipped when genotype columns
are not referenced in the query, making metadata-only queries fast.

### `read_pgen(path [, pvar, psam, samples])`

Read PLINK 2 `.pgen` binary genotype files. Auto-discovers companion `.pvar` and `.psam` files.

```sql
-- Basic usage (auto-discovers .pvar and .psam)
SELECT * FROM read_pgen('path/to/file.pgen');

-- Explicit companion file paths
SELECT * FROM read_pgen('file.pgen', pvar := 'other.pvar', psam := 'other.psam');

-- Only variant metadata (skips genotype decoding — fast)
SELECT chrom, pos, id FROM read_pgen('file.pgen');

-- Subset samples by 0-based index
SELECT id, genotypes FROM read_pgen('file.pgen', samples := [0, 2]);
```

**Output columns:** CHROM, POS, ID, REF, ALT (from `.pvar`), genotypes (from `.pgen`).

**Genotype encoding:** `LIST(TINYINT)` with one entry per sample:
0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate, NULL = missing.

**Projection pushdown:** If the `genotypes` column is not referenced in the query,
genotype decoding is skipped entirely, making metadata-only queries very fast.

**Optional .psam:** The `.psam` file is optional for `read_pgen`. Without it,
genotype lists are sized from the `.pgen` header and the `samples` parameter
only accepts integer indices.

### `read_psam(path)`

Read PLINK 2 `.psam` (sample information) or PLINK 1 `.fam` files as DuckDB tables.

```sql
-- Read a .psam file
SELECT * FROM read_psam('samples.psam');
```

**Supported formats:**

| Format | Header | Columns |
|--------|--------|---------|
| `.psam` (`#FID` header) | `#FID  IID  SEX  [phenotypes...]` | FID + IID + dynamic columns |
| `.psam` (`#IID` header) | `#IID  SEX  [phenotypes...]` | IID + dynamic columns (no FID) |
| `.fam` (no header) | Fixed 6 columns | FID, IID, PAT, MAT, SEX, PHENO1 |

**Type mapping:**

| Column | DuckDB Type | Notes |
|--------|-------------|-------|
| SEX | INTEGER | 1=male, 2=female; `0`/`NA`/`.` → NULL |
| PAT, MAT | VARCHAR | `0`/`.`/`NA` → NULL (unknown parent) |
| All others | VARCHAR | `.`/`NA` → NULL |

**Examples:**

```sql
-- Count samples per family
SELECT FID, COUNT(*) AS n FROM read_psam('cohort.psam') GROUP BY FID;

-- Filter by sex
SELECT IID FROM read_psam('cohort.psam') WHERE SEX = 2;

-- Join sample info with phenotype data
SELECT s.IID, s.SEX, p.BMI
FROM read_psam('cohort.psam') s
JOIN 'phenotypes.csv' p ON s.IID = p.IID;

-- Read a legacy .fam file
SELECT * FROM read_psam('cohort.fam');
```

## Running the extension
To run the extension code, simply start the shell with `./build/release/duckdb`.

Now we can use the features from the extension directly in DuckDB:
```sql
D SELECT CHROM, POS, ID, REF, ALT FROM read_pvar('test/data/example.pvar') LIMIT 3;
┌─────────┬───────┬─────────┬─────────┬─────────┐
│  CHROM  │  POS  │   ID    │   REF   │   ALT   │
│ varchar │ int32 │ varchar │ varchar │ varchar │
├─────────┼───────┼─────────┼─────────┼─────────┤
│ 1       │ 10000 │ rs1     │ A       │ G       │
│ 1       │ 20000 │ rs2     │ C       │ T       │
│ 1       │ 30000 │ rs3     │ G       │ A,C     │
└─────────┴───────┴─────────┴─────────┴─────────┘
```


## Running the tests
Different tests can be created for DuckDB extensions. The primary way of testing DuckDB extensions should be the SQL tests in `./test/sql`. These SQL tests can be run using:
```sh
make test
```

### Installing the deployed binaries
To install your extension binaries from S3, you will need to do two things. Firstly, DuckDB should be launched with the
`allow_unsigned_extensions` option set to true. How to set this will depend on the client you're using. Some examples:

CLI:
```shell
duckdb -unsigned
```

Python:
```python
con = duckdb.connect(':memory:', config={'allow_unsigned_extensions' : 'true'})
```

NodeJS:
```js
db = new duckdb.Database(':memory:', {"allow_unsigned_extensions": "true"});
```

Secondly, you will need to set the repository endpoint in DuckDB to the HTTP url of your bucket + version of the extension
you want to install. To do this run the following SQL query in DuckDB:
```sql
SET custom_extension_repository='bucket.s3.eu-west-1.amazonaws.com/<your_extension_name>/latest';
```
Note that the `/latest` path will allow you to install the latest extension version available for your current version of
DuckDB. To specify a specific version, you can pass the version instead.

After running these steps, you can install and load your extension using the regular INSTALL/LOAD commands in DuckDB:
```sql
INSTALL plinking_duck;
LOAD plinking_duck;
```

## Setting up CLion

### Opening project
Configuring CLion with this extension requires a little work. Firstly, make sure that the DuckDB submodule is available.
Then make sure to open `./duckdb/CMakeLists.txt` (so not the top level `CMakeLists.txt` file from this repo) as a project in CLion.
Now to fix your project path go to `tools->CMake->Change Project Root`([docs](https://www.jetbrains.com/help/clion/change-project-root-directory.html)) to set the project root to the root dir of this repo.

### Debugging
To set up debugging in CLion, there are two simple steps required. Firstly, in `CLion -> Settings / Preferences -> Build, Execution, Deploy -> CMake` you will need to add the desired builds (e.g. Debug, Release, RelDebug, etc). There's different ways to configure this, but the easiest is to leave all empty, except the `build path`, which needs to be set to `../build/{build type}`, and CMake Options to which the following flag should be added, with the path to the extension CMakeList:

```
-DDUCKDB_EXTENSION_CONFIGS=<path_to_the_exentension_CMakeLists.txt>
```

The second step is to configure the unittest runner as a run/debug configuration. To do this, go to `Run -> Edit Configurations` and click `+ -> Cmake Application`. The target and executable should be `unittest`. This will run all the DuckDB tests. To specify only running the extension specific tests, add `--test-dir ../../.. [sql]` to the `Program Arguments`. Note that it is recommended to use the `unittest` executable for testing/development within CLion. The actual DuckDB CLI currently does not reliably work as a run target in CLion.
