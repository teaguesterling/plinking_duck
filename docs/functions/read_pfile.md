# read_pfile

Read a complete PLINK 2 fileset (`.pgen` + `.pvar` + `.psam`) from a common prefix. Supports multiple orient modes, sample subsetting, variant filtering, and region filtering.

## Synopsis

```sql
read_pfile(prefix VARCHAR | LIST(VARCHAR) [, pgen := ..., pvar := ..., psam := ...,
           orient := ..., genotypes := ..., phased := ..., dosages := ...,
           samples := ..., variants := ..., region := ...,
           af_range := ..., ac_range := ...,
           include_genotypes := ..., genotype_range := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `prefix` | `VARCHAR` or `LIST(VARCHAR)` | *(required)* | Common file prefix (e.g., `'data/cohort'`), or a list of prefixes to read as one table (see [Multi-File Input](#multi-file-input)) |
| `pgen` | `VARCHAR` | `prefix.pgen` | Explicit `.pgen` path |
| `pvar` | `VARCHAR` | `prefix.pvar` | Explicit `.pvar` or `.bim` path |
| `psam` | `VARCHAR` | `prefix.psam` | Explicit `.psam` or `.fam` path |
| `orient` | `VARCHAR` | `'variant'` | Output orientation: `'variant'`, `'genotype'`, or `'sample'` |
| `genotypes` | `VARCHAR` | `'array'` | Genotype output format: `'array'`, `'list'`, or `'columns'` |
| `phased` | `BOOLEAN` | `false` | Output phased haplotype pairs instead of dosage counts |
| `dosages` | `BOOLEAN` | `false` | Output dosage values instead of hard calls |
| `samples` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All | Filter to specific samples |
| `variants` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All | Filter to specific variants |
| `region` | `VARCHAR` | All | Filter to genomic region (`chr:start-end`) |
| `af_range` | `STRUCT(min DOUBLE, max DOUBLE)` | All | Filter variants by allele frequency range |
| `ac_range` | `STRUCT(min INTEGER, max INTEGER)` | All | Filter variants by allele count range |
| `include_genotypes` | `LIST(VARCHAR)` | All | Filter by hardcall category: `'hom_ref'`, `'het'`, `'hom_alt'`, `'missing'` (any subset; canonical genotype filter) |
| `genotype_range` | `STRUCT(min TINYINT, max TINYINT [, include_missing BOOLEAN])` | All | Numeric alias of `include_genotypes` for contiguous ranges over [0, 2] |

See [Common Parameters](../common-parameters.md) for details on `samples`, `region`, and filter parameters.

### Variant Filtering

The `variants` parameter accepts:

- **`LIST(VARCHAR)`** -- variant IDs (matched against ID in `.pvar`)
- **`LIST(INTEGER)`** -- 0-based variant indices

When both `region` and `variants` are specified, the result is the intersection (variants must match both filters).

## Output Columns

### Variant Orient Mode (Default)

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | `VARCHAR` | Chromosome |
| `POS` | `INTEGER` | Base-pair position |
| `ID` | `VARCHAR` | Variant identifier |
| `REF` | `VARCHAR` | Reference allele |
| `ALT` | `VARCHAR` | Alternate allele |
| `genotypes` | `ARRAY(TINYINT, N)` | Genotype calls, one per sample |

### Genotype Orient Mode (`orient := 'genotype'`)

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | `VARCHAR` | Chromosome |
| `POS` | `INTEGER` | Base-pair position |
| `ID` | `VARCHAR` | Variant identifier |
| `REF` | `VARCHAR` | Reference allele |
| `ALT` | `VARCHAR` | Alternate allele |
| `FID` | `VARCHAR` | Family ID (from `.psam`) |
| `IID` | `VARCHAR` | Individual ID (from `.psam`) |
| `SEX` | `INTEGER` | Sex (from `.psam`) |
| *additional* | varies | Any extra `.psam` columns |
| `genotype` | `TINYINT` | Single genotype value |

In genotype orient mode, each row represents one variant-sample combination. The sample metadata columns come from the `.psam` file and vary depending on what columns are present.

### Sample Orient Mode (`orient := 'sample'`)

| Column | Type | Description |
|--------|------|-------------|
| `FID` | `VARCHAR` | Family ID (from `.psam`) |
| `IID` | `VARCHAR` | Individual ID (from `.psam`) |
| `SEX` | `INTEGER` | Sex (from `.psam`) |
| *additional* | varies | Any extra `.psam` columns |
| `genotypes` | `ARRAY(TINYINT, M)` | Genotype calls, one per variant |

In sample orient mode, each row represents one sample. The genotypes array contains one value per variant (M = effective variant count after filtering). All genotypes are pre-read at bind time.

## Description

`read_pfile` is the most feature-rich reader function. It combines genotype, variant, and sample data with flexible filtering.

### File Discovery

By default, `read_pfile` constructs file paths by appending extensions to the prefix: `prefix.pgen`, `prefix.pvar`, `prefix.psam`. Each can be overridden individually.

### Multi-File Input

The first argument may be a single prefix or a `LIST(VARCHAR)` of prefixes, mirroring the way `read_csv` and `read_json` accept a list. A list **row-concatenates** the variants from each file, in file order, into one logical table:

```sql
SELECT * FROM read_pfile(['/data/chr1', '/data/chr2', '/data/chr3']);
```

This targets the biobank-scale layout where a fileset is **sharded by variant** — many `.pgen`/`.pvar` files that all share one identical `.psam`. All files in a list must carry the **same sample set: the same IIDs in the same order**. Alignment is **positional** and trust-the-caller — there is no per-file sample-identity check by design (adding one would reintroduce the `.psam`-parse cost on the real workload). The only guard is a free one: each `.pgen` header reports its sample count, and a file whose sample count differs from the first is rejected as a hard error. Sample order is *not* verified; that remains the caller's contract.

Sample and variant metadata are taken from the first file only (identical by contract). Row order across files is not guaranteed (same as today's single-file parallel scan) — add `ORDER BY` if you need a stable order.

Current scope and limitations:

- `variant` and `genotype` orient support multi-file input.
- `orient := 'sample'` with **multiple** files is **not yet supported** (a single file still works).
- `region :=` works across a list (each shard is filtered to the range; the union is returned).
- `variants := [...]` is **not yet supported** with a multi-file list (variant IDs/indices don't resolve globally across shards yet) — use `region :=` for selection across files. It still works on a single file.
- Explicit `pgen` / `pvar` / `psam` overrides are single-file only; combining them with a multi-element list is an error (they cannot disambiguate per shard).

The list of prefixes can come from anywhere, including a catalog macro that returns a `VARCHAR[]` of shard prefixes (e.g. via the scalarfs `pathmacro:` protocol) fed directly to `read_pfile([...])`.

### Genotype Output Formats

The `genotypes` parameter controls how genotype data is structured:

- **`'array'`** (default): fixed-size `ARRAY(TINYINT, N)` — fastest, requires compile-time size
- **`'list'`**: variable-length `LIST(TINYINT)` — use for large cohorts (>100K samples)
- **`'columns'`**: one scalar `TINYINT` column per sample (variant orient) or per variant (sample orient); column names from IIDs or variant IDs. Incompatible with `orient := 'genotype'`.

### Phased Output

With `phased := true`, genotype elements change from `TINYINT` to `ARRAY(TINYINT, 2)` representing `[allele1, allele2]` haplotype pairs. `[0, 0]` = hom ref, `[0, 1]` = REF|ALT, `[1, 0]` = ALT|REF, `[1, 1]` = hom alt, NULL = missing. Incompatible with `dosages := true`. See [Genotype Encoding](../genotype-encoding.md).

### Filter Pushdown

`af_range` and `ac_range` filter variants by allele frequency or count using fast genotype counting (no decompression). `include_genotypes` (and its numeric alias `genotype_range`) filters by hardcall category — in `variant` orient it sets non-matching values to NULL; in `genotype` and `sample` orient it drops non-matching rows, so a carrier query in `sample` orient materializes only the matching subjects. See [Common Parameters](../common-parameters.md).

### Projection Pushdown

Genotype decoding is skipped when genotype columns (`genotypes` or `genotype`) are not referenced in the query.

### Parallel Scan

Default mode uses multi-threaded parallel scanning. Genotype orient mode is single-threaded because it expands each variant into multiple rows (one per sample). Sample orient mode is single-threaded (all genotypes are pre-read at bind time).

## Examples

```sql
-- Read all variants in default mode
SELECT * FROM read_pfile('data/example');
```

```sql
-- Multi-file: read several variant-sharded filesets (identical .psam) as one table
SELECT COUNT(*) FROM read_pfile(['data/chr1', 'data/chr2', 'data/chr3']);
```

```sql
-- Genotype orient mode: one row per variant-sample pair
SELECT chrom, pos, iid, genotype
FROM read_pfile('data/example', orient := 'genotype');
```

```sql
-- Sample orient: one row per sample
SELECT IID, genotypes
FROM read_pfile('data/example', orient := 'sample');
```

```sql
-- Sample subset by name
SELECT * FROM read_pfile('data/example',
    orient := 'genotype',
    samples := ['SAMPLE1', 'SAMPLE3']);
```

```sql
-- Region filter
SELECT * FROM read_pfile('data/example',
    region := '22:1-50000000');
```

```sql
-- Variant filter by ID
SELECT * FROM read_pfile('data/example',
    orient := 'genotype',
    variants := ['rs1', 'rs2']);
```

```sql
-- Combined filters
SELECT iid, genotype
FROM read_pfile('data/example',
    orient := 'genotype',
    region := '1:10000-50000',
    samples := ['SAMPLE1']);
```

```sql
-- Phased haplotype output
SELECT ID, genotypes
FROM read_pfile('data/example', phased := true);
```

```sql
-- Filter by allele frequency (common variants only)
SELECT * FROM read_pfile('data/example',
    af_range := {min: 0.01, max: 0.5});
```

```sql
-- Columns pivot mode
SELECT * FROM read_pfile('data/example',
    genotypes := 'columns');
```

## See Also

- [read_pgen](read_pgen.md) -- lower-level genotype reader
- [read_pvar](read_pvar.md) -- variant metadata only
- [read_psam](read_psam.md) -- sample metadata only
- [Genotype Encoding](../genotype-encoding.md) -- encoding and phased output reference
