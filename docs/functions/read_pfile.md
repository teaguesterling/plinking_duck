# read_pfile

Read a complete PLINK 2 fileset (`.pgen` + `.pvar` + `.psam`) from a common prefix. Supports multiple orient modes, sample subsetting, variant filtering, and region filtering.

## Synopsis

```sql
read_pfile(prefix VARCHAR [, pgen := ..., pvar := ..., psam := ...,
           orient := ..., genotypes := ..., phased := ..., dosages := ...,
           samples := ..., variants := ..., region := ...,
           af_range := ..., ac_range := ..., genotype_range := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `prefix` | `VARCHAR` | *(required)* | Common file prefix (e.g., `'data/cohort'`) |
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
| `genotype_range` | `STRUCT(min TINYINT, max TINYINT)` | All | Filter individual genotype values |

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

### Genotype Output Formats

The `genotypes` parameter controls how genotype data is structured:

- **`'array'`** (default): fixed-size `ARRAY(TINYINT, N)` — fastest, requires compile-time size
- **`'list'`**: variable-length `LIST(TINYINT)` — use for large cohorts (>100K samples)
- **`'columns'`**: one scalar `TINYINT` column per sample (variant orient) or per variant (sample orient); column names from IIDs or variant IDs. Incompatible with `orient := 'genotype'`.

### Phased Output

With `phased := true`, genotype elements change from `TINYINT` to `ARRAY(TINYINT, 2)` representing `[allele1, allele2]` haplotype pairs. `[0, 0]` = hom ref, `[0, 1]` = REF|ALT, `[1, 0]` = ALT|REF, `[1, 1]` = hom alt, NULL = missing. Incompatible with `dosages := true`. See [Genotype Encoding](../genotype-encoding.md).

### Filter Pushdown

`af_range` and `ac_range` filter variants by allele frequency or count using fast genotype counting (no decompression). `genotype_range` filters individual genotype values, setting non-matching values to NULL. See [Common Parameters](../common-parameters.md).

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
