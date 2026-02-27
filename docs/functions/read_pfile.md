# read_pfile

Read a complete PLINK 2 fileset (`.pgen` + `.pvar` + `.psam`) from a common prefix. Supports tidy output mode, sample subsetting, variant filtering, and region filtering.

## Synopsis

```sql
read_pfile(prefix VARCHAR [, pgen := ..., pvar := ..., psam := ...,
           tidy := ..., samples := ..., variants := ..., region := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `prefix` | `VARCHAR` | *(required)* | Common file prefix (e.g., `'data/cohort'`) |
| `pgen` | `VARCHAR` | `prefix.pgen` | Explicit `.pgen` path |
| `pvar` | `VARCHAR` | `prefix.pvar` | Explicit `.pvar` or `.bim` path |
| `psam` | `VARCHAR` | `prefix.psam` | Explicit `.psam` or `.fam` path |
| `tidy` | `BOOLEAN` | `false` | Output one row per variant-sample pair |
| `samples` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All | Filter to specific samples |
| `variants` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All | Filter to specific variants |
| `region` | `VARCHAR` | All | Filter to genomic region (`chr:start-end`) |

See [Common Parameters](../common-parameters.md) for details on `samples` and `region`.

### Reserved Parameters (Not Yet Implemented)

| Name | Type | Description |
|------|------|-------------|
| `dosages` | `BOOLEAN` | Include dosage data (raises error if `true`) |
| `phased` | `BOOLEAN` | Include phasing data (raises error if `true`) |

### Variant Filtering

The `variants` parameter accepts:

- **`LIST(VARCHAR)`** -- variant IDs (matched against ID in `.pvar`)
- **`LIST(INTEGER)`** -- 0-based variant indices

When both `region` and `variants` are specified, the result is the intersection (variants must match both filters).

## Output Columns

### Default Mode (`tidy := false`)

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | `VARCHAR` | Chromosome |
| `POS` | `INTEGER` | Base-pair position |
| `ID` | `VARCHAR` | Variant identifier |
| `REF` | `VARCHAR` | Reference allele |
| `ALT` | `VARCHAR` | Alternate allele |
| `genotypes` | `LIST(TINYINT)` | Genotype calls, one per sample |

### Tidy Mode (`tidy := true`)

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

In tidy mode, each row represents one variant-sample combination. The sample metadata columns come from the `.psam` file and vary depending on what columns are present.

## Description

`read_pfile` is the most feature-rich reader function. It combines genotype, variant, and sample data with flexible filtering.

### File Discovery

By default, `read_pfile` constructs file paths by appending extensions to the prefix: `prefix.pgen`, `prefix.pvar`, `prefix.psam`. Each can be overridden individually.

### Projection Pushdown

Genotype decoding is skipped when genotype columns (`genotypes` or `genotype`) are not referenced in the query.

### Parallel Scan

Default mode uses multi-threaded parallel scanning. Tidy mode is single-threaded because it expands each variant into multiple rows (one per sample).

## Examples

```sql
-- Read all variants in default mode
SELECT * FROM read_pfile('data/example');
```

```sql
-- Tidy mode: one row per variant-sample pair
SELECT chrom, pos, iid, genotype
FROM read_pfile('data/example', tidy := true);
```

```sql
-- Sample subset by name
SELECT * FROM read_pfile('data/example',
    tidy := true,
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
    tidy := true,
    variants := ['rs1', 'rs2']);
```

```sql
-- Combined filters
SELECT iid, genotype
FROM read_pfile('data/example',
    tidy := true,
    region := '1:10000-50000',
    samples := ['SAMPLE1']);
```

## See Also

- [read_pgen](read_pgen.md) -- lower-level genotype reader
- [read_pvar](read_pvar.md) -- variant metadata only
- [read_psam](read_psam.md) -- sample metadata only
- [Genotype Encoding](../genotype-encoding.md) -- encoding reference
