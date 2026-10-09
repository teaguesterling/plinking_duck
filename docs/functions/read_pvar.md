# read_pvar

Read PLINK 2 `.pvar` or legacy `.bim` variant metadata files.

## Synopsis

```sql
read_pvar(path VARCHAR | LIST(VARCHAR) [, with_ordinal := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `path` | `VARCHAR` or `LIST(VARCHAR)` | *(required)* | Path to a `.pvar` or `.bim` file, or a list of paths to read as one table (see [Multi-File Input](#multi-file-input)) |
| `with_ordinal` | `BOOLEAN` | `false` | When `true`, append a `__ordinal` `BIGINT` column (see [Named Parameters](#named-parameters)) |

The file format is auto-detected from the file contents. See [Common Parameters](../common-parameters.md#with_ordinal) for `with_ordinal`, shared with `read_pfile`.

### Named Parameters

#### `with_ordinal`

With `with_ordinal := true`, an extra `__ordinal` `BIGINT` column is appended to the schema, carrying the **0-based row ordinal in natural (file) order**. In a path list the ordinal **restarts at 0 for each file** (it is per-source, matching `read_pfile`'s per-shard variant ordinal), so a multi-file read produces per-file ordinals, not a single global sequence. This retires any need to assume `read_pvar` returns rows in file order: the ordinal is emitted explicitly.

### Multi-File Input

The argument may be a single path or a `LIST(VARCHAR)` of paths, mirroring `read_csv`/`read_json`. A list **row-concatenates** the variants from each file, in file order, into one logical table (each file's header is handled independently):

```sql
SELECT * FROM read_pvar(['data/chr1.pvar', 'data/chr2.pvar', 'data/chr3.pvar']);
```

This is the pure-variant analogue of [`read_pfile`'s multi-file input](read_pfile.md#multi-file-input), for the variant-sharded layout. `.pvar` and `.bim` files may be mixed within a list.

## Output Columns

### `.pvar` Format

Columns are determined by the file header. The standard columns are:

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | `VARCHAR` | Chromosome |
| `POS` | `INTEGER` | Base-pair position (1-based) |
| `ID` | `VARCHAR` | Variant identifier |
| `REF` | `VARCHAR` | Reference allele |
| `ALT` | `VARCHAR` | Alternate allele(s) |

Optional columns that may be present in the header:

| Column | Type | Description |
|--------|------|-------------|
| `QUAL` | `FLOAT` | Variant quality score |
| `FILTER` | `VARCHAR` | Filter status |
| `INFO` | `VARCHAR` | Additional information |
| `CM` | `DOUBLE` | Centimorgan position |

### `.bim` Format

Output is normalized to match `.pvar` column ordering:

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | `VARCHAR` | Chromosome |
| `POS` | `INTEGER` | Base-pair position (1-based) |
| `ID` | `VARCHAR` | Variant identifier |
| `REF` | `VARCHAR` | Reference allele |
| `ALT` | `VARCHAR` | Alternate allele |
| `CM` | `DOUBLE` | Centimorgan position |

!!! note
    The `.bim` file stores columns in a different order (CHROM, ID, CM, POS, ALT, REF). PlinkingDuck normalizes this so queries work identically on both formats.

## Description

`read_pvar` parses text-based variant metadata files. It auto-detects the format:

- If the first non-comment line starts with `#CHROM`, the file is parsed as `.pvar` (tab-delimited, dynamic columns from header).
- Otherwise, the file is parsed as `.bim` (whitespace-delimited, 6 fixed columns).

Lines starting with `##` are skipped as comment/meta lines. Dot values (`.`) in any column are returned as SQL `NULL`.

Projection pushdown is supported -- only columns referenced in the query are parsed.

## Examples

```sql
-- Read all variants
SELECT * FROM read_pvar('data/example.pvar');
```

```sql
-- Filter by chromosome
SELECT CHROM, POS, ID, REF, ALT
FROM read_pvar('data/example.pvar')
WHERE CHROM = '22';
```

```sql
-- Count variants per chromosome
SELECT CHROM, COUNT(*) AS n_variants
FROM read_pvar('data/example.pvar')
GROUP BY CHROM
ORDER BY CHROM;
```

```sql
-- Read a legacy .bim file (same interface)
SELECT * FROM read_pvar('data/example.bim');
```

```sql
-- Append the natural-order 0-based ordinal
SELECT __ordinal, CHROM, POS, ID
FROM read_pvar('data/example.pvar', with_ordinal := true);
```

## See Also

- [read_psam](read_psam.md) -- read sample metadata
- [read_pgen](read_pgen.md) -- read genotype data (uses `.pvar` as companion)
- [File Formats](../file-formats.md) -- format specifications
