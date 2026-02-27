# read_pvar

Read PLINK 2 `.pvar` or legacy `.bim` variant metadata files.

## Synopsis

```sql
read_pvar(path VARCHAR) -> TABLE
```

## Parameters

| Name | Type | Description |
|------|------|-------------|
| `path` | `VARCHAR` | Path to a `.pvar` or `.bim` file |

No named parameters. The file format is auto-detected from the file contents.

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

## See Also

- [read_psam](read_psam.md) -- read sample metadata
- [read_pgen](read_pgen.md) -- read genotype data (uses `.pvar` as companion)
- [File Formats](../file-formats.md) -- format specifications
