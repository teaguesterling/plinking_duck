# read_psam

Read PLINK 2 `.psam` or legacy `.fam` sample metadata files.

## Synopsis

```sql
read_psam(path VARCHAR) -> TABLE
```

## Parameters

| Name | Type | Description |
|------|------|-------------|
| `path` | `VARCHAR` | Path to a `.psam` or `.fam` file |

No named parameters. The file format is auto-detected from the file contents.

## Output Columns

### `.psam` Format (with `#FID` header)

| Column | Type | Description |
|--------|------|-------------|
| `FID` | `VARCHAR` | Family ID |
| `IID` | `VARCHAR` | Individual ID |
| `SEX` | `INTEGER` | Sex (1=male, 2=female; 0/NA/. = NULL) |
| *additional* | `VARCHAR` | Any extra columns from the header (phenotypes, etc.) |

### `.psam` Format (with `#IID` header)

| Column | Type | Description |
|--------|------|-------------|
| `IID` | `VARCHAR` | Individual ID |
| `SEX` | `INTEGER` | Sex (1=male, 2=female; 0/NA/. = NULL) |
| *additional* | `VARCHAR` | Any extra columns from the header |

### `.fam` Format

| Column | Type | Description |
|--------|------|-------------|
| `FID` | `VARCHAR` | Family ID |
| `IID` | `VARCHAR` | Individual ID |
| `PAT` | `VARCHAR` | Paternal ID (0/. = NULL) |
| `MAT` | `VARCHAR` | Maternal ID (0/. = NULL) |
| `SEX` | `INTEGER` | Sex (1=male, 2=female; 0 = NULL) |
| `PHENO1` | `VARCHAR` | Primary phenotype (as string) |

## Description

`read_psam` parses text-based sample metadata files. It auto-detects the format:

- If the first line starts with `#FID`, the file is parsed as `.psam` with family ID.
- If the first line starts with `#IID`, the file is parsed as `.psam` without family ID.
- Otherwise, the file is parsed as `.fam` (whitespace-delimited, 6 fixed columns).

### Missing Value Handling

| Column | Missing Values | Mapped To |
|--------|---------------|-----------|
| `SEX` | `0`, `NA`, `.`, empty | `NULL` |
| `PAT`, `MAT` | `0`, `NA`, `.`, empty | `NULL` |
| All others | `NA`, `.`, empty | `NULL` |

!!! note
    In `.fam` files, PLINK conventionally uses `-9` as a missing phenotype sentinel. PlinkingDuck leaves this as the string `"-9"` rather than mapping it to NULL, since this is a PLINK-specific convention. Filter with `WHERE PHENO1 != '-9'` if needed.

Projection pushdown is supported.

## Examples

```sql
-- Read all sample information
SELECT * FROM read_psam('data/example.psam');
```

```sql
-- Count samples per family
SELECT FID, COUNT(*) AS n
FROM read_psam('data/cohort.psam')
GROUP BY FID;
```

```sql
-- Filter by sex (female samples)
SELECT IID FROM read_psam('data/cohort.psam')
WHERE SEX = 2;
```

```sql
-- Join with external phenotype data
SELECT s.IID, s.SEX, p.BMI
FROM read_psam('data/cohort.psam') s
JOIN 'phenotypes.csv' p ON s.IID = p.IID;
```

```sql
-- Read a legacy .fam file
SELECT * FROM read_psam('data/example.fam');
```

## See Also

- [read_pvar](read_pvar.md) -- read variant metadata
- [read_pgen](read_pgen.md) -- read genotype data (uses `.psam` as companion)
- [File Formats](../file-formats.md) -- format specifications
