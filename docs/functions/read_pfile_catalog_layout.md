# read_pfile_catalog_layout

Return the **interpreted chunk map** of a `.pfile.json` catalog: **one row per grid cell, including
declared holes**.

This is the *structural* view. The function computes the grid (the variants-major product of the two
axes' blocks) and resolves each cell's paths from the templates. It does **no pruning and no
verification, and it opens no data file**. It does not check that a path exists, and it reads no
`.pvar` / `.psam` / `.pgen`. A catalog whose cell files are absent still returns its full layout.
For cells you can read, use [`read_pfile_catalog`](read_pfile_catalog.md). For the literal declared
schema, use [`read_pfile_catalog_metadata`](read_pfile_catalog_metadata.md).

## Synopsis

```sql
read_pfile_catalog_layout(path VARCHAR) -> TABLE
```

## Output Columns

| Column | Type | Description |
|--------|------|-------------|
| `keys` | `MAP(VARCHAR, VARCHAR)` | Every key of both axes and its value for this cell. Empty for a version-1 catalog. |
| `variant_start`, `variant_stop` | `BIGINT` | The variant-block's half-open ordinal span `[start, stop)`. |
| `sample_start`, `sample_stop` | `BIGINT` | The sample-block's half-open ordinal span `[start, stop)`. |
| `pos_min`, `pos_max` | `BIGINT` | The variant-block's inclusive POS bound from the variants-axis `ordinal_interval` over `POS` key. `NULL` when the variants axis has no such key. |
| `present` | `BOOLEAN` | `false` for a declared hole, `true` otherwise. |
| `prefix` | `VARCHAR` | The cell's `.pgen` path, no extension. Absolute, normalized, **unchecked**. `NULL` for a hole. |
| `pvar` | `VARCHAR` | The variant-block's shared `.pvar` path. Template-resolved, **unchecked**. `NULL` for a hole. |
| `psam` | `VARCHAR` | The sample-block's shared `.psam` path. Template-resolved, **unchecked**. `NULL` for a hole. |

Rows are in variants-major order. A hole keeps its keys and spans but has no paths: the reader never
resolves a path for a hole.

## Examples

The whole grid, with the holes marked:

```sql
SELECT keys, variant_start, variant_stop, sample_start, sample_stop, present
FROM read_pfile_catalog_layout('cohort.pfile.json');
```

Only the declared holes:

```sql
SELECT keys FROM read_pfile_catalog_layout('cohort.pfile.json') WHERE NOT present;
```

Which cells cover a POS window on chromosome 1:

```sql
SELECT keys, prefix FROM read_pfile_catalog_layout('cohort.pfile.json')
WHERE keys['CHROM'] = '1' AND pos_min <= 20000 AND pos_max >= 10000 AND present;
```

## Version 1 and single catalogs

The function does not raise on an older catalog. Each shard is one **present** cell with an empty
`keys` map. The span on the concat axis is the shard's `partition`, and the other axis span is `NULL`
(it would need a read of the files). A shared companion is folded in, as for `read_pfile_catalog`.
`pos_min` / `pos_max` come from the shard's declared `range.POS`. A single-fileset catalog is one
cell with no spans, because its extent is not declared.

## Errors

A malformed catalog raises `catalog parse error`, as in `read_pfile_catalog`. That includes two blocks
that resolve to the same companion path, because the layout would be ambiguous. A missing catalog file
raises `file not found`. No error is ever raised for a missing **cell** file.

## Metadata or layout?

| | [`read_pfile_catalog_metadata`](read_pfile_catalog_metadata.md) | `read_pfile_catalog_layout` |
|---|---|---|
| View | Literal declared schema | Interpreted chunk map |
| One row per | Declared key | Grid cell, holes included |
| Computes the grid | No | Yes (structural only) |
| Opens data files | No | No |
