# read_pfile_catalog_metadata

Read a `.pfile.json` catalog **literally** and return its **declared schema**: one row per
declared partition key.

This is the *literal* view. It reports what the catalog file says and nothing else. It does not
tile the grid, prune, verify, or open any `.pgen` / `.pvar` / `.psam`. It reads the `.pfile.json`
only. For the computed grid, use [`read_pfile_catalog_layout`](read_pfile_catalog_layout.md). For the
shards or cells themselves, use [`read_pfile_catalog`](read_pfile_catalog.md).

## Synopsis

```sql
read_pfile_catalog_metadata(path VARCHAR) -> TABLE
```

## Output Columns

| Column | Type | Description |
|--------|------|-------------|
| `dataset` | `VARCHAR` | The catalog's `dataset` name. |
| `format` | `VARCHAR` | `'matrix'` (version 2), `'variants'` or `'samples'` (a version-1 multi-file catalog, named by its `concat_dim`), or `'single'` (a single-fileset catalog). |
| `axis` | `VARCHAR` | `'variants'` or `'samples'`: the axis the key belongs to. For a version-1 multi-file catalog, the concat axis. `NULL` for `'single'`. |
| `key` | `VARCHAR` | The key name as declared, for example `CHROM`. `NULL` for version 1. |
| `type` | `VARCHAR` | `'enum'`, `'ordinal_interval'` or `'value_bin'`. `NULL` for version 1. |
| `over` | `VARCHAR` | What an interval key ranges over: `'POS'` or `'sample_index'` for `ordinal_interval`, the binned value for `value_bin`, `'CHROM'` for an `enum` key marked as the contig. `NULL` when the key declares no `over`. |
| `axis_position` | `INTEGER` | 0-based position of the key within its axis, in declared order (outer to inner). `NULL` for version 1. |

## Examples

The declared keys of a version-2 catalog:

```sql
SELECT axis, key, type, over, axis_position
FROM read_pfile_catalog_metadata('cohort.pfile.json');
```

| axis | key | type | over | axis_position |
|------|-----|------|------|---------------|
| variants | CHROM | enum | NULL | 0 |
| variants | CHUNK | ordinal_interval | POS | 1 |
| samples | GROUP | enum | NULL | 0 |
| samples | SLICE | ordinal_interval | sample_index | 1 |
| samples | age | value_bin | age | 2 |

Find the keys that carry an interval or a bin:

```sql
SELECT axis, key FROM read_pfile_catalog_metadata('cohort.pfile.json') WHERE type <> 'enum';
```

## Version 1 and single catalogs

A version-1 catalog has no partition keys. The function does not raise. It returns **one row**:
`dataset` and `format` are set (`'variants'`, `'samples'` or `'single'`), `axis` is the concat axis
(`NULL` for a single fileset), and `key`, `type`, `over` and `axis_position` are `NULL`.

## Metadata or layout?

| | `read_pfile_catalog_metadata` | [`read_pfile_catalog_layout`](read_pfile_catalog_layout.md) |
|---|---|---|
| View | Literal declared schema | Interpreted chunk map |
| One row per | Declared key | Grid cell, holes included |
| Computes the grid | No | Yes (structural only) |
| Opens data files | No | No |
