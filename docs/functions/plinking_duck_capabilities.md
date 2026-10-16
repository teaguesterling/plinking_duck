# plinking_duck_capabilities

Return the list of **capability tokens** this build of PlinkingDuck supports, for SQL feature detection. This is a **scalar** function (not a table function).

## Synopsis

```sql
plinking_duck_capabilities() -> VARCHAR[]
```

## Description

The result is a `LIST(VARCHAR)` of capability tokens, computed from the **live catalog at bind time** — not a hardcoded list. A token is emitted only when its underlying registration fact is genuinely present in the build, so the function can never claim a feature the readers do not actually honor. Tokens are spelled exactly as the underlying named-parameter key or table-function name.

On a build carrying the phase-2 partitions work the tokens are:

| Token | Meaning |
|-------|---------|
| `variant_range` | [`read_pfile`](read_pfile.md#variant-range-ordinal-space) accepts `variant_range := [lo, hi)` |
| `with_ordinal` | [`read_pfile`](read_pfile.md) and [`read_pvar`](read_pvar.md) accept `with_ordinal := true` |
| `pfile_partitions` | the [`pfile_partitions`](pfile_partitions.md) table function exists |

A build that predates these features returns a list that omits the corresponding tokens (an older build with no capability surface returns an empty list).

## Feature detection

Test for a single capability with `list_contains`:

```sql
SELECT list_contains(plinking_duck_capabilities(), 'variant_range');
```

A consumer that may be running against an **older build that lacks this function entirely** should first check the catalog (a plain read, safe inside an open transaction) before evaluating the scalar — calling an absent function raises a `BinderException`:

```sql
-- Present only when this detection surface itself exists:
SELECT count(*) > 0 AS has_capability_surface
FROM duckdb_functions()
WHERE function_name = 'plinking_duck_capabilities';
```

If `has_capability_surface` is false, treat the build as capability-less and fall back. This two-step shape (catalog read, then — only if present — evaluate) is exactly how the `ducks_in_a_row`/`pfile_catalog` consumers detect PlinkingDuck features at rewrite time.

## Optional: behavioral probing via scalarfs (third-party pattern)

> This is an **optional** pattern for a third party; the PlinkingDuck readers do not depend on it, and neither do `ducks_in_a_row`/`pfile_catalog`. `scalarfs` is **not** a dependency.

When a build is old enough to lack `plinking_duck_capabilities()`, a consumer can still probe a capability *behaviorally* — hermetically, in SQL, with no temp file — if the [`scalarfs`](https://github.com/teaguesterling/duckdb_scalarfs) extension is also loaded. PlinkingDuck's readers go through DuckDB's **FileSystem**, so a scalarfs inline-content protocol (which serves a literal payload where a path is expected) can stand in for a file path. Schematically — pass an inline `.pvar` literal to a reader and ask for the feature:

```
-- schematic (exact scalarfs protocol spelling omitted):
read_pvar('<scalarfs-inline-pvar-literal>', with_ordinal := true)
```

A build lacking the parameter errors at bind; a build that honors it returns the emitted ordinal, with no file created. The exact inline-protocol syntax is scalarfs's to document; the point here is only that the readers' use of the FileSystem makes such a hermetic probe possible. `scalarfs`'s `chore/bump-duckdb-1.5.6` branch pins the **same** DuckDB submodule commit (`069cc9f9b5`) as this PlinkingDuck branch, so the two are on the same ABI line and can be co-loaded once both are built. (Prefer `plinking_duck_capabilities()` whenever it is present — it is exact and needs no second extension.)

## See Also

- [pfile_partitions](pfile_partitions.md), [read_pfile](read_pfile.md), [read_pvar](read_pvar.md) — the features these tokens gate
- [Function Reference](index.md)
