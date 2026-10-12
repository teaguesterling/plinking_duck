# plink_pca

Principal component analysis of the sample genetic relationship matrix (GRM) —
the standard tool for inferring population structure and deriving ancestry
covariates for GWAS.

!!! warning "Conditional on Eigen3"
    `plink_pca` is the **only** function that is not always present. It needs
    [Eigen3](../development.md#optional-dependency-eigen3), and CMake compiles it
    out (with a configure-time warning) when Eigen3 is not found. Binaries
    released from this repo are built through vcpkg, which supplies Eigen3, so
    they always have it; a local `make` on a machine without `libeigen3-dev` does
    not. Calling it on a build that lacks it raises
    `Catalog Error: Table Function with name plink_pca does not exist!`. To check
    a build before relying on it:

    ```sql
    SELECT count(*) > 0 AS has_pca
    FROM duckdb_functions()
    WHERE function_name = 'plink_pca';
    ```

## Synopsis

```sql
plink_pca(path VARCHAR [, mode := ..., n_pcs := ...,
          pvar := ..., psam := ..., samples := ..., region := ...,
          selection := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `path` | `VARCHAR` | *(required)* | Path to the `.pgen` file |
| `mode` | `VARCHAR` | `'samples'` | Output shape: `'samples'`, `'pcs'`, or `'both'` (see [Output](#output-modes)) |
| `n_pcs` | `INTEGER` | `10` | Number of principal components to compute (must be `>= 1` and strictly less than the sample count) |
| `pvar` | `VARCHAR` | Auto-discovered | Path to `.pvar` or `.bim` file |
| `psam` | `VARCHAR` | Auto-discovered | Path to `.psam` or `.fam` file |
| `samples` | `LIST(VARCHAR)` or `LIST(INTEGER)` | All | Subset to specific samples |
| `region` | `VARCHAR` | All | Restrict to a genomic region (`chr:start-end`) |
| `selection` | `STRUCT` | None | Composable per-axis [selection](../common-parameters.md#selection) struct |

See [Common Parameters](../common-parameters.md) for details on `pvar`, `psam`,
`samples`, `region`, and the composable `selection` struct.

## Output modes

The `mode` parameter selects one of three output shapes.

### `mode := 'samples'` (default)

One row per sample — the per-sample projections onto each PC. This is the shape
you join back to phenotypes and pass as GWAS covariates.

| Column | Type | Description |
|--------|------|-------------|
| `FID` | `VARCHAR` | Family ID |
| `IID` | `VARCHAR` | Individual ID |
| `PC1` … `PCk` | `DOUBLE` | Projection onto each of the `n_pcs` components |

### `mode := 'pcs'`

One row per principal component — the spectrum and how much variance each PC
explains.

| Column | Type | Description |
|--------|------|-------------|
| `PC` | `INTEGER` | Component index (1-based) |
| `EIGENVALUE` | `DOUBLE` | Eigenvalue of the component |
| `VARIANCE_PROPORTION` | `DOUBLE` | Fraction of total variance this PC explains |
| `CUMULATIVE_VARIANCE` | `DOUBLE` | Running sum of `VARIANCE_PROPORTION` |

### `mode := 'both'`

A single row carrying both results as nested values, for callers that want the
projections and the spectrum from one scan.

| Column | Type | Description |
|--------|------|-------------|
| `EIGENVEC` | `LIST(STRUCT(FID VARCHAR, IID VARCHAR, PC1 DOUBLE, …, PCk DOUBLE))` | Per-sample projections |
| `EIGENVAL` | `LIST(DOUBLE)` | The `n_pcs` eigenvalues |

## Description

`plink_pca` computes the top `n_pcs` eigenvectors/eigenvalues of the sample GRM
using **randomized subspace iteration** (an approximate solver), mirroring
`plink2 --pca approx`. The numerics are not bit-identical to plink2 — if you need
exact agreement, validate against `plink2 --pca approx` on a fixture before
relying on it. Monomorphic variants are dropped; each remaining variant is
mean-centered and variance-normalized before the GRM is formed.

### Constraints

- `n_pcs` must be `>= 1` and **strictly less than** the (post-subset) sample
  count; too few variants or samples for the requested `n_pcs` raises a clear
  error (the approximate solver needs headroom above `n_pcs`).
- The intermediate projection matrix is bounded by
  [`plinking_max_matrix_elements`](read_pfile.md) — a very large
  `n_pcs` × sample count can exceed it, with an error that tells you to reduce
  `n_pcs`/variant count or raise the limit:

  ```sql
  SET plinking_max_matrix_elements = 20_000_000_000;
  ```

## Examples

```sql
-- Default: per-sample projections onto 10 PCs
SELECT * FROM plink_pca('data/cohort.pgen');
```

```sql
-- First 5 PCs, as GWAS covariates joined to phenotypes
SELECT p.IID, p.PC1, p.PC2, p.PC3, p.PC4, p.PC5
FROM plink_pca('data/cohort.pgen', n_pcs := 5) AS p;
```

```sql
-- The variance each component explains (scree data)
SELECT PC, VARIANCE_PROPORTION, CUMULATIVE_VARIANCE
FROM plink_pca('data/cohort.pgen', mode := 'pcs', n_pcs := 20)
ORDER BY PC;
```

```sql
-- PCA within an LD-pruned region and a sample subset
SELECT * FROM plink_pca('data/cohort.pgen',
    n_pcs := 10,
    selection := struct_pack(region := '1:1-50000000',
                             samples := ['S1', 'S2', 'S3']));
```

## See Also

- [Common Parameters](../common-parameters.md) — `samples`, `region`, `selection`
- [Function Reference](index.md) — why `plink_pca` is conditional on Eigen3
- [Development](../development.md#optional-dependency-eigen3) — building with Eigen3
