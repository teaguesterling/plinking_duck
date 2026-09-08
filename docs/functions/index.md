# Function Reference

PlinkingDuck provides 12 SQL table functions: 5 file readers and 7 analysis functions. Eleven of them are always present. `plink_pca` is the exception — it is registered only when the build found Eigen3.

## File Readers

| Function | Input | Description |
|----------|-------|-------------|
| [`read_pvar(path)`](read_pvar.md) | `.pvar` / `.bim` (single or list) | Variant metadata (CHROM, POS, ID, REF, ALT) |
| [`read_psam(path)`](read_psam.md) | `.psam` / `.fam` | Sample metadata (FID, IID, SEX, phenotypes) |
| [`read_pgen(path)`](read_pgen.md) | `.pgen` | Binary genotypes as `ARRAY(TINYINT, N)` |
| [`read_pfile(prefix)`](read_pfile.md) | `.pgen` + `.pvar` + `.psam` (single or list) | Complete fileset with orient mode support |
| [`read_plink_vcf(path)`](read_plink_vcf.md) | `.vcf` / `.vcf.gz` | Fast biallelic genotype extraction from VCF |

## Analysis Functions

| Function | Input | Description |
|----------|-------|-------------|
| [`plink_freq(path)`](plink_freq.md) | `.pgen` | Per-variant allele frequencies and genotype counts |
| [`plink_hardy(path)`](plink_hardy.md) | `.pgen` | Hardy-Weinberg equilibrium exact test p-values |
| [`plink_missing(path)`](plink_missing.md) | `.pgen` | Per-variant or per-sample missingness rates |
| [`plink_ld(path)`](plink_ld.md) | `.pgen` | Pairwise linkage disequilibrium (r², D') |
| [`plink_score(path)`](plink_score.md) | `.pgen` | Polygenic risk scoring |
| [`plink_glm(prefix)`](plink_glm.md) | pfile prefix | Per-variant GWAS regression (linear, logistic, Firth) |
| `plink_pca(path)` ⚠️ | `.pgen` | Principal component analysis. **Conditional — see below.** |

### `plink_pca` is conditional on Eigen3

`plink_pca` is the only function that is not always there. It needs
[Eigen3](../development.md#optional-dependency-eigen3), and CMake compiles it
out (with a configure-time warning) when Eigen3 is not found. Binaries released
from this repo are built through vcpkg, which supplies Eigen3, so they always
have it; a local `make` on a machine without `libeigen3-dev` does not.

Calling it on a build that lacks it raises:

```
Catalog Error: Table Function with name plink_pca does not exist!
```

To check a build before relying on it:

```sql
SELECT count(*) = 1 AS has_pca
FROM duckdb_functions()
WHERE function_name = 'plink_pca';
```

## Common Features

All functions support **projection pushdown**: columns not referenced in the query are not computed. For genotype-based functions, this means metadata-only queries skip genotype decoding entirely.

Analysis functions (except `plink_ld`) share a common set of named parameters documented in [Common Parameters](../common-parameters.md): `pvar`, `psam`, `samples`, and `region`.
