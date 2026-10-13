# Function Reference

PlinkingDuck provides 13 core SQL table functions: 6 file readers and 7 analysis functions. (`plink_pca` additionally requires Eigen3 at build time — a hard dependency supplied by `vcpkg.json`, so standard builds always include it; see [its page](plink_pca.md).) In addition there is a utility table function ([`pfile_partitions`](pfile_partitions.md)) and a scalar feature-detection function ([`plinking_duck_capabilities`](plinking_duck_capabilities.md)).

## File Readers

| Function | Input | Description |
|----------|-------|-------------|
| [`read_pvar(path)`](read_pvar.md) | `.pvar` / `.bim` (single or list) | Variant metadata (CHROM, POS, ID, REF, ALT) |
| [`read_psam(path)`](read_psam.md) | `.psam` / `.fam` | Sample metadata (FID, IID, SEX, phenotypes) |
| [`read_pgen(path)`](read_pgen.md) | `.pgen` | Binary genotypes as `ARRAY(TINYINT, N)` |
| [`read_pfile(prefix)`](read_pfile.md) | `.pgen` + `.pvar` + `.psam` (single or list) | Complete fileset with orient mode support |
| [`read_plink_vcf(path)`](read_plink_vcf.md) | `.vcf` / `.vcf.gz` | Fast biallelic genotype extraction from VCF |
| [`read_howard_manifest(path)`](read_howard_manifest.md) | `.howard.json` catalog | One row per surviving shard of a sharded pfile set (verify + prune) |

## Analysis Functions

| Function | Input | Description |
|----------|-------|-------------|
| [`plink_freq(path)`](plink_freq.md) | `.pgen` | Per-variant allele frequencies and genotype counts |
| [`plink_hardy(path)`](plink_hardy.md) | `.pgen` | Hardy-Weinberg equilibrium exact test p-values |
| [`plink_missing(path)`](plink_missing.md) | `.pgen` | Per-variant or per-sample missingness rates |
| [`plink_ld(path)`](plink_ld.md) | `.pgen` | Pairwise linkage disequilibrium (r², D') |
| [`plink_score(path)`](plink_score.md) | `.pgen` | Polygenic risk scoring |
| [`plink_glm(prefix)`](plink_glm.md) | pfile prefix | Per-variant GWAS regression (linear, logistic, Firth) |
| [`plink_pca(path)`](plink_pca.md) | `.pgen` | Principal component analysis (requires Eigen3 at build time) |

## Utility Functions

| Function | Input | Description |
|----------|-------|-------------|
| [`pfile_partitions(prefix)`](pfile_partitions.md) | pfile prefix (single or list) | Half-open `[lo, hi)` ordinal-range partitions over the variants dimension, for chunked/parallel reads |

## Scalar Functions

| Function | Returns | Description |
|----------|---------|-------------|
| [`plinking_duck_capabilities()`](plinking_duck_capabilities.md) | `VARCHAR[]` | Capability tokens this build supports, for SQL feature detection |

!!! note "`plink_pca` requires Eigen3 at build time"
    Eigen3 is a hard build dependency (declared in `vcpkg.json`), so every
    standard build — including the released binaries — includes `plink_pca`. A
    build that cannot find Eigen3 fails at configure time rather than silently
    omitting it. See [its page](plink_pca.md).

## Writing filesets

Beyond reading, PlinkingDuck can **write** PLINK 2 filesets back out of SQL
through DuckDB **COPY** formats (not table functions, so they are not in the table
above):

```sql
COPY (SELECT * FROM read_pfile('in', genotypes := 'struct')) TO 'out' (FORMAT pfile);
```

`FORMAT pfile` writes a whole fileset (`.pvar` + `.pgen` + `.psam`); `FORMAT pvar`
/ `pgen` / `psam` write the individual files. See
[Writing filesets](write_pfile.md) for the options, the sample-identity contract,
local-only / atomic-overwrite behavior, and the current limitations (biallelic
hardcalls only).

The writer can also **shard** a cohort into one fileset per key with `shard_by`,
emitting a [`.howard.json` catalog](../formats/howard.md) that
[`read_howard_manifest`](read_howard_manifest.md) reads back as a single logical
dataset — see [Writing filesets](write_pfile.md#sharded-output-shard_by).

## Common Features

All functions support **projection pushdown**: columns not referenced in the query are not computed. For genotype-based functions, this means metadata-only queries skip genotype decoding entirely.

Analysis functions share a common set of named parameters documented in [Common Parameters](../common-parameters.md): `pvar`, `psam`, `samples`, `region`, and the composable `selection` struct. (`plink_ld` takes only the sample axis from `selection`, and reads its region from the `region` parameter directly.)
