# Function Reference

PlinkingDuck provides 9 SQL table functions: 4 file readers and 5 analysis functions.

## File Readers

| Function | Input | Description |
|----------|-------|-------------|
| [`read_pvar(path)`](read_pvar.md) | `.pvar` / `.bim` | Variant metadata (CHROM, POS, ID, REF, ALT) |
| [`read_psam(path)`](read_psam.md) | `.psam` / `.fam` | Sample metadata (FID, IID, SEX, phenotypes) |
| [`read_pgen(path)`](read_pgen.md) | `.pgen` | Binary genotypes as `LIST(TINYINT)` |
| [`read_pfile(prefix)`](read_pfile.md) | `.pgen` + `.pvar` + `.psam` | Complete fileset with tidy mode support |

## Analysis Functions

| Function | Input | Description |
|----------|-------|-------------|
| [`plink_freq(path)`](plink_freq.md) | `.pgen` | Per-variant allele frequencies and genotype counts |
| [`plink_hardy(path)`](plink_hardy.md) | `.pgen` | Hardy-Weinberg equilibrium exact test p-values |
| [`plink_missing(path)`](plink_missing.md) | `.pgen` | Per-variant or per-sample missingness rates |
| [`plink_ld(path)`](plink_ld.md) | `.pgen` | Pairwise linkage disequilibrium (r², D') |
| [`plink_score(path)`](plink_score.md) | `.pgen` | Polygenic risk scoring |

## Common Features

All functions support **projection pushdown**: columns not referenced in the query are not computed. For genotype-based functions, this means metadata-only queries skip genotype decoding entirely.

Analysis functions (except `plink_ld`) share a common set of named parameters documented in [Common Parameters](../common-parameters.md): `pvar`, `psam`, `samples`, and `region`.
