# PlinkingDuck

A DuckDB extension for reading [PLINK 2](https://www.cog-genomics.org/plink/2.0/) genomics file formats and running common genetic analyses directly in SQL.

PlinkingDuck brings PLINK genotype, variant, and sample data into DuckDB, letting you query genomics datasets with standard SQL instead of format-specific command-line tools. Read files, filter variants, compute allele frequencies, test Hardy-Weinberg equilibrium, measure missingness, calculate linkage disequilibrium, and run polygenic scoring — all from a SQL prompt.

!!! tip "Built on the shoulders of [DuckHTS](https://github.com/RGenomicsETL/duckhts)"

    DuckHTS pioneered the idea of querying genomics file formats directly in DuckDB using htslib. It supports VCF, BCF, BAM, CRAM, FASTA, FASTQ, GTF, GFF, and tabix-indexed files — the sequencing side of the house. PlinkingDuck picks up where DuckHTS leaves off, covering the PLINK genotype file formats used in population genetics and GWAS. If you work with both sequencing and genotype data, you'll want both extensions.

## Quick Example

```sql
-- Load the extension
LOAD plinking_duck;

-- Read variant metadata
SELECT CHROM, POS, ID, REF, ALT
FROM read_pvar('cohort.pvar')
WHERE CHROM = '22';

-- Read genotypes with sample info
SELECT chrom, pos, iid, genotype
FROM read_pfile('cohort', tidy := true)
WHERE genotype = 2;

-- Compute allele frequencies
SELECT ID, ALT_FREQ
FROM plink_freq('cohort.pgen')
WHERE ALT_FREQ > 0.01;

-- Hardy-Weinberg equilibrium test
SELECT ID, P_HWE
FROM plink_hardy('cohort.pgen')
WHERE P_HWE < 1e-6;
```

## Functions

PlinkingDuck provides **4 file readers** and **5 analysis functions**:

### File Readers

| Function | Description |
|----------|-------------|
| [`read_pvar(path)`](functions/read_pvar.md) | Read `.pvar` or `.bim` variant metadata |
| [`read_psam(path)`](functions/read_psam.md) | Read `.psam` or `.fam` sample metadata |
| [`read_pgen(path)`](functions/read_pgen.md) | Read `.pgen` binary genotype files |
| [`read_pfile(prefix)`](functions/read_pfile.md) | Read a complete PLINK fileset with tidy mode |

### Analysis Functions

| Function | Description |
|----------|-------------|
| [`plink_freq(path)`](functions/plink_freq.md) | Per-variant allele frequencies |
| [`plink_hardy(path)`](functions/plink_hardy.md) | Hardy-Weinberg equilibrium exact test |
| [`plink_missing(path)`](functions/plink_missing.md) | Per-variant or per-sample missingness |
| [`plink_ld(path)`](functions/plink_ld.md) | Pairwise linkage disequilibrium |
| [`plink_score(path)`](functions/plink_score.md) | Polygenic risk scoring |

## Features

- **Native SQL access** to PLINK 2 binary genotype data
- **Projection pushdown** -- skip genotype decoding when only metadata is queried
- **Parallel scan** -- multi-threaded variant processing for large files
- **Sample subsetting** -- filter to specific samples by IID or index
- **Region filtering** -- restrict to genomic regions (`chr:start-end`)
- **Legacy format support** -- read PLINK 1 `.bim` and `.fam` files
- **VFS integration** -- works with DuckDB's filesystem abstraction
