# File Formats

PlinkingDuck reads PLINK 2 and legacy PLINK 1 genomics file formats. This page describes the file formats and what PlinkingDuck supports.

## PLINK 2 Fileset

PLINK 2 uses a three-file fileset to represent genotyped samples. All three files share a common prefix (e.g., `cohort.pgen`, `cohort.pvar`, `cohort.psam`).

| File | Contents | Description |
|------|----------|-------------|
| `.pgen` | Binary genotypes | Compressed genotype matrix (variants x samples) |
| `.pvar` | Variant metadata | Tab-delimited: `#CHROM`, POS, ID, REF, ALT, plus optional columns |
| `.psam` | Sample metadata | Tab-delimited: `#FID`/`#IID`, IID, SEX, plus optional phenotypes |

### `.pvar` Format

The `.pvar` file is a tab-delimited text file with a header line starting with `#CHROM`.

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	10000	rs1	A	G	.	.	.
1	20000	rs2	C	T	30	PASS	DP=100
```

Required columns: `CHROM`, `POS`, `ID`, `REF`, `ALT`. Optional columns include `QUAL`, `FILTER`, `INFO`, and `CM`.

### `.psam` Format

The `.psam` file is a tab-delimited text file with a header line starting with either `#FID` or `#IID`.

```
#FID	IID	SEX
FAM001	SAMPLE1	1
FAM001	SAMPLE2	2
FAM002	SAMPLE3	0
```

The `#FID` header indicates both family ID and individual ID columns are present. The `#IID` header indicates no family ID column. Additional columns after `SEX` are treated as phenotype columns.

### `.pgen` Format

The `.pgen` file is a binary format storing compressed genotype data. Each entry represents a biallelic genotype encoded as a 2-bit value. PlinkingDuck reads `.pgen` files using the pgenlib library from the [plink-ng](https://github.com/chrchang/plink-ng) project.

## Legacy PLINK 1 Formats

PlinkingDuck supports reading PLINK 1 text-based metadata formats.

| PLINK 1 | PLINK 2 equivalent | Supported |
|---------|-------------------|-----------|
| `.bim` | `.pvar` | Yes -- via `read_pvar()` |
| `.fam` | `.psam` | Yes -- via `read_psam()` |
| `.bed` | `.pgen` | No -- use `plink2 --make-pgen` to convert |

### `.bim` Format

The `.bim` file is a whitespace-delimited text file with no header and 6 fixed columns:

```
1	rs1	0.5	10000	G	A
1	rs2	0.0	20000	T	C
```

File column order: CHROM, ID, CM, POS, ALT, REF.

PlinkingDuck normalizes `.bim` output to match `.pvar` column ordering: CHROM, POS, ID, REF, ALT, CM. This means queries work identically on both formats.

### `.fam` Format

The `.fam` file is a whitespace-delimited text file with no header and 6 fixed columns:

```
FAM001 SAMPLE1 0 0 1 -9
FAM001 SAMPLE2 0 0 2 -9
```

Column order: FID, IID, PAT (paternal ID), MAT (maternal ID), SEX, PHENO1.

Special missing value conventions:

- `0` in PAT/MAT columns means unknown parent (mapped to NULL)
- `0` in SEX means unknown sex (mapped to NULL)
- `-9` in PHENO1 is left as the string `"-9"` (not mapped to NULL, since this is a PLINK-specific convention)

## Auto-Detection

PlinkingDuck auto-detects file formats:

- **`.pvar` vs `.bim`**: If the first non-comment line starts with `#CHROM`, the file is parsed as `.pvar`. Otherwise, it's parsed as `.bim`.
- **`.psam` vs `.fam`**: If the first line starts with `#FID` or `#IID`, the file is parsed as `.psam`. Otherwise, it's parsed as `.fam`.

## Companion File Discovery

Functions that read `.pgen` files (`read_pgen`, `read_pfile`, `plink_freq`, etc.) automatically discover companion `.pvar` and `.psam` files by replacing the `.pgen` extension. The search order is:

1. `.pvar` then `.bim` (for variant metadata)
2. `.psam` then `.fam` (for sample metadata)

You can override auto-discovery with the `pvar` and `psam` named parameters.
