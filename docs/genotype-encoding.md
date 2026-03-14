# Genotype Encoding

All genotype data in PlinkingDuck uses `TINYINT` values with a simple additive encoding.

## Encoding Table

| Value | Meaning | VCF Equivalent |
|-------|---------|----------------|
| `0` | Homozygous reference | `0/0` |
| `1` | Heterozygous | `0/1` |
| `2` | Homozygous alternate | `1/1` |
| `NULL` | Missing genotype | `./.` |

The value represents the **count of alternate alleles**: 0 copies, 1 copy, or 2 copies.

## Contexts

Genotypes appear in two contexts depending on the function and mode:

### List Context

In `read_pgen` and `read_pfile` default mode, genotypes are returned as an `ARRAY(TINYINT, N)` column named `genotypes`, where N is the sample count. Each array element corresponds to one sample, in the same order as the `.psam` file.

```sql
-- Each row has a list of genotypes, one per sample
SELECT ID, genotypes FROM read_pgen('data.pgen');
-- rs1  [0, 1, 2, NULL]
-- rs2  [0, 0, 1, 1]
```

### Scalar Context

In `read_pfile` genotype orient mode (`orient := 'genotype'`), each genotype is a scalar `TINYINT` column named `genotype`. There is one row per variant-sample combination.

```sql
-- Each row is one variant-sample pair
SELECT ID, IID, genotype
FROM read_pfile('data', orient := 'genotype');
-- rs1  SAMPLE1  0
-- rs1  SAMPLE2  1
-- rs1  SAMPLE3  2
-- rs1  SAMPLE4  NULL
```

### Columns Context

With `genotypes := 'columns'`, each sample (in variant orient) or each variant (in sample orient) gets its own scalar `TINYINT` column. Column names come from sample IIDs or variant IDs.

```sql
-- One column per sample
SELECT ID, SAMPLE1, SAMPLE2, SAMPLE3
FROM read_pfile('data', genotypes := 'columns');
-- rs1  0  1  2
-- rs2  1  1  0
```

### Phased Haplotype Output

With `phased := true`, genotype elements change from scalar `TINYINT` to `ARRAY(TINYINT, 2)` representing `[allele1, allele2]` haplotype pairs:

| Value | Meaning | VCF Equivalent |
|-------|---------|----------------|
| `[0, 0]` | Homozygous reference | `0\|0` |
| `[0, 1]` | Ref/Alt (or unphased het) | `0\|1` |
| `[1, 0]` | Alt/Ref | `1\|0` |
| `[1, 1]` | Homozygous alternate | `1\|1` |
| `NULL` | Missing | `.\|.` |

```sql
SELECT ID, genotypes
FROM read_pfile('data', phased := true);
-- rs1  [[0, 0], [0, 1], [1, 1], NULL]
```

Phased output works across all orient modes and genotypes modes (array, list). Incompatible with `dosages := true`.

## Working with Genotypes

### Counting alternate alleles per sample

```sql
SELECT IID, SUM(genotype) AS total_alt_alleles
FROM read_pfile('data', orient := 'genotype')
WHERE genotype IS NOT NULL
GROUP BY IID;
```

### Filtering for heterozygous calls

```sql
SELECT ID, IID
FROM read_pfile('data', orient := 'genotype')
WHERE genotype = 1;
```

### Computing allele frequency from genotypes

```sql
SELECT ID,
    AVG(genotype) / 2.0 AS alt_freq
FROM read_pfile('data', orient := 'genotype')
WHERE genotype IS NOT NULL
GROUP BY ID;
```

!!! note
    For efficient allele frequency computation, prefer [`plink_freq()`](functions/plink_freq.md) which uses pgenlib's fast counting path without full genotype decompression.

## Internal Representation

PLINK 2 `.pgen` files store genotypes as 2-bit packed values (4 genotypes per byte). The raw encoding is:

| 2-bit value | Meaning |
|-------------|---------|
| `00` | Homozygous reference |
| `01` | Heterozygous |
| `10` | Homozygous alternate |
| `11` | Missing |

PlinkingDuck converts these to the user-facing `0/1/2/NULL` encoding during reads. Missing values (`11` = `-9` in pgenlib's byte representation) are mapped to SQL `NULL`.
