# read_plink_vcf

Read genotypes from a VCF file using plink-ng's optimized GT parsing. Designed as a fast path for extracting biallelic genotypes — skips multiallelic variants and outputs in the same format as `read_pfile()`.

For full VCF/BCF parsing (multiallelic, INFO fields, annotations), use [DuckHTS](https://github.com/teaguesterling/duckhts) `read_bcf()` instead.

## Synopsis

```sql
read_plink_vcf(path VARCHAR [, genotypes := ..., phased := ...,
               region := ..., min_gq := ..., min_dp := ..., max_dp := ...,
               halfcall := ...]) -> TABLE
```

## Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `path` | `VARCHAR` | *(required)* | Path to `.vcf` or `.vcf.gz` file |
| `genotypes` | `VARCHAR` | `'auto'` | Output format: `'array'`, `'list'`, or `'columns'` |
| `phased` | `BOOLEAN` | `false` | Output phased haplotype pairs (`ARRAY(TINYINT, 2)`) |
| `region` | `VARCHAR` | All | Filter to genomic region (`chr` or `chr:start-end`) |
| `min_gq` | `INTEGER` | -1 (disabled) | Minimum genotype quality (GQ); samples below → NULL |
| `min_dp` | `INTEGER` | -1 (disabled) | Minimum read depth (DP); samples below → NULL |
| `max_dp` | `INTEGER` | -1 (disabled) | Maximum read depth (DP); samples above → NULL |
| `halfcall` | `VARCHAR` | `'missing'` | Half-call handling: `'missing'`, `'reference'`, `'haploid'`, or `'error'` |

### Genotype Output Modes

- **`'auto'`** — `'array'` if sample count ≤ 100,000, else `'list'`
- **`'array'`** — `ARRAY(TINYINT, N)` where N = sample count (fixed-size, fastest)
- **`'list'`** — `LIST(TINYINT)` (variable-size, no sample count limit)
- **`'columns'`** — one `TINYINT` column per sample, named by sample ID

### Quality Filtering

Quality filters apply per-sample: if a sample's GQ or DP falls outside the specified range, that sample's genotype is set to NULL (missing) for that variant. The variant itself is still output.

### Phased Mode

When `phased := true`, genotype elements become `ARRAY(TINYINT, 2)` representing `[allele1, allele2]`:
- `[0, 0]` = hom ref
- `[0, 1]` = ref|alt (standard het)
- `[1, 0]` = alt|ref (phase-flipped het)
- `[1, 1]` = hom alt
- `NULL` = missing

## Output Columns

| Column | Type | Description |
|--------|------|-------------|
| `CHROM` | `VARCHAR` | Chromosome |
| `POS` | `INTEGER` | Position (1-based) |
| `ID` | `VARCHAR` | Variant identifier (NULL if `.`) |
| `REF` | `VARCHAR` | Reference allele |
| `ALT` | `VARCHAR` | Alternate allele |
| `genotypes` | `ARRAY(TINYINT, N)` | Per-sample genotype values |

In `'columns'` mode, the `genotypes` column is replaced by one column per sample.

### Genotype Encoding

| Value | Meaning |
|-------|---------|
| 0 | Homozygous reference |
| 1 | Heterozygous |
| 2 | Homozygous alternate |
| NULL | Missing genotype |

## Behavior

- **Multiallelic variants are skipped** — a warning reports the count at scan completion
- **GT must be the first FORMAT subfield** — throws an error otherwise
- **Compressed VCF** (`.vcf.gz`, bgzf) is transparently decompressed via DuckDB VFS
- **Projection pushdown** — genotype parsing is skipped entirely when only metadata columns are referenced

## Examples

```sql
-- Basic variant count
SELECT COUNT(*) FROM read_plink_vcf('cohort.vcf.gz');

-- Genotype distribution
SELECT
  genotypes[1] as gt,
  COUNT(*) as n
FROM read_plink_vcf('sample.vcf')
GROUP BY gt;

-- High-quality variants on chr1
SELECT *
FROM read_plink_vcf('wgs.vcf.gz', region := 'chr1', min_gq := 30, min_dp := 10)
LIMIT 100;

-- Phased haplotypes
SELECT id, genotypes
FROM read_plink_vcf('phased.vcf', phased := true);

-- Per-sample columns for easy pivoting
SELECT SAMPLE1, SAMPLE2
FROM read_plink_vcf('trio.vcf', genotypes := 'columns');

-- Cross-validate VCF against pfile
SELECT COUNT(*) FROM (
  SELECT v.id
  FROM read_plink_vcf('data.vcf') v
  JOIN read_pfile('data') p ON v.id = p.id
  WHERE v.genotypes != p.genotypes
);
```

## Limitations

- Only biallelic variants (multiallelic are skipped with a warning)
- No INFO/QUAL/FILTER field access (use DuckHTS for that)
- Single-threaded sequential scan (no parallel chunking)
- `genotypes := 'struct'`, `'counts'`, and `'stats'` modes are not supported
- BCF (binary VCF) is not supported — only text VCF and gzipped VCF
