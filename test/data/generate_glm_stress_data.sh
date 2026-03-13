#!/bin/bash
# Generate large synthetic PLINK2 dataset for plink_glm stress testing
# Creates: glm_stress.pgen/pvar/psam (100K variants × 1000 samples)
# Also creates: glm_stress_phenotypes.parquet, glm_stress_covariates.parquet
# Requires: plink2, python3, and a built DuckDB with plinking_duck extension
set -euo pipefail
cd "$(dirname "$0")"

PLINK2="${PLINK2:-/mnt/aux-data/teague/Dev/spack/opt/spack/linux-zen3/plink2-2.0.0-a.6.9-e2l3wx22vy6tkmdwzhaexlml4rs3okx3/bin/plink2}"
DUCKDB="${DUCKDB:-../../build/release/duckdb}"

NUM_VARIANTS=100000
NUM_SAMPLES=1000
PREFIX="glm_stress"

echo "Generating ${NUM_VARIANTS} variants × ${NUM_SAMPLES} samples..."

# --- Generate VCF using Python for speed ---
python3 -c "
import random
random.seed(42)

num_variants = ${NUM_VARIANTS}
num_samples = ${NUM_SAMPLES}
chroms = [(str(c), num_variants // 5) for c in range(1, 6)]  # 5 chroms, 20K each

samples = [f'SAMP{i:04d}' for i in range(1, num_samples + 1)]
refs = ['A', 'C', 'G', 'T']
alts = ['G', 'T', 'A', 'C']
gts = ['0/0', '0/1', '1/1', './.']
gt_weights = [0.25, 0.40, 0.30, 0.05]  # hom-ref, het, hom-alt, missing

print('##fileformat=VCFv4.3')
for c, _ in chroms:
    print(f'##contig=<ID={c},length=100000000>')
print('##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">')
print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples))

v = 1
for chrom, count in chroms:
    for i in range(1, count + 1):
        pos = i * 100
        ai = (v - 1) % 4
        genotypes = random.choices(gts, weights=gt_weights, k=num_samples)
        line = f'{chrom}\t{pos}\tvar{v}\t{refs[ai]}\t{alts[ai]}\t.\t.\t.\tGT\t' + '\t'.join(genotypes)
        print(line)
        v += 1
" > "${PREFIX}.vcf"

echo "Converting VCF to PGEN..."
"$PLINK2" --vcf "${PREFIX}.vcf" --make-pgen --out "$PREFIX" --allow-extra-chr
rm -f "${PREFIX}.vcf"

echo "Generating parquet files..."

# --- Generate phenotype and covariate parquet files using DuckDB ---
"$DUCKDB" -unsigned <<SQL
-- Phenotype file: IID, pheno_quant (continuous), pheno_binary (case/control)
COPY (
    SELECT
        'SAMP' || LPAD(CAST(i AS VARCHAR), 4, '0') AS IID,
        -- Continuous phenotype: normal-ish distribution
        ROUND(2.5 + 1.5 * (random() + random() + random() - 1.5), 4) AS pheno_quant,
        -- Binary phenotype: ~40% cases
        CASE WHEN random() < 0.4 THEN 1 ELSE 0 END AS pheno_binary
    FROM generate_series(1, ${NUM_SAMPLES}) AS t(i)
) TO '${PREFIX}_phenotypes.parquet' (FORMAT PARQUET);

-- Covariate file: IID, age, bmi, pc1-pc5
COPY (
    SELECT
        'SAMP' || LPAD(CAST(i AS VARCHAR), 4, '0') AS IID,
        ROUND(20.0 + 50.0 * random(), 1) AS age,
        ROUND(18.0 + 15.0 * random(), 1) AS bmi,
        ROUND(random() * 2 - 1, 6) AS pc1,
        ROUND(random() * 2 - 1, 6) AS pc2,
        ROUND(random() * 2 - 1, 6) AS pc3,
        ROUND(random() * 2 - 1, 6) AS pc4,
        ROUND(random() * 2 - 1, 6) AS pc5
    FROM generate_series(1, ${NUM_SAMPLES}) AS t(i)
) TO '${PREFIX}_covariates.parquet' (FORMAT PARQUET);

SELECT 'Phenotypes:' AS info, COUNT(*) AS rows FROM '${PREFIX}_phenotypes.parquet';
SELECT 'Covariates:' AS info, COUNT(*) AS rows FROM '${PREFIX}_covariates.parquet';
SQL

echo ""
echo "Generated files:"
ls -lh "${PREFIX}".{pgen,pvar,psam} "${PREFIX}_phenotypes.parquet" "${PREFIX}_covariates.parquet"
echo ""
echo "Done! ${NUM_VARIANTS} variants × ${NUM_SAMPLES} samples"
