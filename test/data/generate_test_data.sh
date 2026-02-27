#!/bin/bash
# Generate PLINK2 test data from example.vcf
# Requires: plink2 in PATH (or set PLINK2 environment variable)
set -euo pipefail
cd "$(dirname "$0")"

PLINK2="${PLINK2:-plink2}"

# Generate main test files: pgen_example.pgen, pgen_example.pvar, pgen_example.psam
"$PLINK2" --vcf example.vcf --make-pgen --out pgen_example --allow-extra-chr

# Generate all-missing genotype test files
cat > all_missing.vcf <<'EOF'
##fileformat=VCFv4.3
##contig=<ID=1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
1	1000	rs_miss1	A	G	.	.	.	GT	./.	./.
1	2000	rs_miss2	C	T	.	.	.	GT	./.	./.
EOF
"$PLINK2" --vcf all_missing.vcf --make-pgen --out all_missing --allow-extra-chr

# Generate orphan .pgen (no .psam companion) by copying main files without .psam
cp pgen_example.pgen pgen_orphan.pgen
cp pgen_example.pvar pgen_orphan.pvar

# Generate large dataset: 3000 variants × 8 samples
# Genotype pattern cycles: 0/0, 0/1, 1/1, ./. across samples and variants
# This exercises multi-batch scanning (>2048 variants) and parallel scan (MaxThreads > 1)
LARGE_VCF="large_example.vcf"
NUM_VARIANTS=3000
NUM_SAMPLES=8
SAMPLES=""
for s in $(seq 1 $NUM_SAMPLES); do
    SAMPLES="${SAMPLES}\tSAMP${s}"
done

{
    echo "##fileformat=VCFv4.3"
    for c in $(seq 1 3); do
        echo "##contig=<ID=${c},length=10000000>"
    done
    echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT${SAMPLES}\n"

    # Genotype cycle: 0/0=0, 0/1=1, 1/1=2, ./.=missing
    GT_VALS=("0/0" "0/1" "1/1" "./.")
    REFS=("A" "C" "G" "T")
    ALTS=("G" "T" "A" "C")

    # plink2 requires contiguous chromosomes — group 1000 variants per chrom
    v=1
    for chrom in 1 2 3; do
        for i in $(seq 1 1000); do
            pos=$(( i * 100 ))
            allele_idx=$(( (v - 1) % 4 ))
            id="var${v}"
            line="${chrom}\t${pos}\t${id}\t${REFS[$allele_idx]}\t${ALTS[$allele_idx]}\t.\t.\t.\tGT"
            for s in $(seq 0 $(( NUM_SAMPLES - 1 )) ); do
                gt_idx=$(( (v + s) % 4 ))
                line="${line}\t${GT_VALS[$gt_idx]}"
            done
            printf "${line}\n"
            v=$(( v + 1 ))
        done
    done
} > "$LARGE_VCF"

"$PLINK2" --vcf "$LARGE_VCF" --make-pgen --out large_example --allow-extra-chr
rm -f "$LARGE_VCF"

echo "Test data generated successfully."
echo "Generated files:"
ls -la pgen_example.pgen pgen_example.pvar pgen_example.psam \
       all_missing.pgen all_missing.pvar all_missing.psam \
       pgen_orphan.pgen pgen_orphan.pvar \
       large_example.pgen large_example.pvar large_example.psam
