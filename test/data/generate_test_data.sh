#!/bin/bash
# Generate PLINK2 test data from example.vcf
# Requires: plink2 in PATH (or set PLINK2 environment variable)
set -euo pipefail
cd "$(dirname "$0")"

PLINK2="${PLINK2:-plink2}"

# Generate main test files: example.pgen, example.pvar, example.psam
"$PLINK2" --vcf example.vcf --make-pgen --out example --allow-extra-chr

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

echo "Test data generated successfully."
echo "Generated files:"
ls -la example.pgen example.pvar example.psam all_missing.pgen all_missing.pvar all_missing.psam
