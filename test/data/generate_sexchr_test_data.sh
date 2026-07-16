#!/bin/bash
# Generate a small sex-chromosome test dataset exercising ploidy-aware MAF/HWE.
# Requires: plink2 in PATH (or set PLINK2 environment variable).
#
# Cohort: 3 males (M1,M2,M3), 3 females (F1,F2,F3).
# Variants span autosomal, chrX PAR (diploid both sexes), chrX non-PAR
# (males haploid), chrY (males only), and chrMT (haploid everyone).
#
# See test/sql/plink_sexchr.test for the hand-computed expected values.
set -euo pipefail
cd "$(dirname "$0")"

PLINK2="${PLINK2:-plink2}"

cat > sexchr_example.vcf <<'EOF'
##fileformat=VCFv4.3
##contig=<ID=1>
##contig=<ID=X>
##contig=<ID=Y>
##contig=<ID=MT>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	M1	M2	M3	F1	F2	F3
1	10000	a1	A	G	.	.	.	GT	0/1	0/0	1/1	0/0	0/1	1/1
X	1000000	par1	A	G	.	.	.	GT	0/1	0/0	1/1	0/0	0/1	0/0
X	50000000	x1	A	G	.	.	.	GT	1	0	0	0/1	0/0	1/1
Y	5000000	y1	C	T	.	.	.	GT	1	1	0	.	.	.
MT	500	mt1	A	G	.	.	.	GT	1	1	1	1	0	0
EOF

# Sex is required at import so plink2 can apply male-haploid ploidy on chrX/Y.
cat > sexchr_sex.txt <<'EOF'
#IID	SEX
M1	1
M2	1
M3	1
F1	2
F2	2
F3	2
EOF

# Import to .pgen/.pvar/.psam.
#  --split-par b38: relabel chrX pseudo-autosomal region variants to PAR1/PAR2
#                   so they stay diploid in males (GRCh38 boundaries).
#  --update-sex:    supply sex so male chrX non-PAR / chrY calls are haploid.
"$PLINK2" --vcf sexchr_example.vcf --update-sex sexchr_sex.txt --split-par b38 \
    --make-pgen --sort-vars --out sexchr_example --allow-extra-chr

# Ensure the .psam carries the SEX column (the ploidy-aware stats read it).
cat > sexchr_example.psam <<'EOF'
#IID	SEX
M1	1
M2	1
M3	1
F1	2
F2	2
F3	2
EOF

# Alternate .pvar that relabels the PAR1 variant back to chrom X at a PAR1
# coordinate (GRCh38 X:1000000). Reuses the same .pgen (par1 is stored diploid).
# This exercises the *coordinate*-based PAR detection and the `build` parameter:
# under GRCh38 X:1000000 is PAR (diploid); under build:='none' it is treated as
# chrX non-PAR (males haploid).
awk 'BEGIN{OFS="\t"} /^#/{print;next} $3=="par1"{$1="X"; print; next} {print}' \
    sexchr_example.pvar > sexchr_xpar.pvar

echo "Generated sexchr_example.{pgen,pvar,psam} and sexchr_xpar.pvar"
