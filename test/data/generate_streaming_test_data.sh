#!/bin/bash
# Generate large PLINK2 test data for streaming and threading tests
# Requires: plink2 binary
# Creates streaming_example.{pgen,pvar,psam} with 50K variants × 8 samples
set -euo pipefail
cd "$(dirname "$0")"

PLINK2="${PLINK2:-/mnt/aux-data/teague/Dev/spack/opt/spack/linux-zen3/plink2-2.0.0-a.6.9-e2l3wx22vy6tkmdwzhaexlml4rs3okx3/bin/plink2}"

echo "Generating streaming test data (50K variants × 8 samples)..."

# Generate VCF with 50K variants across 3 chromosomes
python3 -c "
import random
random.seed(42)
print('##fileformat=VCFv4.3')
for c in [1, 2, 3]:
    print(f'##contig=<ID={c},length=100000000>')
print('##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">')
samples = '\t'.join(f'SAMPLE{i}' for i in range(1, 9))
print(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}')
vid = 0
for chrom in [1, 2, 3]:
    n_vars = [20000, 15000, 15000][chrom-1]
    for i in range(n_vars):
        vid += 1
        pos = (i + 1) * 100
        ref, alt = random.choice([('A','G'), ('C','T'), ('G','A'), ('T','C')])
        gts = []
        for _ in range(8):
            r = random.random()
            if r < 0.05: gts.append('./.')
            elif r < 0.30: gts.append('0/0')
            elif r < 0.70: gts.append('0/1')
            else: gts.append('1/1')
        gt_str = '\t'.join(gts)
        print(f'{chrom}\t{pos}\tvar{vid}\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt_str}')
" > streaming_example.vcf

"$PLINK2" --vcf streaming_example.vcf --make-pgen --out streaming_example
rm -f streaming_example.vcf streaming_example.log

echo "Generated files:"
ls -la streaming_example.pgen streaming_example.pvar streaming_example.psam
echo "Done."
