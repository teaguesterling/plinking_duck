#!/bin/bash
# Generate variant-sharded pfile fixtures for read_pfile() LIST(VARCHAR) tests
# (P8-001 multi-file reads).
#
# Splits the committed large_example (3000 variants x 8 samples, chroms 1/2/3)
# into 3 variant-disjoint shards via round-robin --extract on variant ID:
#   shard1 <- var1, var4, var7, ...   (index % 3 == 0)
#   shard2 <- var2, var5, var8, ...   (index % 3 == 1)
#   shard3 <- var3, var6, var9, ...   (index % 3 == 2)
#
# Round-robin (not per-chromosome) is deliberate: every chromosome is spread
# across all 3 shards, so a `region := 'chr:start-end'` filter genuinely must
# pull variants from multiple shards -- exercising the cross-source region path.
#
# All three shards keep ALL 8 samples, so their .psam is identical (same IIDs,
# same order) -- the positional-alignment contract read_pfile([...]) relies on.
# Reading read_pfile([shard1,shard2,shard3]) must reproduce the SAME variants and
# genotypes as reading the whole large_example.
#
# Requires: plink2 in PATH (or set PLINK2 environment variable).
set -euo pipefail
cd "$(dirname "$0")"

PLINK2="${PLINK2:-plink2}"

if [[ ! -f large_example.pvar ]]; then
    echo "ERROR: large_example.pvar not found; run generate_test_data.sh first" >&2
    exit 1
fi

# Build round-robin variant-ID extract lists from large_example.pvar (skip header lines).
rm -f shard1.extract shard2.extract shard3.extract
i=0
while IFS=$'\t' read -r chrom pos id rest; do
    case $chrom in \#*) continue;; esac
    shard=$(( (i % 3) + 1 ))
    echo "$id" >> "shard${shard}.extract"
    i=$(( i + 1 ))
done < large_example.pvar

for s in 1 2 3; do
    "$PLINK2" --pfile large_example --extract "shard${s}.extract" \
              --make-pgen --out "shard${s}" --allow-extra-chr
done

rm -f shard1.extract shard2.extract shard3.extract
rm -f shard1.log shard2.log shard3.log

echo "Shard fixtures generated:"
for s in 1 2 3; do
    n=$(grep -vc '^#' "shard${s}.pvar")
    echo "  shard${s}: ${n} variants"
    ls -la "shard${s}.pgen" "shard${s}.pvar" "shard${s}.psam"
done
