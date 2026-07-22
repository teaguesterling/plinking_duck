#!/bin/bash
# Generate a RARE-variant, REF-major pfile fixture for the sparse (pgen difflist)
# sample-orient counts path (plinking_sample_counts_sparse).
#
# The committed AF=0.5 fixtures (pgen_example, wes_chr10, ...) are all
# genovec-encoded, so they never exercise the difflist fast path. This makes a
# fixture whose variants are RARE (so plink2 --make-pgen stores them as
# difflists) and REF-major (hom_ref is the common genotype, so common_geno==0 and
# the sparse fast path actually fires).
#
# Default output: rare_small.{pgen,pvar,psam} -- 256 samples x 400 ultra-rare
# variants, seed 777. Committed and used by test/sql/read_pfile_sample_counts_sparse.test.
#   ./generate_rare_test_data.sh                 # -> rare_small.* (the committed fixture)
#   ./generate_rare_test_data.sh N M MISS PREFIX # -> a custom (e.g. benchmark) fixture
#
# GOTCHA that cost real time: plink2 reads a plink1 .bed's A1 as the ALT allele
# (A2 = REF). So to get plink2 to report OUR alt-allele count g we must store
# g=0(hom_ref)->hom-A2=0b11, g=1(het)->0b10, g=2(hom_alt)->hom-A1=0b00, missing->0b01.
# The naive [0b00,0b10,0b11,...] makes ALT the major allele (hom_alt-major), which
# silently defeats the difflist fast path (common_geno becomes 2, never 0). ALWAYS
# verify sum(hom_ref) >> sum(hom_alt) on a supposedly-rare fixture.
#
# Requires: plink2 in PATH (or set PLINK2), python3 + numpy.
set -euo pipefail
PLINK2="${PLINK2:-plink2}"
N="${1:-256}"; M="${2:-400}"; MISS="${3:-0.01}"; PREFIX="${4:-rare_small}"
cd "$(dirname "$0")"

python3 - "$N" "$M" "$MISS" "$PREFIX" <<'PY'
import sys, numpy as np
n_samp, n_var, miss, prefix = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
rng = np.random.default_rng(777)
lo = 1.0 / (2 * n_samp)
# 95% rare (MAF in [1/2N, 1%]), 5% common -- exome-like carrier regime. The 1%
# cap keeps the range valid for small N (where 1/2N can exceed a lower cap) and
# still yields a ~4x sparse speedup at scale; ultra-rare (steeper cap) gives more.
is_common = rng.random(n_var) < 0.05
logp = np.where(is_common, rng.uniform(np.log10(0.01), np.log10(0.5), n_var),
                rng.uniform(np.log10(lo), np.log10(0.01), n_var))
pvec = 10.0 ** logp
# g = ALT count {0,1,2,3=missing} -> plink1 code, ALT=A1 (see header gotcha).
g2code = np.array([0b11, 0b10, 0b00, 0b01], dtype=np.uint8)
n_bytes = (n_samp + 3) // 4; pad = n_bytes * 4 - n_samp
with open(prefix + ".bed", "wb") as bed:
    bed.write(bytes([0x6c, 0x1b, 0x01]))
    for v in range(n_var):
        g = rng.binomial(2, pvec[v], n_samp).astype(np.uint8)
        if miss > 0: g[rng.random(n_samp) < miss] = 3
        codes = g2code[g]
        if pad: codes = np.concatenate([codes, np.zeros(pad, np.uint8)])
        c = codes.reshape(n_bytes, 4)
        bed.write((c[:,0] | (c[:,1]<<2) | (c[:,2]<<4) | (c[:,3]<<6)).astype(np.uint8).tobytes())
with open(prefix + ".bim", "w") as f:
    for v in range(n_var): f.write(f"1\tv{v}\t0\t{v+1}\tA\tG\n")
with open(prefix + ".fam", "w") as f:
    for i in range(n_samp): f.write(f"0\tS{i}\t0\t0\t0\t-9\n")
PY

"$PLINK2" --bed "$PREFIX.bed" --bim "$PREFIX.bim" --fam "$PREFIX.fam" --make-pgen --out "$PREFIX"
rm -f "$PREFIX.bed" "$PREFIX.bim" "$PREFIX.fam" "$PREFIX.log"
echo "Wrote $PREFIX.{pgen,pvar,psam} ($M variants x $N samples, rare + REF-major)."
