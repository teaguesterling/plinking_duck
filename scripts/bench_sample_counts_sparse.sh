#!/bin/bash
# Benchmark: dense vs sparse (pgen difflist) sample-orient counts, and the
# effect of sample subsetting. Generates a large RARE, REF-major fixture and
# A/B-times `read_pfile(orient:='sample', genotypes:='counts')` under
# plinking_sample_counts_sparse = false/true, using DuckDB's own `.timer on`.
#
# NOT a pass/fail test (timings are machine-dependent) -- a reproducible harness
# to measure the sparse win on your own data/hardware. Correctness (sparse ==
# dense) is asserted by test/sql/read_pfile_sample_counts_sparse.test.
#
#   DUCKDB=./build/release/duckdb ./scripts/bench_sample_counts_sparse.sh [N_SAMP] [N_VAR]
#
# Requires: a built duckdb with the extension (DUCKDB=path), plink2 (PLINK2=path),
# python3 + numpy. Example result (100k x 30k, REF-major rare, this box):
#   dense  t1 ~4.5s  -> sparse t1 ~1.1s (~4x)     ultra-rare pushes it to ~6x
#   dense  t8 ~0.74s -> sparse t8 ~0.24s (~3x)
#   10% sample subset: ~5.8x on top (dense); 1%: ~15x (sublinear, fixed per-variant floor)
set -euo pipefail
DUCKDB="${DUCKDB:-./build/release/duckdb}"
PLINK2="${PLINK2:-plink2}"
N="${1:-100000}"; M="${2:-30000}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

echo "Generating $M variants x $N samples (rare, REF-major) ..."
PLINK2="$PLINK2" "$ROOT/test/data/generate_rare_test_data.sh" "$N" "$M" 0.005 "$TMP/bench" >/dev/null

P="$TMP/bench"
run() { # threads sparse
  "$DUCKDB" -c ".timer on" -c \
    "SET threads=$1; SET plinking_sample_counts_sparse=$2; SELECT sum(genotypes.het) FROM read_pfile('$P', orient:='sample', genotypes:='counts');" \
    2>&1 | grep -i "Run Time" | sed 's/^/    /'
}
for t in 1 8; do
  echo "threads=$t:"
  echo -n "  dense  "; run "$t" false
  echo -n "  sparse "; run "$t" true
done

echo "10% sample subset (dense, threads=8):"
"$DUCKDB" -c ".timer on" -c \
  "SET threads=8; SET VARIABLE s=(SELECT list('S'||(i*10)) FROM range(0,$((N/10))) t(i)); SELECT sum(genotypes.het) FROM read_pfile('$P', orient:='sample', genotypes:='counts', samples:=getvariable('s'));" \
  2>&1 | grep -i "Run Time" | sed 's/^/    /'
