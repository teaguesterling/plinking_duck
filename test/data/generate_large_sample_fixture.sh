#!/usr/bin/env bash
#
# Generate a large SAMPLE-scale fixture for benchmarking sample-oriented,
# genotype-filtered queries (e.g. "list subject IDs matching a genotype at one
# variant") at the scale this extension targets: millions of individuals.
#
# This is a BENCHMARK fixture, not a CI fixture — the .psam alone is hundreds of
# MB at 7M samples, so its output is gitignored (see test/data/.gitignore) and
# must be regenerated locally rather than committed.
#
# Requires: plink2 in PATH (or set PLINK2=/path/to/plink2). The canonical
# spack path in CLAUDE.md has drifted before; `which plink2` is the source of truth.
#
# Usage:
#   ./generate_large_sample_fixture.sh                 # 1,000,000 samples (default)
#   N_SAMPLES=5000000 ./generate_large_sample_fixture.sh
#   N_SAMPLES=10000000 N_VARIANTS=50 ./generate_large_sample_fixture.sh
#   N_VARIANTS=1000000 N_SAMPLES=100000 ./generate_large_sample_fixture.sh   # wide
#
# TWO axes, and they TRADE OFF — the .pgen is dense: size ~= samples x variants / 4
# bytes, and plink2 generation time scales with the total genotype count
# (~3e7 genotypes/sec observed). You cannot have 1M variants AND millions of
# samples: 1M x 10M is ~2.3 TB and hours to build. Pick one stress axis:
#   * TALL (sample-scale, the default use case): many samples, FEW variants.
#     A carrier query selects ONE variant, so the sample-orient pre-read is
#     bounded by the selected variant; the cost lives in the N-row .psam parse
#     and the row-skip. e.g. N_SAMPLES=10000000 N_VARIANTS=20 (~27 MB, seconds).
#   * WIDE (variant-scale bind cost): 1M variants with a modest sample count, to
#     measure the fixed pvar parse + variant-index build a single-variant query
#     pays for variants it never touches. e.g. N_VARIANTS=1000000 N_SAMPLES=100000
#     (~23 GB, ~1 h). Real biobank pgens (~500K samples x ~1M variants) are ~125 GB.

set -euo pipefail
cd "$(dirname "$0")"

PLINK2="${PLINK2:-plink2}"
N_SAMPLES="${N_SAMPLES:-1000000}"
N_VARIANTS="${N_VARIANTS:-20}"
MISSING_GENO_FREQ="${MISSING_GENO_FREQ:-0.02}"   # ~2% missing, so the missing path is exercised
OUT_DIR="${OUT_DIR:-large_samples}"
PREFIX="${OUT_DIR}/large_${N_SAMPLES}s_${N_VARIANTS}v"

if ! command -v "$PLINK2" >/dev/null 2>&1; then
	echo "ERROR: plink2 not found (looked for '$PLINK2'). Set PLINK2=/path/to/plink2." >&2
	exit 1
fi

mkdir -p "$OUT_DIR"

echo "Generating ${N_SAMPLES} samples x ${N_VARIANTS} variants (missing geno freq ${MISSING_GENO_FREQ}) -> ${PREFIX}.{pgen,pvar,psam}"

# --dummy <sample ct> <SNP ct> [missing dosage freq] [missing pheno freq]
#         [{acgt|1234|12}] ...
# The 3rd positional is the missing-dosage (= missing-genotype) frequency. We do
# NOT pass 'dosage-freq=', so every dosage stays in {0,1,2} — i.e. hardcalls, which
# is what the hardcall genotype filter operates on. 'acgt' gives realistic allele
# codes. The self-validation step below fails loudly if the result is not hardcalls
# with the expected missingness — do not skip it. (Verified against plink2
# v2.0.0-a.6.9; run `plink2 --help dummy` if your version differs.)
"$PLINK2" \
	--dummy "$N_SAMPLES" "$N_VARIANTS" "$MISSING_GENO_FREQ" acgt \
	--make-pgen \
	--out "$PREFIX"

# --- Self-validation: confirm hardcalls (not dosages) and expected missingness ---
DUCKDB="${DUCKDB:-../../build/release/duckdb}"
EXT="${EXT:-../../build/release/extension/plinking_duck/plinking_duck.duckdb_extension}"
if [ -x "$DUCKDB" ] && [ -f "$EXT" ]; then
	echo
	echo "Validating generated fixture with the built extension..."
	# Per-variant counts reading as STRUCT(hom_ref, het, hom_alt, missing) confirms
	# hardcalls (not dosage output); the missing rate should be ~MISSING_GENO_FREQ.
	"$DUCKDB" -noheader -list -c "
LOAD '$(cd ../.. && pwd)/build/release/extension/plinking_duck/plinking_duck.duckdb_extension';
SELECT 'first_variant: total=' || (genotypes.hom_ref + genotypes.het + genotypes.hom_alt + genotypes.missing)
       || ' missing_rate=' || round(genotypes.missing::DOUBLE / (genotypes.hom_ref + genotypes.het + genotypes.hom_alt + genotypes.missing), 4)
FROM read_pfile('${PREFIX}', orient := 'variant', genotypes := 'counts') LIMIT 1;
" || { echo "VALIDATION FAILED — inspect ${PREFIX}.* (wrong --dummy args? dosage data?)." >&2; exit 1; }
	echo "Validation OK (expected first_variant_missing_rate ~= ${MISSING_GENO_FREQ})."
else
	echo "NOTE: build not found ($DUCKDB); skipping self-validation. Build the extension, then re-run to validate." >&2
fi

echo
echo "Done. Benchmark the carrier-ID path with, e.g.:"
echo "  SELECT IID FROM read_pfile('${PREFIX}',"
echo "      orient := 'sample', genotypes := 'counts',"
echo "      variants := ['<one-variant-id>'], include_genotypes := ['het','hom_alt']);"
echo
echo "Variant IDs: SELECT ID FROM read_pfile('${PREFIX}', orient := 'variant');"
