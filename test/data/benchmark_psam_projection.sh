#!/usr/bin/env bash
#
# Benchmark the carrier-lookup / sample-orient path against a large fixture, and
# show that a WIDE psam (many covariate columns) costs ~the same as a NARROW one
# thanks to the projection-aware parquet psam load (v0.6.1+). This is a manual
# benchmark helper, NOT a CI test — timings are machine-dependent.
#
# Prereqs:
#   1. Build the extension:  make -j 4
#   2. A large fixture:      N_SAMPLES=10000000 N_VARIANTS=20 ./generate_large_sample_fixture.sh
#
# Usage:
#   ./benchmark_psam_projection.sh [fixture_prefix]
#   ./benchmark_psam_projection.sh large_samples/large_10000000s_20v      # default
#   REPEAT=5 ./benchmark_psam_projection.sh large_samples/large_5000000s_20v
#
# It builds a narrow (IID/FID/SEX) and a wide (+ N_COVARIATES DOUBLE cols) parquet
# psam companion from the fixture's .psam, then times a few representative queries.

set -euo pipefail
cd "$(dirname "$0")"
ROOT="$(cd ../.. && pwd)"

DUCKDB="${DUCKDB:-$ROOT/build/release/duckdb}"
EXT="${EXT:-$ROOT/build/release/extension/plinking_duck/plinking_duck.duckdb_extension}"
FIXTURE="${1:-large_samples/large_10000000s_20v}"
REPEAT="${REPEAT:-3}"
N_COVARIATES="${N_COVARIATES:-20}"
TMP="${TMPDIR:-/tmp}"

if [ ! -x "$DUCKDB" ] || [ ! -f "$EXT" ]; then
	echo "ERROR: build not found. Run 'make -j 4' first (looked for $DUCKDB / $EXT)." >&2
	exit 1
fi
if [ ! -f "$FIXTURE.pgen" ]; then
	echo "ERROR: fixture '$FIXTURE.pgen' not found. Generate one with generate_large_sample_fixture.sh." >&2
	exit 1
fi

LOAD="LOAD '$EXT';"
NARROW="$TMP/bench_narrow.psam.parquet"
WIDE="$TMP/bench_wide.psam.parquet"

# Narrow companion = the psam as-is; wide = same + N_COVARIATES DOUBLE columns.
covs=""
for i in $(seq 0 $((N_COVARIATES - 1))); do covs="$covs, random() AS c$i"; done
echo "Building parquet psam companions (narrow + wide with $N_COVARIATES covariates)..."
"$DUCKDB" -c "$LOAD COPY (SELECT * FROM read_psam('$FIXTURE.psam')) TO '$NARROW' (FORMAT parquet);" >/dev/null
"$DUCKDB" -c "$LOAD COPY (SELECT *$covs FROM read_psam('$FIXTURE.psam')) TO '$WIDE' (FORMAT parquet);" >/dev/null

n_samples=$("$DUCKDB" -noheader -list -c "$LOAD SELECT COUNT(*) FROM read_psam('$FIXTURE.psam');")
v0=$("$DUCKDB" -noheader -list -c "$LOAD SELECT ID FROM read_pfile('$FIXTURE', orient := 'variant') LIMIT 1;")

# best-of-REPEAT wall seconds for one SQL statement
bench() {
	local sql="$1" best=99 t
	for _ in $(seq 1 "$REPEAT"); do
		t=$( { /usr/bin/time -f "%e" "$DUCKDB" -noheader -list -c "$LOAD $sql" >/dev/null; } 2>&1 | tail -1 )
		awk "BEGIN{exit !($t < $best)}" && best=$t
	done
	echo "$best"
}

# __PSAM__ is substituted with the narrow/wide companion path.
declare -a LABELS=(
	"carrier SELECT IID (het/hom_alt)"
	"SELECT IID all samples (no filter)"
	"counts-only, NO psam column projected"
	"IID + SEX (2 psam cols), 1 row"
)
declare -a QUERIES=(
	"SELECT count(*) FROM (SELECT IID FROM read_pfile('$FIXTURE', orient:='sample', genotypes:='counts', variants:=['$v0'], include_genotypes:=['het','hom_alt'], psam:='__PSAM__'));"
	"SELECT count(*) FROM (SELECT IID FROM read_pfile('$FIXTURE', orient:='sample', genotypes:='counts', variants:=['$v0'], psam:='__PSAM__'));"
	"SELECT genotypes FROM read_pfile('$FIXTURE', orient:='sample', genotypes:='counts', variants:=['$v0'], psam:='__PSAM__') LIMIT 1;"
	"SELECT IID, SEX FROM read_pfile('$FIXTURE', orient:='sample', genotypes:='counts', variants:=['$v0'], psam:='__PSAM__') LIMIT 1;"
)

echo
echo "Fixture: $FIXTURE  ($n_samples samples, variant '$v0')   best of $REPEAT, wall seconds"
printf '%-46s %10s %14s\n' "query" "narrow(3)" "wide($((3 + N_COVARIATES)))"
printf '%-46s %10s %14s\n' "----------------------------------------------" "---------" "-------------"
for i in "${!QUERIES[@]}"; do
	qn="${QUERIES[$i]//__PSAM__/$NARROW}"
	qw="${QUERIES[$i]//__PSAM__/$WIDE}"
	printf '%-46s %10s %14s\n' "${LABELS[$i]}" "$(bench "$qn")" "$(bench "$qw")"
done
echo
echo "Projection-aware load: the wide column reflects only the columns each query"
echo "selects, so extra covariate columns add ~no cost."
rm -f "$NARROW" "$WIDE"
