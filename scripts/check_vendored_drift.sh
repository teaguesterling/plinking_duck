#!/usr/bin/env bash
# check_vendored_drift.sh — drift canary for verbatim-vendored plink2 extracts.
#
# Two source files in this repo are hand-copied ("vendored") extracts of upstream
# plink-ng (plink2) functions, because their home translation units are huge and
# entangled with plink2's g_bigstack arena / ThreadGroup / CLI model and cannot be
# linked wholesale (see docs/planning/plink2-wrapping-audit.md, YELLOW #2/#3):
#
#   src/plink2_glm_logistic_math.cpp  <- third_party/plink-ng/2.0/plink2_glm_logistic.cc
#   src/vcf_genotype_parse.cpp        <- third_party/plink-ng/2.0/plink2_import.cc
#
# The plink-ng submodule is pinned to a fixed commit, so upstream only moves when
# someone bumps the submodule. This canary exists to catch exactly that: a submodule
# bump that alters any vendored function body. When that happens the vendored copies
# must be re-extracted and the pin comments in the .cpp headers updated.
#
# HOW IT DETECTS DRIFT (function-marker extraction, not raw line numbers):
#   For each named vendored function, the script locates the *definition* in both the
#   vendored copy and the CURRENT submodule source by a column-0 "<ret> Name(" marker
#   (call sites are indented, so column-0 anchoring skips them), then captures the body
#   from the first '{' to the terminating column-0 '}' (top-level closing brace; every
#   nested brace is indented). Signatures are intentionally dropped: the vendored copies
#   are clang-format-reflowed and some gain a 'static' qualifier, none of which changes
#   behavior. Bodies are normalized (strip /* */ and // comments, then strip ALL
#   whitespace) and sha256'd. A function may have multiple column-0 definitions guarded
#   by #ifdef __LP64__ / #else (SIMD vs scalar); all occurrences are compared in order,
#   and an occurrence-count mismatch is itself a failure.
#
#   The GLM extract's bodies are byte-identical to upstream. The VCF extract applies two
#   documented renamings (VcfHalfCall type -> uint32_t, kVcfHalfCall* constants -> kHalfCall*,
#   via a locally-defined enum in src/include/vcf_genotype_parse.hpp); the script applies the
#   same renamings to the upstream side (VCF_UPSTREAM_RENAMES) so the comparison stays an
#   exact body match. Any drift beyond those two known renamings is reported.
#
# LIMITATIONS (documented on purpose):
#   - Signature-only upstream changes (param type/name, added/removed qualifier) are NOT
#     caught — only bodies are compared. Body changes, added/removed occurrences, and
#     outright removal of a vendored function ARE caught.
#   - The function name lists below are the source of truth for "what we vendored". A
#     brand-new upstream function we never copied is (correctly) out of scope.
#   - Only function-body TEXT is compared. The vendored bodies also depend on macros and
#     inline helpers resolved from live submodule headers at build time (e.g. PAIR_TABLE16,
#     kBitsPerWordD2, AdvToNthDelimChecked). If a submodule bump changes one of THOSE
#     definitions, behavior can drift while the compared body text stays identical and this
#     canary still prints OK. Header-level changes are out of scope by design.
#
# USAGE:  scripts/check_vendored_drift.sh
#   exit 0 + "OK: <file> matches pinned upstream"  when every body is byte-identical
#   exit 1 + a message naming the file + re-extract instructions  on any drift
#
# CI: this is a clean, dependency-light check (bash + awk + perl + sha256sum). To wire
# it into .github/workflows/MainDistributionPipeline.yml, add a step BEFORE the build,
# e.g. in a job that has checked out submodules:
#     - name: Check vendored plink2 drift
#       run: bash scripts/check_vendored_drift.sh
# (Left as documentation rather than auto-added to avoid perturbing the release matrix.)

set -u

# Resolve repo root from this script's location so it runs from anywhere.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT" || { echo "cannot cd to repo root"; exit 2; }

UPSTREAM_GLM="third_party/plink-ng/2.0/plink2_glm_logistic.cc"
UPSTREAM_VCF="third_party/plink-ng/2.0/plink2_import.cc"
VENDORED_GLM="src/plink2_glm_logistic_math.cpp"
VENDORED_VCF="src/vcf_genotype_parse.cpp"

# Functions vendored from plink2_glm_logistic.cc (the entire L147-1012 solver block,
# minus LogisticRegressionResidualizedF which we do not use).
GLM_FUNCS=(
  GenoarrToFloatsRemoveMissing
  MultMatrixDxnVectNF
  LogisticSseF
  ComputeVAndPMinusYF
  ComputeVF
  TripleProductF
  ComputeTwoDiagTripleProductF
  ComputeThreeTripleProductF
  ComputeTwoPlusOneTripleProductF
  CopyAndMeanCenterF
  ComputeLoglikCheckedF
  ComputeHessianF
  CholeskyDecompositionF
  SolveLinearSystemF
  LogisticRegressionF
  FirthComputeHdiagWeightsF
  FirthComputeSecondWeightsF
  FirthRegressionF
)

# Functions vendored from plink2_import.cc (the biallelic GT parsers + qual check).
# NOTE: the other functions in vcf_genotype_parse.cpp (FillBaseContext, TranslateResult,
# ParseUnphasedBiallelicGT, ParsePhasedBiallelicGT) are PlinkingDuck glue, not upstream
# copies, and are deliberately absent from this list.
VCF_FUNCS=(
  VcfCheckQuals
  VcfConvertUnphasedBiallelicLine
  VcfConvertPhasedBiallelicLine
)

# Documented adaptations applied to the VENDORED copies (see src/vcf_genotype_parse.cpp
# and src/include/vcf_genotype_parse.hpp): the extract uses a locally-defined kHalfCall*
# enum and a plain uint32_t in place of upstream's VcfHalfCall type/constants. These sed
# rules are applied to the UPSTREAM side so the two are compared on equal footing; they
# are pure renamings, not logic. Order matters (kVcfHalfCall before the bare type name).
# GLM extract has no such adaptations (bodies are byte-identical) -> empty rule set.
VCF_UPSTREAM_RENAMES='s/kVcfHalfCall/kHalfCall/g; s/VcfHalfCall/uint32_t/g'

# Emit the normalized sha256 of each column-0 definition body of $name in $file,
# one hash per line, in file order. Prints nothing if the function is absent.
# $3 (optional) is an extra sed program applied after whitespace stripping (rename rules).
body_hashes() {
  local name="$1" file="$2" extra_sed="${3:-}"
  awk -v name="$name" '
    BEGIN { re = "^[A-Za-z_].*[ *]" name "\\(|^" name "\\(" }
    state==0 { if ($0 ~ re) { state=1; buf=$0 "\n" } next }
    state==1 {
      buf = buf $0 "\n"
      if ($0 ~ /^\}/) {
        i = index(buf, "{")
        printf "%s\x01", substr(buf, i+1)   # body only (drop signature up to first {)
        state=0
      }
    }
  ' "$file" \
  | perl -0777 -ne '
      for my $body (split /\x01/, $_) {
        next unless $body =~ /\S/;
        $body =~ s{/\*.*?\*/}{}gs;   # block comments
        $body =~ s{//[^\n]*}{}g;     # line comments
        $body =~ s/\s+//g;           # all whitespace
        print "$body\n";
      }
    ' \
  | { if [[ -n "$extra_sed" ]]; then sed "$extra_sed"; else cat; fi; } \
  | while IFS= read -r line; do printf '%s' "$line" | sha256sum | cut -d" " -f1; done
}

check_file() {
  local label="$1" vendored="$2" upstream="$3" upstream_renames="$4"; shift 4
  local funcs=("$@")
  local file_ok=1

  if [[ ! -f "$vendored" ]]; then echo "MISSING vendored file: $vendored"; return 1; fi
  if [[ ! -f "$upstream" ]]; then
    echo "MISSING upstream file: $upstream (is the plink-ng submodule checked out?)"; return 1; fi

  for fn in "${funcs[@]}"; do
    local vh uh
    vh="$(body_hashes "$fn" "$vendored")"
    uh="$(body_hashes "$fn" "$upstream" "$upstream_renames")"
    if [[ -z "$vh" ]]; then
      echo "DRIFT: $vendored — vendored function '$fn' not found (marker changed?)"; file_ok=0; continue
    fi
    if [[ -z "$uh" ]]; then
      echo "DRIFT: upstream $upstream no longer defines '$fn' (vendored in $vendored)"; file_ok=0; continue
    fi
    local vc uc
    vc="$(printf '%s\n' "$vh" | grep -c .)"
    uc="$(printf '%s\n' "$uh" | grep -c .)"
    if [[ "$vc" != "$uc" ]]; then
      echo "DRIFT: '$fn' definition count differs (vendored=$vc upstream=$uc)"; file_ok=0; continue
    fi
    if [[ "$vh" != "$uh" ]]; then
      echo "DRIFT: '$fn' body differs between $vendored and $upstream"; file_ok=0
    fi
  done

  if [[ "$file_ok" == 1 ]]; then
    echo "OK: $vendored matches pinned upstream ($upstream)"
    return 0
  fi
  echo "----"
  echo "FAIL: $label — vendored/upstream mismatch (expected: submodule bumped); re-extract the listed function(s)"
  echo "      from $upstream and update the pin comment (upstream commit + line range)"
  echo "      at the top of $vendored."
  return 1
}

echo "Vendored-drift canary — pinned plink-ng commit: $(git -C third_party/plink-ng rev-parse --short HEAD 2>/dev/null || echo '?')"
rc=0
check_file "GLM logistic/Firth solver" "$VENDORED_GLM" "$UPSTREAM_GLM" "" "${GLM_FUNCS[@]}" || rc=1
check_file "VCF biallelic GT parser"   "$VENDORED_VCF" "$UPSTREAM_VCF" "$VCF_UPSTREAM_RENAMES" "${VCF_FUNCS[@]}" || rc=1
exit $rc
