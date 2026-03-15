#!/bin/bash
# Generate PCA test dataset (250 samples x 500 variants)
# Requires plink2 in PATH or set PLINK2 environment variable

set -euo pipefail

PLINK2="${PLINK2:-plink2}"
DIR="$(cd "$(dirname "$0")" && pwd)"

echo "Generating pca_example (250 samples x 500 variants)..."
"$PLINK2" --dummy 250 500 --make-pgen --out "${DIR}/pca_example" 2>&1
rm -f "${DIR}/pca_example-temporary."* "${DIR}/pca_example.log"

echo "Done. Files created:"
ls -la "${DIR}"/pca_example.{pgen,pvar,psam}
