#!/usr/bin/env bash
set -euo pipefail
PART="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ROOT="$(cd "$PART/.." && pwd)"
OPENFHE_ROOT="${OPENFHE_ROOT:-$HOME/openfhe-development}"
CXX="${CXX:-g++}"
OUT="$PART/uy_reconstruction_ptct"
"$CXX" -std=c++17 -O3 -DNDEBUG -DMATHBACKEND=4 -fopenmp \
  -I"$OPENFHE_ROOT/src/core/include" \
  -I"$OPENFHE_ROOT/src/pke/include" \
  -I"$OPENFHE_ROOT/src/binfhe/include" \
  -I"$OPENFHE_ROOT/build/src/core" \
  -I"$OPENFHE_ROOT/build/src/pke" \
  -I"$OPENFHE_ROOT/build/src/binfhe" \
  -I"$OPENFHE_ROOT/third-party/cereal/include" \
  -I"$ROOT/include" \
  "$PART/src/uy_reconstruction_ptct.cpp" \
  -L"$OPENFHE_ROOT/build/lib" -Wl,-rpath,"$OPENFHE_ROOT/build/lib" \
  -lOPENFHEpke -lOPENFHEcore -lOPENFHEbinfhe \
  -o "$OUT"
echo "built $OUT"
