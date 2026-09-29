#!/usr/bin/env bash
set -euo pipefail
PART="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ROOT=/Users/nico/All/uy-reconstruction
OPENFHE_ROOT="${OPENFHE_ROOT:-$HOME/openfhe-development}"
CXX="${CXX:-g++}"
OUT=/private/tmp/claude-501/-Users-nico-All-uy-reconstruction/674bcf46-31fe-490a-bfd7-5fd224fa8ad4/scratchpad/ptct_imag
"$CXX" -std=c++17 -O3 -DNDEBUG -DMATHBACKEND=4 -fopenmp \
  -I"$OPENFHE_ROOT/src/core/include" \
  -I"$OPENFHE_ROOT/src/pke/include" \
  -I"$OPENFHE_ROOT/src/binfhe/include" \
  -I"$OPENFHE_ROOT/build/src/core" \
  -I"$OPENFHE_ROOT/build/src/pke" \
  -I"$OPENFHE_ROOT/build/src/binfhe" \
  -I"$OPENFHE_ROOT/third-party/cereal/include" \
  -I"$ROOT/include" -I"$PART/src" \
  "/private/tmp/claude-501/-Users-nico-All-uy-reconstruction/674bcf46-31fe-490a-bfd7-5fd224fa8ad4/scratchpad/ptct_imag.cpp" \
  -L"$OPENFHE_ROOT/build/lib" -Wl,-rpath,"$OPENFHE_ROOT/build/lib" \
  -lOPENFHEpke -lOPENFHEcore -lOPENFHEbinfhe \
  -o "$OUT"
echo "built $OUT"
