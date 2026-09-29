#!/usr/bin/env bash
set -euo pipefail
PART="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ROOT="$(cd "$PART/.." && pwd)"
PROFILE="${1:-stable}"
BLOCKS="${2:-50}"
OUT="${3:-$PART/results/$PROFILE}"
mkdir -p "$OUT"
"$PART/uy_reconstruction_ctct" --profile "$PROFILE" --blocks "$BLOCKS" --results-dir "$OUT" --disable-dynamic-scaling | tee "$OUT/run.log"
python3 "$ROOT/tools/plot_basic_results.py" "$OUT"
