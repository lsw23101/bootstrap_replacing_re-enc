#!/usr/bin/env bash
# usage: run.sh PROFILE K STEPS [OUT] [extra options, e.g. --d 1 --keep-imag]
set -euo pipefail
PART="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ROOT="$(cd "$PART/.." && pwd)"
PROFILE="${1:-stable}"
K="${2:-4}"
STEPS="${3:-200}"
OUT="${4:-$PART/results/${PROFILE}_K${K}}"
mkdir -p "$OUT"
"$PART/arx_window_ctct" --profile "$PROFILE" --K "$K" --steps "$STEPS" --results-dir "$OUT" "${@:5}" | tee "$OUT/run.log"
python3 "$ROOT/tools/plot_basic_results.py" "$OUT"
