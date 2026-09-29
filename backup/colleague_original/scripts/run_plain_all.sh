#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
for PROFILE in stable unstable_low unstable_high; do
  echo "===== 01 TAC: $PROFILE ====="
  python3 "$ROOT/01_plain_tac/tac_plain.py" --profile "$PROFILE" --samples 200
  echo "===== 02 u/y reconstruction: $PROFILE ====="
  python3 "$ROOT/02_plain_uy_reconstruction/uy_reconstruction_plain.py" --profile "$PROFILE" --blocks 50 --period 4 --nu 4
done
