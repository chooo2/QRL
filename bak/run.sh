#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

MODE="${1:-run}"
[ "$#" -gt 0 ] && shift

case "$MODE" in
  hpo)
    FAMILIES="${HPO_FAMILIES:-grid xtree falcon eagle aspen_11 aspen_m}"
    for f in $FAMILIES; do
      echo "[HPO] ===== tuning config: $f ====="
      python src/bayesopt.py --family "$f" "$@"
      python scripts/apply_hpo_result.py "$f"
    done
    echo "[HPO] running src/main.py with per-family tuned TrainingConfig..."
    python src/main.py
    ;;
  run)
    echo "Running src/main.py..."
    python src/main.py "$@"
    ;;
  *)
    echo "Usage: $0 {hpo [bayesopt.py args...]|run [main.py args...]}" >&2
    exit 1
    ;;
esac
