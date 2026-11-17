#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
run_dir="${1:-runs/synthetic_demo}"
python examples/make_synthetic_input.py "$run_dir/input"
python -m nmr_cohort run \
  --bin "$run_dir/input/spectra.tsv" --format header \
  --batch-map "$run_dir/input/batch_map.tsv" --groups "$run_dir/input/groups.json" \
  --reference-stats "$run_dir/input/reference_stats.tsv" \
  --out "$run_dir/result"
python scripts/check_run.py "$run_dir/result/quantified" "$run_dir/result/filtered"
