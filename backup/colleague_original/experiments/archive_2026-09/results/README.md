# Included results

The repository includes selected raw CSV/log files so that the figures and summary values can be reproduced.

- `range_100/`: fixed-range 100-block ct-ct run.
- `mode_100/`: same-measurement controller-internal error diagnostic.
- `budget_sweep_100/`: L3/L5/L7 post-bootstrap tower sweep and summary CSV.
- `slot_50/`: 50-block bootstrap/control-representation diagnostic.
- `summary.csv`: compact top-level comparison generated from the included runs.

The figures in the root `figures/` directory are selected presentation-ready outputs. Raw diagnostic CSVs remain in each result folder.


- `mode_1000/`: requested 1000-block long-horizon run; stopped after 149 full blocks because CKKS decode reported excessive approximation error.
