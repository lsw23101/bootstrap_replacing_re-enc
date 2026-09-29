import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/Users/nico/All/uy-reconstruction")
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT
PARTS = [
    ("03 PT-CT", "03_openfhe_ptct_uy_reconstruction", "epsilon_u_abs"),
    ("04 CT-CT", "04_openfhe_ctct_uy_reconstruction", "e_u_nominal_abs"),
]
PROFILES = ["stable", "unstable_low", "unstable_high"]
S1, S2 = "#2a78d6", "#eb6834"
TEXT2, GRID = "#52514e", "#e4e3df"

plt.rcParams.update({"font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": TEXT2,
                     "xtick.color": TEXT2, "ytick.color": TEXT2})
fig, axes = plt.subplots(2, 3, figsize=(13, 6.5), constrained_layout=True)
rows = []
for r, (label, part, col) in enumerate(PARTS):
    for c, prof in enumerate(PROFILES):
        ax = axes[r, c]
        f = ROOT / part / "results" / f"{prof}_500" / "sample_trace.csv"
        if not f.exists():
            ax.set_title(f"{label} / {prof}: missing"); continue
        d = pd.read_csv(f)
        t = d.global_sample.to_numpy(); e = d[col].abs().to_numpy()
        # skip block 0 (fresh encryption, no bootstrap yet)
        m = d.block.to_numpy() > 0
        slope, ic = np.polyfit(t[m], np.log(e[m]), 1)
        q = len(e[m]) // 4
        early, late = e[m][:q].mean(), e[m][-q:].mean()
        rho = pd.Series(t[m]).corr(pd.Series(e[m]), method="spearman")
        rows.append(dict(part=label, profile=prof, steps=len(d), max=e.max(),
                         mean_first_quarter=early, mean_last_quarter=late,
                         late_over_early=late / early, growth_per_block=np.exp(slope*4),
                         
                         spearman=rho))
        ax.plot(t, e, color=S1, lw=1, label="|u_ckks − u_plain|")
        pass
        ax.set_title(f"{label} / {prof}", loc="left", fontsize=10)
        ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
        ax.set_yscale("log")
        if r == 1: ax.set_xlabel("step")
        if c == 0: ax.set_ylabel("|u error|")
axes[0, 0].legend(frameon=False, fontsize=8)
fig.suptitle("u error vs step, 500 steps (K=4, 125 blocks)", x=0.01, ha="left")
fig.savefig(OUT / "u_error_500.png", dpi=150)
s = pd.DataFrame(rows)
s.to_csv(OUT / "u_error_500_summary.csv", index=False)
pd.set_option("display.width", 200)
print(s.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
