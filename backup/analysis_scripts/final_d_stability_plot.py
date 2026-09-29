"""Figure for tex: 128-bit, K=64, per-step |u_ckks - u_plain| for d=64 (diverges) vs d=61 (bounded)."""
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "05_openfhe_ctct_arx_window/results"
ARCH = ROOT / "backup/05_results_archive"
S1, S2, S3, T2, GRID = "#2a78d6", "#1baf7a", "#eb6834", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": T2,
                     "xtick.color": T2, "ytick.color": T2})
fig, ax = plt.subplots(figsize=(10, 3.6), constrained_layout=True)
runs = [("unstable_high_K64_s128_1000", "unstable_high, d=64 (d > K-n+1)", S3),
        ("unstable_high_K64_d61_s128", "unstable_high, d=61", S1),
        ("unstable_low_K64_d61_s128", "unstable_low, d=61", S2)]
for name, lab, c in runs:
    f = R / name / "sample_trace.csv"
    if not f.exists():
        f = ARCH / name / "sample_trace.csv"
    d = pd.read_csv(f)
    ax.plot(d.global_sample, d.e_u_abs.clip(lower=1e-18), color=c, lw=1.1, label=lab)
for b in range(64, 1000, 64):
    ax.axvline(b, color=GRID, lw=0.8, zorder=0)
ax.set_yscale("log"); ax.set_xlabel("step (vertical lines: block starts, K=64)")
ax.set_ylabel(r"$|u_{ckks}-u_{plain}|$"); ax.grid(axis="y", color=GRID, lw=0.6)
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.set_title("128-bit (N=65536), K=64, packed y + window", loc="left", fontsize=10)
fig.savefig(ROOT / "tex/figures/final_d_stability.png", dpi=160)
