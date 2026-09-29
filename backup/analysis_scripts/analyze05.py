import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/Users/nico/All/uy-reconstruction")
R = ROOT / "backup/05_results_archive"  # archived comparison runs
FIG = ROOT / "tex/figures"
runs = [("stable", 4), ("stable", 8), ("stable", 16), ("unstable_low", 4), ("unstable_low", 8),
        ("unstable_low", 16), ("unstable_high", 4), ("unstable_high", 8), ("unstable_high", 16)]
extra = [("unstable_low_K4_d1", "unstable_low", 4, 1), ("unstable_low_K4_keepimag", "unstable_low", 4, 4),
         ("unstable_low_K16_pack_500", "unstable_low", 16, 16), ("unstable_low_K4_packyu_500", "unstable_low", 4, 4)]

rows = []
def summarize(name, prof, K, d):
    d_ = R / name
    if not (d_ / "sample_trace.csv").exists():
        return None
    tr = pd.read_csv(d_ / "sample_trace.csv")
    tm = pd.read_csv(d_ / "timing.csv")
    log = (d_ / "run.log").read_text()
    ts = re.search(r"min_feasible_sampling_period_ms=([\d.e+-]+)", log)
    done = re.search(r"completed_steps=(\d+)", log)
    win = tr[tr.phase == "window"]
    eu = win.e_u_abs.to_numpy()
    q = max(len(eu) // 4, 1)
    pe = np.max(np.abs(tr[[f"plant_ckks_{i}" for i in range(4)]].to_numpy()
                       - tr[[f"plant_plain_{i}" for i in range(4)]].to_numpy()), axis=1)
    boots = tm[tm.bootstrap_ms > 0].bootstrap_ms
    steady = tm[tm.global_sample.isin(win.global_sample)].copy()
    steady["srv_ms"] = steady.u_eval_ms + (steady.register_ms if "register_ms" in steady else 0)
    tbt = boots.mean()
    bound = max(steady.srv_ms.max(), tbt * (1 if name.endswith("packyu") or "packyu_" in name else min(4, K)) / K, tbt / d)
    return dict(run=name, profile=prof, K=K, d=d, steps=int(done.group(1)) if done else len(tr),
                max_eu=tr.e_u_abs.max(), eu_first_q=eu[:q].mean(), eu_last_q=eu[-q:].mean(),
                max_plant=pe.max(), max_local=tr.e_local_abs.max(),
                boots_per_step=len(boots) / len(tr), boot_ms=tbt,
                u_ms_mean=steady.srv_ms.mean(), u_ms_max=steady.srv_ms.max(),
                mults_max=steady.ctct_mults.max(), rot_max=steady.rotations.max() if "rotations" in steady else 0,
                leak=tr.u_inactive_leakage.max() if "u_inactive_leakage" in tr else np.nan,
                Ts_model_ms=float(ts.group(1)) if ts else np.nan, Ts_steady_bound_ms=bound)

for prof, K in runs:
    for suf in ["", "_pack", "_packyu"]:
        r = summarize(f"{prof}_K{K}{suf}", prof, K, K)
        if r: rows.append(r)
for name, prof, K, d in extra:
    r = summarize(name, prof, K, d)
    if r: rows.append(r)
s = pd.DataFrame(rows)
s_ = s
pd.set_option("display.width", 250)
print(s.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
s.to_csv(R / "summary.csv", index=False)

# figure: u error vs step, one panel per profile, K = 4/8/16 (+ keep-imag on unstable_low)
S = ["#2a78d6", "#1baf7a", "#4a3aa7", "#eb6834"]
T2, GRID = "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": T2,
                     "xtick.color": T2, "ytick.color": T2})
fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), constrained_layout=True, sharey=True)
for ax, prof in zip(axes, ["stable", "unstable_low", "unstable_high"]):
    for c, K in zip(S, [4, 8, 16]):
        f = R / f"{prof}_K{K}_packyu/sample_trace.csv"
        if f.exists():
            tr = pd.read_csv(f)
            ax.plot(tr.global_sample, tr.e_u_abs.clip(lower=1e-18), color=c, lw=1.1, label=f"K=d={K}")
    if prof == "unstable_low":
        f = R / "unstable_low_K4_keepimag/sample_trace.csv"
        if f.exists():
            tr = pd.read_csv(f)
            ax.plot(tr.global_sample, tr.e_u_abs.clip(lower=1e-18), color=S[3], lw=1.1,
                    label="K=d=4, imag kept")
    ax.set_yscale("log"); ax.set_title(prof, loc="left", fontsize=10); ax.set_xlabel("step")
    ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
axes[0].set_ylabel(r"$|u_{ckks}-u_{plain}|$")
fig.suptitle("05 CT-CT delayed-window ARX, packed y + packed window, 200 steps", x=0.01, ha="left")
fig.savefig(FIG / "arx_packyu_u_error.png", dpi=160)

# figure: min feasible sampling period vs K for the three implementations
fig, axes = plt.subplots(1, 3, figsize=(13, 3.4), constrained_layout=True, sharey=True)
variants = [("", "scalar", S[3]), ("_pack", "packed y", S[1]), ("_packyu", "packed y + window", S[0])]
for ax, prof in zip(axes, ["stable", "unstable_low", "unstable_high"]):
    for suf, lab, c in variants:
        sub = s_[(s_.profile == prof) & (s_.run.str.fullmatch(rf"{prof}_K\d+{suf}"))].sort_values("K")
        if len(sub):
            ax.plot(sub.K, sub.Ts_model_ms, "o-", color=c, lw=2, ms=8, label=lab)
    ax.set_xscale("log", base=2); ax.set_xticks([4, 8, 16]); ax.set_xticklabels(["4", "8", "16"])
    ax.set_title(prof, loc="left", fontsize=10); ax.set_xlabel("K = d")
    ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
axes[0].set_ylabel("min feasible $T_s$ [ms]"); axes[0].legend(frameon=False, fontsize=8)
fig.suptitle("Real-time schedule model (N=4096): minimum sampling period", x=0.01, ha="left")
fig.savefig(FIG / "arx_min_Ts.png", dpi=160)
