#!/usr/bin/env python3
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import place_poles

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config.controller_profiles import A, B, C, XP0, get_profile, profile_metrics  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description="Plain recursive controller vs TAC reformulation")
    p.add_argument("--profile", default="stable", help="stable | unstable_low | unstable_high (hyphens also accepted)")
    p.add_argument("--samples", type=int, default=200)
    p.add_argument("--results-dir", default=None)
    p.add_argument("--tac-poles", nargs=4, type=float, default=[0.2, 0.3, 0.4, 0.5])
    return p.parse_args()


def simulate(profile, samples: int, tac_poles):
    F, G, H = profile.F, profile.G, profile.H
    # Dual pole placement: F - R H has the requested Schur poles.
    R = place_poles(F.T, H.reshape(-1, 1), tac_poles).gain_matrix.T.reshape(-1)
    Ftac = F - np.outer(R, H)

    xp_ref = XP0.copy(); xc_ref = profile.xc0.copy()
    xp_tac = XP0.copy(); xc_tac = profile.xc0.copy()

    rows = []
    for k in range(samples):
        y_ref = float(C @ xp_ref)
        u_ref = float(H @ xc_ref)
        y_tac = float(C @ xp_tac)
        u_tac = float(H @ xc_tac)

        rows.append({
            "sample": k,
            "u_reference": u_ref,
            "u_test": u_tac,
            "u_error": u_tac - u_ref,
            **{f"xp_reference_{i}": xp_ref[i] for i in range(4)},
            **{f"xp_test_{i}": xp_tac[i] for i in range(4)},
            "plant_error_inf": float(np.linalg.norm(xp_tac - xp_ref, np.inf)),
        })

        xp_ref_next = A @ xp_ref + B * u_ref
        xc_ref_next = F @ xc_ref + G * y_ref

        # TAC identity: x+ = (F-RH)x + Gy + Ru.
        xp_tac_next = A @ xp_tac + B * u_tac
        xc_tac_next = Ftac @ xc_tac + G * y_tac + R * u_tac

        xp_ref, xc_ref = xp_ref_next, xc_ref_next
        xp_tac, xc_tac = xp_tac_next, xc_tac_next

    return rows, R, Ftac


def save_results(rows, out: Path, profile_name: str):
    out.mkdir(parents=True, exist_ok=True)
    with (out / "sample_trace.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    k = np.array([r["sample"] for r in rows])
    u_ref = np.array([r["u_reference"] for r in rows]); u = np.array([r["u_test"] for r in rows])
    eu = np.abs(u-u_ref)
    ep = np.array([r["plant_error_inf"] for r in rows])
    xp_ref = np.array([[r[f"xp_reference_{i}"] for i in range(4)] for r in rows])
    xp = np.array([[r[f"xp_test_{i}"] for i in range(4)] for r in rows])

    plt.figure(); plt.plot(k, u_ref, label="reference"); plt.plot(k, u, "--", label="TAC"); plt.xlabel("sample"); plt.ylabel("u"); plt.legend(); plt.tight_layout(); plt.savefig(out/"u.png", dpi=160); plt.close()
    plt.figure();
    for i in range(4): plt.plot(k, xp_ref[:,i], label=f"ref x_p[{i}]"); plt.plot(k, xp[:,i], "--", label=f"TAC x_p[{i}]")
    plt.xlabel("sample"); plt.ylabel("plant state"); plt.legend(ncol=2, fontsize=8); plt.tight_layout(); plt.savefig(out/"plant_state.png", dpi=160); plt.close()
    plt.figure(); plt.semilogy(k, np.maximum(eu, 1e-18)); plt.xlabel("sample"); plt.ylabel("|u error|"); plt.tight_layout(); plt.savefig(out/"u_error.png", dpi=160); plt.close()
    plt.figure(); plt.semilogy(k, np.maximum(ep, 1e-18)); plt.xlabel("sample"); plt.ylabel("plant error inf-norm"); plt.tight_layout(); plt.savefig(out/"plant_state_error.png", dpi=160); plt.close()


def main():
    args = parse_args(); profile = get_profile(args.profile)
    out = Path(args.results_dir or (Path(__file__).resolve().parent / "results" / args.profile))
    rows, R, Ftac = simulate(profile, args.samples, args.tac_poles)
    save_results(rows, out, args.profile)
    m = profile_metrics(profile)
    print(f"profile={profile.name}: {profile.label}")
    print(f"rho(F)={m['rho_F']:.9f}, rho(Acl)={m['rho_Acl']:.9f}")
    print(f"TAC poles={np.linalg.eigvals(Ftac)}")
    print(f"R={R}")
    print(f"max |u_error|={max(abs(r['u_error']) for r in rows):.3e}")
    print(f"max plant_error_inf={max(r['plant_error_inf'] for r in rows):.3e}")
    print(f"saved: {out}")

if __name__ == "__main__": main()
