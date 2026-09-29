#!/usr/bin/env python3
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config.controller_profiles import (A, B, C, XP0, get_profile, profile_metrics,
                                        reconstruction_matrices)  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description="Plain periodic u/y state reconstruction")
    p.add_argument("--profile", default="stable", help="stable | unstable_low | unstable_high (hyphens also accepted)")
    p.add_argument("--blocks", type=int, default=50)
    p.add_argument("--period", type=int, default=4, help="reconstruction period K")
    p.add_argument("--nu", type=int, default=4, help="reconstruction horizon nu (current 4-state example: 4)")
    p.add_argument("--results-dir", default=None)
    return p.parse_args()


def simulate(profile, blocks: int, K: int, nu: int):
    if K < nu: raise ValueError("period K must satisfy K >= nu")
    F,G,H = profile.F, profile.G, profile.H
    O,T,Mu,My = reconstruction_matrices(F,G,H,nu)

    xp_ref=XP0.copy(); xc_ref=profile.xc0.copy()
    xp=XP0.copy(); xc_block=profile.xc0.copy()
    rows=[]
    for b in range(blocks):
        ys=[]; us=[]
        for i in range(K):
            # Reference recursive controller.
            y_ref=float(C@xp_ref); u_ref=float(H@xc_ref)

            # Lifted controller evaluated from the block-start reconstructed state.
            y=float(C@xp)
            u=float(H @ np.linalg.matrix_power(F,i) @ xc_block)
            for j in range(i):
                u += float(H @ np.linalg.matrix_power(F,i-1-j) @ G) * ys[j]
            ys.append(y); us.append(u)

            rows.append({
                "block":b,"within_block":i,"sample":b*K+i,
                "u_reference":u_ref,"u_test":u,"u_error":u-u_ref,
                **{f"xp_reference_{r}":xp_ref[r] for r in range(4)},
                **{f"xp_test_{r}":xp[r] for r in range(4)},
                "plant_error_inf":float(np.linalg.norm(xp-xp_ref,np.inf)),
            })

            xp_ref=A@xp_ref+B*u_ref; xc_ref=F@xc_ref+G*y_ref
            xp=A@xp+B*u

        # Reconstruct x_c(t+K) from the last nu u/y values.
        U=np.asarray(us[-nu:]); Y=np.asarray(ys[-nu:])
        xc_block=Mu@U+My@Y

    return rows, O, Mu, My


def save(rows, out:Path):
    out.mkdir(parents=True,exist_ok=True)
    with (out/"sample_trace.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
    k=np.array([r['sample'] for r in rows]);ur=np.array([r['u_reference'] for r in rows]);u=np.array([r['u_test'] for r in rows]);eu=np.abs(u-ur)
    ep=np.array([r['plant_error_inf'] for r in rows]);xr=np.array([[r[f'xp_reference_{i}'] for i in range(4)] for r in rows]);x=np.array([[r[f'xp_test_{i}'] for i in range(4)] for r in rows])
    plt.figure();plt.plot(k,ur,label='reference');plt.plot(k,u,'--',label='u/y reconstruction');plt.xlabel('sample');plt.ylabel('u');plt.legend();plt.tight_layout();plt.savefig(out/'u.png',dpi=160);plt.close()
    plt.figure();
    for i in range(4):plt.plot(k,xr[:,i],label=f'ref x_p[{i}]');plt.plot(k,x[:,i],'--',label=f'recon x_p[{i}]')
    plt.xlabel('sample');plt.ylabel('plant state');plt.legend(ncol=2,fontsize=8);plt.tight_layout();plt.savefig(out/'plant_state.png',dpi=160);plt.close()
    plt.figure();plt.semilogy(k,np.maximum(eu,1e-18));plt.xlabel('sample');plt.ylabel('|u error|');plt.tight_layout();plt.savefig(out/'u_error.png',dpi=160);plt.close()
    plt.figure();plt.semilogy(k,np.maximum(ep,1e-18));plt.xlabel('sample');plt.ylabel('plant error inf-norm');plt.tight_layout();plt.savefig(out/'plant_state_error.png',dpi=160);plt.close()


def main():
    a=parse_args();p=get_profile(a.profile);out=Path(a.results_dir or (Path(__file__).resolve().parent/'results'/a.profile))
    rows,O,Mu,My=simulate(p,a.blocks,a.period,a.nu);save(rows,out);m=profile_metrics(p,a.period)
    print(f"profile={p.name}: {p.label}")
    print(f"K={a.period}, nu={a.nu}, rho(F)={m['rho_F']:.9f}, rho(Acl)={m['rho_Acl']:.9f}")
    print(f"rank(O4)={np.linalg.matrix_rank(O)}, cond(O4)={np.linalg.cond(O):.4f}")
    print(f"||Mu||_inf={np.linalg.norm(Mu,np.inf):.6g}, ||My||_inf={np.linalg.norm(My,np.inf):.6g}")
    print(f"max |u_error|={max(abs(r['u_error']) for r in rows):.3e}")
    print(f"max plant_error_inf={max(r['plant_error_inf'] for r in rows):.3e}")
    print(f"saved: {out}")
if __name__=='__main__':main()
