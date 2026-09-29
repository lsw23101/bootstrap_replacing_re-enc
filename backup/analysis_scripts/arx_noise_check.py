"""Plain (double) closed loop with the delayed-window ARX formula plus small injected noise.
Noise models CKKS error: n_u on every applied u, n_bt on the bootstrapped window copy.
If the error grows geometrically, the window scheme itself has an unstable error mode."""
import sys
import numpy as np

sys.path.insert(0, "/Users/nico/All/uy-reconstruction")
from config.controller_profiles import A, B, C, XP0, get_profile, observability_matrix

mp = np.linalg.matrix_power
rng = np.random.default_rng(0)


def run(profile, K, d, steps, sig_u=1e-12, sig_bt=1e-12):
    p = get_profile(profile)
    F, G, H = p.F, p.G, p.H
    n = 4
    R = mp(F, n) @ np.linalg.inv(observability_matrix(F, H, n)) @ np.eye(n)[:, -1]
    Fb = F - np.outer(R, H)
    HFm = [H @ mp(F, m) for m in range(d + K + 2 * n + 70)]
    kstar = -(-(n - 1 + d) // K)
    xp, xc = XP0.copy(), p.xc0.copy()
    xpr, xcr = XP0.copy(), p.xc0.copy()
    u_app, u_bt, y_h = {}, {}, {}
    blk_err = {}
    for t in range(steps):
        y = C @ xp; y_h[t] = y
        yr = C @ xpr; ur = H @ xcr
        k = t // K
        if k < kstar:
            u = HFm[t] @ p.xc0 + sum(HFm[l] @ G * y_h[t - 1 - l] for l in range(t))
        else:
            tau = k * K - d; m = t - tau - 1
            u = sum(HFm[m] @ mp(Fb, j) @ R * u_bt[tau - j] + HFm[m] @ mp(Fb, j) @ G * y_h[tau - j] for j in range(n))
            u += sum(HFm[l] @ G * y_h[t - 1 - l] for l in range(m))
        u += sig_u * rng.standard_normal()
        u_app[t] = u
        u_bt[t] = u + sig_bt * rng.standard_normal()   # bootstrapped copy (imag already cleared)
        blk_err[k] = max(blk_err.get(k, 0.0), abs(u - ur))
        xp = A @ xp + B * u
        xpr = A @ xpr + B * ur; xcr = F @ xcr + G * yr
    e = np.array([blk_err[k] for k in sorted(blk_err)])
    return e


for prof in ["unstable_low", "unstable_high", "stable"]:
    for K in [4, 8, 16, 32, 64]:
        steps = max(12 * K, 200)
        e = run(prof, K, K, steps)
        tail = e[-4:]
        ratio = (tail[-1] / tail[0]) ** (1 / 3) if tail[0] > 0 else np.nan
        print(f"{prof:14s} K=d={K:3d} steps={steps:4d} block-max err: first={e[2]:.1e} last={e[-1]:.1e} "
              f"growth/block (last 4)={ratio:.2f}")
