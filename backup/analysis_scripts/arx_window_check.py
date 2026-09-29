"""Check: u(t) computed only from a delayed, block-reused u-window plus recent y's
reproduces the original controller output exactly (plain arithmetic, closed loop).

Block k covers t = kK .. kK+K-1.  Its u-window is W_k = [u(tau-n+1), ..., u(tau)],
tau = kK - d.  For t in block k, with m = t - tau - 1:
  u(t) = sum_j H F^m Fb^j R u(tau-j) + sum_j H F^m Fb^j G y(tau-j)
         + sum_{s=tau+1}^{t-1} H F^{t-1-s} G y(s),       Fb = F - R H, Fb^n = 0.
"""
import sys
import numpy as np

sys.path.insert(0, "/Users/nico/All/uy-reconstruction")
from config.controller_profiles import A, B, C, XP0, get_profile, observability_matrix

mp = np.linalg.matrix_power


def deadbeat_R(F, H):
    n = F.shape[0]
    return mp(F, n) @ np.linalg.inv(observability_matrix(F, H, n)) @ np.eye(n)[:, -1]


def run(profile, K, d, steps=400):
    p = get_profile(profile)
    F, G, H = p.F, p.G, p.H
    n = F.shape[0]
    R = deadbeat_R(F, H)
    Fb = F - np.outer(R, H)
    # reference closed loop with the original controller
    xp, xc = XP0.copy(), p.xc0.copy()
    u_ref, y_ref = [], []
    for _ in range(steps):
        y = C @ xp; u = H @ xc
        y_ref.append(y); u_ref.append(u)
        xp = A @ xp + B * u; xc = F @ xc + G * y
    u_ref, y_ref = np.array(u_ref), np.array(y_ref)
    # window-based ARX: the plant is driven by the formula output, so any mismatch
    # would propagate through the closed loop.
    xp = XP0.copy()
    u_hist, y_hist = {}, {}
    start = d + n  # warm-up: use the reference controller until a full window exists
    coef_max = 0.0
    xc = p.xc0.copy()
    for t in range(steps):
        y = C @ xp; y_hist[t] = y
        if t < start:
            u = H @ xc; xc = F @ xc + G * y
        else:
            k0 = (t // K) * K
            tau = k0 - d
            if tau - n + 1 < 0:
                u = H @ xc; xc = F @ xc + G * y
                u_hist[t] = u; xp = A @ xp + B * u
                continue
            m = t - tau - 1
            u = 0.0
            for j in range(n):
                cu = H @ mp(F, m) @ mp(Fb, j) @ R
                cy = H @ mp(F, m) @ mp(Fb, j) @ G
                coef_max = max(coef_max, abs(cu))
                u += cu * u_hist[tau - j] + cy * y_hist[tau - j]
            for s in range(tau + 1, t):
                u += H @ mp(F, t - 1 - s) @ G * y_hist[s]
            xc = F @ xc + G * y  # kept only so warm-up branch stays consistent
        u_hist[t] = u
        xp = A @ xp + B * u
    u_arx = np.array([u_hist[t] for t in range(steps)])
    return np.abs(u_arx - u_ref).max(), coef_max, np.abs(mp(Fb, n)).max()


print(f"{'profile':14s} {'K':>3s} {'d':>3s} {'max|u_arx-u_ref|':>17s} {'max|coef_u|':>12s} {'||Fb^n||':>9s}")
for prof in ["stable", "unstable_low", "unstable_high"]:
    for K, d in [(4, 1), (4, 4), (8, 8), (16, 16)]:
        e, c, fb = run(prof, K, d)
        print(f"{prof:14s} {K:3d} {d:3d} {e:17.2e} {c:12.3e} {fb:9.1e}")
