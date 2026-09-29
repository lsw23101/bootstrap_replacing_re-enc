#!/usr/bin/env python3
"""Verify the algebra used by the sparse-scalar OpenFHE layouts.

This script does not use OpenFHE.  It emulates only the slot operations:
- controller state: [x1,x2,x3,x4]
- each scalar y/u: [s,0,0,0]
- positive control rotation: left rotation
- reconstruction scatter: right rotation of a sparse scalar
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config.controller_profiles import get_profile, reconstruction_matrices


def rot_left(v: np.ndarray, k: int) -> np.ndarray:
    return np.roll(v, -k)


def rot_right(v: np.ndarray, k: int) -> np.ndarray:
    return np.roll(v, k)


def sparse_inner(row: np.ndarray, state: np.ndarray) -> np.ndarray:
    out = np.zeros(4)
    for j in range(4):
        mask = np.zeros(4)
        mask[0] = row[j]
        out += mask * rot_left(state, j)
    return out


def sparse_scale(scalar_ct: np.ndarray, coefficient: float) -> np.ndarray:
    mask = np.array([coefficient, 0.0, 0.0, 0.0])
    return mask * scalar_ct


def scatter(column: np.ndarray, scalar_ct: np.ndarray) -> np.ndarray:
    out = np.zeros(4)
    for r in range(4):
        mask = np.zeros(4)
        mask[r] = column[r]
        out += mask * rot_right(scalar_ct, r)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", default="unstable_low")
    args = ap.parse_args()

    p = get_profile(args.profile)
    _, T, Mu, My = reconstruction_matrices(p.F, p.G, p.H, 4)

    # Use a deterministic controller state and four deterministic measurements.
    x0 = p.xc0.copy()
    y = np.array([0.13, -0.07, 0.035, -0.02])

    # Lifted P and q for u(t+i) from the block-start controller state.
    P = np.vstack([p.H @ np.linalg.matrix_power(p.F, i) for i in range(4)])
    q = np.zeros((4, 4))
    for i in range(4):
        for j in range(i):
            q[i, j] = p.H @ np.linalg.matrix_power(p.F, i - 1 - j) @ p.G

    y_ct = [np.array([value, 0.0, 0.0, 0.0]) for value in y]
    u = np.zeros(4)
    u_ct = []
    for i in range(4):
        ct = sparse_inner(P[i], x0)
        for j in range(i):
            ct += sparse_scale(y_ct[j], q[i, j])
        expected = P[i] @ x0 + q[i] @ y
        u[i] = expected
        u_ct.append(ct)
        np.testing.assert_allclose(ct, [expected, 0.0, 0.0, 0.0], atol=1e-14, rtol=0.0)

    x_recon_slots = np.zeros(4)
    for j in range(4):
        x_recon_slots += scatter(Mu[:, j], u_ct[j])
        x_recon_slots += scatter(My[:, j], y_ct[j])

    x_recon_matrix = Mu @ u + My @ y
    np.testing.assert_allclose(x_recon_slots, x_recon_matrix, atol=1e-13, rtol=0.0)

    print(f"profile={p.name}")
    print("u sparse layout checks: PASS")
    print("reconstruction scatter check: PASS")
    print("u =", u)
    print("x_reconstruction =", x_recon_slots)
    print("max reconstruction mismatch =", np.max(np.abs(x_recon_slots - x_recon_matrix)))


if __name__ == "__main__":
    main()
