# Controller profiles

All examples use the same 4-state SISO plant

```math
x_p(k+1)=Ax_p(k)+Bu(k),\qquad y(k)=Cx_p(k),
```

and the same controller form

```math
x_c(k+1)=Fx_c(k)+Gy(k),\qquad u(k)=Hx_c(k).
```

Select a preset with `--profile` in Python or with the first argument of the OpenFHE run scripts.

## Presets

| profile | $\rho(F)$ | $\rho(A_{cl})$ | $g_F(4)$ | $g_{HF}(4)$ | $\mathrm{cond}(\mathcal O_4)$ | purpose |
|---|---:|---:|---:|---:|---:|---|
| `stable` | 0.6000 | 0.6893 | 0.600 | 0.700 | 29.27 | numerical baseline |
| `unstable_low` | 1.0460 | 0.7852 | 1.197 | 3.841 | 25.60 | unstable $F$, mild internal transient |
| `unstable_high` | 1.0460 | 0.7852 | 32.658 | 24.931 | 767.22 | unstable $F$, large internal transient / cancellation stress test |

The transient indicators are

```math
g_F(K)=\max_{1\le i\le K}\|F^i\|_\infty,
\qquad
g_{HF}(K)=\max_{0\le i<K}\|HF^i\|_1.
```

They are **diagnostics of the controller realization**, not closed-loop stability criteria.

## Why `unstable_high` is useful

`unstable_high` is not a different ideal input/output controller. It is constructed from `unstable_low` by a non-normal similarity transformation

```math
F_h=T F_\ell T^{-1},\qquad
G_h=T G_\ell,\qquad
H_h=H_\ell T^{-1},
```

with the initial controller state transformed as

```math
x_{c,h}(0)=T x_{c,\ell}(0).
```

Therefore, in exact arithmetic, `unstable_low` and `unstable_high` have the same $u$, $y$, and plant trajectory. They also have the same eigenvalues and the same $\rho(A_{cl})$. What changes is the numerical representation: the high-transient realization has much larger intermediate coefficients, larger $\|F^i\|$, and poorer observability conditioning.

That makes the low/high pair a controlled experiment for the question:

> Can two exactly equivalent controller realizations behave differently after CKKS approximation, multiplication, rescaling, and bootstrapping?

This is why the repository uses **internal transient** rather than plant overshoot to define low/high transient.

## Recommended order

Run every implementation first with `stable`, then `unstable_low`, then `unstable_high`.

For the OpenFHE examples, start with 3 or 10 blocks before a long run. The high-transient profile is intentionally numerically difficult; an early CKKS failure should not be called physical closed-loop instability unless the actual $u$ and plant trajectories demonstrate it.
