# Background and preliminary observations

This repository focuses on the current 4-state OpenFHE ct-ct experiment, but the experimental direction was motivated by earlier tests.

## Preliminary 2-state DESILO study

A 2-state SISO plant/controller was used to compare unstable controller realizations with very different finite-horizon numerical gains.

### High-transient realization

- $\rho(F)\approx 3.8023$
- $\rho(A_{cl})\approx 0.9000$
- bootstrap-to-control amplification diagnostic $\Gamma_B\approx1.8\times10^{11}$

A small bootstrap perturbation was rapidly amplified through reconstruction and lifted control evaluation, producing large control error and subsequent trajectory deviation.

### Low-transient realization

- $\rho(F)\approx1.001249>1$
- $\rho(A_{cl})=0.9$
- $\Gamma_B\approx1.88$

This realization ran for roughly 700 samples before the WSL process was OOM-killed. The logged control/plant errors remained small up to that environmental termination.

### Preliminary lesson

These tests suggested that

$$
\rho(F)>1
$$

alone is not a useful short-horizon implementability criterion. Finite-horizon quantities such as

$$
HF^i,\qquad M_u,\qquad HF^iM_u
$$

can dominate the immediate effect of a newly injected CKKS/bootstrap perturbation.

The current 4-state study was designed to keep that short-horizon amplification small while investigating the remaining long-horizon numerical dynamics.

## Earlier 4-state OpenFHE observations

Before the current range-designed experiment:

- stable-$F$ realizations could complete long ct-ct runs;
- an unstable-$F$ realization with a stable physical closed loop showed long-horizon production-ciphertext degradation;
- fresh/semantic re-encryption diagnostics indicated that simply resetting ciphertext history could strongly reduce the observed error, motivating a closer study of production ciphertext lineage and bootstrap behavior.

Those observations are treated as motivation rather than final evidence in this repository. The reproducible raw data bundled here starts from the fixed-range ct-ct experiments documented in `results/`.

