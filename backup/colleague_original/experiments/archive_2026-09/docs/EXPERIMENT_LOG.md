# Experiment log and current conclusions

## E0 — Plaintext range analysis

Purpose: determine whether signal magnitude or intermediate overflow is a plausible cause of CKKS failure.

Result: all relevant plaintext signal/intermediate magnitudes remain below 0.8; a fixed working range of 1 is sufficient.

Interpretation: later error growth is not explained by a large nominal dynamic range in this controller.

## E1 — Fixed-range ct-ct long run

Directory: `results/range_100/`

- 100 blocks = 400 samples
- 4 production bootstraps per block
- max semantic control error: 1.60e-5
- max plant-state error: 2.02e-6
- fitted control-error growth after block 20: about 1.1955/block

The physical trajectory remains close for a long time, but a clear long-horizon growth trend appears.

## E2 — Same-measurement controller-internal error

Directory: `results/mode_100/`

The reference controller receives exactly the same physical measurement as the encrypted controller.

Theoretical unstable block factor:

$$
\rho(F)^4=1.19728077.
$$

Measured log-fit factor:

$$
1.19728116/\text{block}.
$$

Interpretation: a controller-internal numerical error mode follows $F^4$ extremely closely. This does not by itself imply that the full physical closed-loop error obeys $F^4$.

## E3 — Post-bootstrap RNS-tower budget sweep

Directory: `results/budget_sweep_100/`

Cases: L3, L5, L7.

The absolute error amplitude changes, but the controller-internal growth factor stays essentially equal to $\rho(F)^4$:

- L3: 1.19723
- L5: 1.19729
- L7: 1.19730

Interpretation: in this tested range, insufficient remaining tower count is not the primary cause of the exponential slope.

## E4 — Replicated scalar-control error decomposition

Directory: `results/slot_50/`

For each bootstrapped control ciphertext, decompose the slotwise error into:

1. equivalent scalar control-input error,
2. replication-inconsistency error.

Over the last ten blocks:

- mean scalar-equivalent error: 2.04e-10
- mean replication-inconsistency error: 4.22e-10

Fitted growth factors from block 10:

- scalar-equivalent: 1.19668/block
- replication inconsistency: 1.19255/block

Interpretation: a non-scalar encrypted-implementation error exists and grows at a rate comparable to the unstable controller-internal mode. Causality has not yet been established.

## What is supported by the current data

- `rho(F) > 1` does not imply immediate encrypted-control failure.
- Short-horizon sensitivity and long-horizon numerical stability are different issues.
- The controller-internal numerical error has a highly reproducible block growth near `rho(F)^4`.
- Increasing the tested post-bootstrap tower budget does not remove that slope.
- A replication-inconsistency component is present in the scalar-control ciphertexts.

## What is not yet established

- The replication-inconsistency component has not yet been shown to be the cause of physical control-error growth.
- The current experiments do not establish an infinite-horizon theorem.
- The current CKKS parameters are for numerical diagnosis; they are not presented as a production security parameter set.



## E5 — Requested 1000-block long-horizon run

Directory: `results/mode_1000/`

The run requested 1000 blocks but completed 149 full blocks (596 samples) and failed during block 149 after 598 completed samples. OpenFHE reported `Decode(): The decryption failed because the approximation error is too high.`

The same-measurement controller-internal unstable component continued to track the theoretical block factor

$$
\rho(F)^4=1.19728077.
$$

Measured log-fit factors were approximately 1.19720 over blocks 20--80, 1.19729 over blocks 40--100, 1.19741 over blocks 80--130, and 1.19743 over blocks 120--148.

The unstable component crossed magnitude 1 at block 148. The largest recorded plant-state error before the failure was about 1.75e-2 and the largest recorded nominal-control error about 1.66e-1.

The control/plaintext magnitude remained inside the order-one design range, while the production bootstrap-added error grew as large as roughly 0.38. This points away from nominal message-range overflow and toward recursive ciphertext-history/numerical-state degradation.

See `docs/LONG_HORIZON_1000.md`.
