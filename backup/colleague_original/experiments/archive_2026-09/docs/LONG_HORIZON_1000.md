# Long-horizon run requested for 1000 blocks

This experiment used the same baseline controller and CKKS configuration as the same-measurement diagnostic, but requested 1000 blocks.

## Outcome

The run did not reach 1000 blocks. It completed 149 full blocks (596 samples) and then failed during the third sample of block 149 after 598 completed samples.

The OpenFHE failure was:

```text
Decode(): The decryption failed because the approximation error is too high.
```

The failure tracker reports the most recent online operation as `q_2_1`; the actual exception occurs at the subsequent ciphertext decode used to recover the control value.

## Internal numerical mode

The theoretical block growth factor is

\[
\rho(F)^4=1.1972807715668585.
\]

The same-measurement controller-internal error fit is:

| fit interval | measured factor/block |
|---|---:|
| blocks 20--80 | 1.1972027 |
| blocks 40--100 | 1.1972852 |
| blocks 80--130 | 1.1974061 |
| blocks 120--148 | 1.1974262 |
| all blocks 20+ | 1.1973239 |

The unstable component exceeded the following magnitudes at approximately:

| threshold | first block |
|---:|---:|
| \(10^{-6}\) | 71 |
| \(10^{-4}\) | 97 |
| \(10^{-3}\) | 109 |
| \(10^{-2}\) | 122 |
| \(10^{-1}\) | 135 |
| \(1\) | 148 |

The final completed block has controller-internal unstable-component error magnitude about 1.11.

## Physical performance before decode failure

The largest observed plant-state error over the completed blocks is approximately

\[
1.75\times10^{-2},
\]

and the largest nominal-control error in the recorded sample trace is approximately

\[
1.66\times10^{-1}.
\]

Thus the cryptographic/numerical failure occurs after the controller-internal numerical mode has become large, but before the recorded plant trajectory has experienced a comparable order-one state error.

## Bootstrap behavior

The decoded plaintext magnitude remains within the designed order-one range, with the maximum control magnitude near 0.50. Nevertheless, production bootstrap error and slot inconsistency grow substantially with ciphertext history. The largest recorded production bootstrap-added error is approximately 0.38.

The RNS-chain position remains structurally unchanged during the late blocks: the reconstructed controller state remains at 3 towers and the pre-bootstrap control path remains at the same configured level/tower pattern. This is therefore not a simple case of consuming one additional tower every block until the chain reaches zero.

## Interpretation

This long-horizon run strengthens two empirical conclusions:

1. the controller-internal numerical mode continues to track the unstable \(F^4\) mode over a much longer interval than the earlier 100-block test;
2. production ciphertext quality degrades with recursive history even though the nominal plaintext range remains small.

The run does **not** establish that \(\rho(F)>1\) alone causes physical closed-loop instability. The plant/controller closed loop remains a separate stable dynamical system; the unresolved issue is how the recursive CKKS numerical state couples into the physical control channel.
