# Next experiment: replicated-scalar projection A/B test

## Objective

Test whether the loss of the replicated scalar-control structure is causally responsible for the long-horizon numerical growth.

## Baseline

For each block position $i$, bootstrap

$$
c_{u_i}^B=Boot(c_{u_i})
$$

and use it directly in the reconstruction.

## Projection case

Project the bootstrapped ciphertext onto the replicated-scalar subspace:

$$
\Pi=\frac14\mathbf1_4\mathbf1_4^T,
\qquad
c_{u_i}^{proj}=\Pi c_{u_i}^B.
$$

Then use $c_{u_i}^{proj}$ in the same reconstruction.

A rotation/add implementation can sum all four slots without an additional ct-ct multiplication. The factor `1/4` may be absorbed into the encrypted reconstruction coefficient $M_u[:,i]$, avoiding an extra online rescale.

## A/B controls

Keep fixed:

- plant/controller matrices
- OpenFHE version
- ring dimension and logical slots
- scaling/first modulus bits
- 2-iteration bootstrap
- post-bootstrap tower budget
- block length
- number of samples

Only the projection step should differ.

## Primary metrics

1. physical control error $|u_{CKKS}-u_{semantic}|$
2. plant-state trajectory error
3. controller-internal same-measurement error growth
4. replication-inconsistency error immediately after bootstrap

## Decision rule

If projection removes or substantially reduces the $\approx1.197^j$ growth while the baseline retains it, the replication-inconsistency component is causally important.

If the growth remains essentially unchanged after projection, the next target should be the scalar-consistent controller-internal numerical dynamics and the production-ciphertext bootstrap history.



## Secondary parameter study: separate K and nu

After the projection A/B test, hold the reconstruction horizon at `nu = 4` and vary the refresh period, for example

$$
K\in\{4,8,12,20\}.
$$

The predicted inherited homogeneous controller-internal factor is

$$
\rho(F)^K,
$$

while the average bootstrap frequency is approximately `nu/K` per sample. This experiment separates the error-propagation role of `K` from the reconstruction-conditioning role of `nu`.
