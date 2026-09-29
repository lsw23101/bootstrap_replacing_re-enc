# CKKS Dynamic Control with u/y Reconstruction

This repository provides a step-by-step implementation of a 4-state SISO dynamic controller, periodic state reconstruction from control/measurement histories, and OpenFHE CKKS realizations.

## Tutorial path

1. `01_plain_tac` — plain TAC reformulation
2. `02_plain_uy_reconstruction` — plain u/y reconstruction
3. `03_openfhe_ptct_uy_reconstruction` — ciphertext signals with plaintext coefficients
4. `04_openfhe_ctct_uy_reconstruction` — ciphertext signals and ciphertext coefficients

## Important ciphertext layout

The controller state is packed as

```math
c_x\leftrightarrow[x_1,x_2,x_3,x_4].
```

Every measurement and every control sample is encrypted **separately** as a sparse scalar:

```math
c_{y_i}\leftrightarrow[y(t+i),0,0,0],
```

```math
c_{u_i}\leftrightarrow[u(t+i),0,0,0].
```

The OpenFHE implementations do **not** use replicated scalar layouts such as `[u,u,u,u]` or `[y,y,y,y]`.

See [`docs/ENCRYPTED_REPRESENTATION.md`](docs/ENCRYPTED_REPRESENTATION.md) for the exact control and reconstruction circuits.

## Profiles

```text
stable
unstable_low
unstable_high
```

The low/high unstable pair has the same ideal control/plant trajectory but different internal numerical sensitivity.

## OpenFHE quick start

```bash
export OPENFHE_ROOT="$HOME/openfhe-development"
./03_openfhe_ptct_uy_reconstruction/scripts/build.sh
./03_openfhe_ptct_uy_reconstruction/scripts/run.sh stable 10

./04_openfhe_ctct_uy_reconstruction/scripts/build.sh
./04_openfhe_ctct_uy_reconstruction/scripts/run.sh stable 10
```

Previous long-horizon results from the older replicated-scalar implementation are not treated as results of the current architecture and have been removed from the shared v4 archive.
