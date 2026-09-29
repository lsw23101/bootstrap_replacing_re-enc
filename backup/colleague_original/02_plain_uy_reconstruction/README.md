# 02 — Plain periodic u/y reconstruction

Purpose: verify the exact reconstruction architecture before encryption.

For $\nu=4$,

```math
U=\mathcal O_4 x + T_4Y,
\qquad
x(t+4)=M_uU+M_yY.
```

The script also allows $K>\nu$: it reconstructs the state at the end of each period from the last $\nu$ controls and measurements.

```bash
python3 uy_reconstruction_plain.py --profile stable --blocks 50 --period 4 --nu 4
python3 uy_reconstruction_plain.py --profile unstable_low --blocks 50 --period 4 --nu 4
python3 uy_reconstruction_plain.py --profile unstable_high --blocks 50 --period 4 --nu 4
```
