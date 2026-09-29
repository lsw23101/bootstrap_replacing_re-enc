# Model and encrypted architecture

## Plant/controller

$$
x_p(k+1)=Ax_p(k)+Bu(k),\qquad y(k)=Cx_p(k),
$$

$$
x_c(k+1)=Fx_c(k)+Gy(k),\qquad u(k)=Hx_c(k).
$$

The exact matrices used in the experiment are stored in `generated_siso_params_unstable_lowtransient.h`.

Key numerical properties:

$$
\rho(F)=1.0460417098754866,
$$

$$
\rho(A_{cl})=0.7851745975760147,
$$

$$
\operatorname{cond}(O_4)=25.5999232,
\qquad
\|M_u\|_\infty=1.2841548,
$$

$$
\Gamma_B=\max_{0\le i<4}\|HF^iM_u\|_1=0.8947301.
$$

## Four-sample lifted controller

For a block starting at $t$,

$$
u(t+i)=P_i x_c(t)+\sum_{j=0}^{i-1}q_{ij}y(t+j),
\qquad i=0,1,2,3.
$$

The four controls are separate replicated-scalar ciphertexts. They are each bootstrapped once and then used in

$$
x_c(t+4)=M_uU_t+M_yY_t.
$$

All entries of $P_i,q_{ij},M_u,M_y$ are encrypted once at setup. The online path therefore uses ct-ct coefficient multiplication.

## Per-block online operations

The main implementation reports approximately:

- 18 ct-ct multiplications
- 18 relinearizations
- 18 product rescales
- 8 rotations
- 21 additions
- 4 production bootstraps

## Why the same-measurement diagnostic is used

A physical closed-loop comparison mixes three effects:

1. controller numerical error,
2. resulting plant trajectory deviation,
3. plant-feedback correction back into the controller.

To identify the controller-internal propagation, the diagnostic reference receives the same physical measurement as the encrypted controller. This gives

$$
\delta x_{c,j+1}\approx F^4\delta x_{c,j}+\eta_j
$$

at block boundaries and lets us test whether an internal numerical mode follows the unstable eigenstructure of $F$.



## Refresh period K versus reconstruction horizon nu

The baseline uses `K = nu = 4`, but these are different design parameters.

- `K`: refresh/reconstruction period in samples.
- `nu`: number of input/output samples used to reconstruct the controller state.

For `K >= nu`, the final `nu` samples of a block can be used for reconstruction. If `L_nu O_nu = I`, then the homogeneous inherited controller-internal numerical error satisfies

$$
\delta x_{j+1}=F^K\delta x_j.
$$

Thus changing `nu` alone does not change the inherited homogeneous factor `F^K`; it changes reconstruction conditioning and the injection gains. Increasing `K` reduces average bootstrap frequency (`nu/K` bootstraps per sample in this architecture) but exposes the recursive state to a longer interval of `F` propagation.

The full derivation is in `docs/ERROR_DYNAMICS.md`.
