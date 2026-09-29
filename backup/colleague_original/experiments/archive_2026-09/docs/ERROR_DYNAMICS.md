# Block error dynamics: what is being analyzed

This note separates the physical closed-loop dynamics from the numerical-error dynamics induced by the recursive CKKS implementation.

## 1. Plant and controller

The nominal system is

\[
x_p(k+1)=Ax_p(k)+Bu(k),\qquad y(k)=Cx_p(k),
\]

\[
x_c(k+1)=Fx_c(k)+Gy(k),\qquad u(k)=Hx_c(k).
\]

The physical augmented closed-loop matrix is

\[
A_{cl}=\begin{bmatrix}A&BH\\GC&F\end{bmatrix}.
\]

For the controller used in this repository,

\[
\rho(F)=1.0460417098754866>1,
\qquad
\rho(A_{cl})=0.7851745975760147<1.
\]

Hence an unstable controller realization does not imply an unstable ideal closed loop.

## 2. Refresh period K and reconstruction horizon nu are different parameters

Let the encrypted controller be refreshed every \(K\) samples. Let \(\nu\) denote the number of input/output samples used for reconstruction.

For a block starting at \(t_j=jK\), define

\[
s_j=t_j+K-\nu.
\]

Using the last \(\nu\) controls and measurements in the block,

\[
U_j=\begin{bmatrix}u(s_j)&\cdots&u(s_j+\nu-1)\end{bmatrix}^{T},
\]

\[
Y_j=\begin{bmatrix}y(s_j)&\cdots&y(s_j+\nu-1)\end{bmatrix}^{T}.
\]

Define the observability matrix

\[
\mathcal O_\nu=
\begin{bmatrix}
H\\HF\\\vdots\\HF^{\nu-1}
\end{bmatrix}.
\]

When \(\mathcal O_\nu\) has full column rank, choose a left inverse \(L_\nu\) such that

\[
L_\nu\mathcal O_\nu=I.
\]

Then

\[
x_c(t_j+K)=M_uU_j+M_yY_j,
\]

with

\[
M_u=F^{\nu}L_\nu,
\]

and the corresponding \(M_y\) obtained from the lifted measurement terms.

In the current baseline,

\[
K=\nu=4,
\]

but this equality is an implementation choice, not a theoretical requirement.

## 3. Inherited controller-internal numerical error

To isolate controller-side numerical propagation, compare the encrypted controller state with a plaintext reference controller driven by the same measured output.

Let the block-boundary controller-state error be

\[
\delta x_j=\operatorname{Dec}(c_x(t_j))-x_c^{ref}(t_j).
\]

Assume temporarily that no new HE error is injected during the block. At the beginning of the reconstruction window,

\[
\delta x(s_j)=F^{K-\nu}\delta x_j.
\]

The induced control-history error is

\[
\delta U_j=\mathcal O_\nu F^{K-\nu}\delta x_j.
\]

Reconstruction gives

\[
\begin{aligned}
\delta x_{j+1}
&=M_u\delta U_j\\
&=F^{\nu}L_\nu\mathcal O_\nu F^{K-\nu}\delta x_j\\
&=F^K\delta x_j.
\end{aligned}
\]

Therefore the homogeneous block-level controller-internal error dynamics are

\[
\boxed{\delta x_{j+1}=F^K\delta x_j.}
\]

For the current baseline \(K=4\),

\[
\rho(F)^4=1.1972807715668585.
\]

This is the factor repeatedly observed in the same-measurement diagnostic runs.

The key interpretation is:

> The observability-based reconstruction is a state representation refresh, not an error-correction reset. If the control history was generated from a perturbed controller state, the reconstruction reproduces the corresponding perturbed state.

## 4. New HE error injection

A useful block-level model is

\[
\boxed{
\delta x_{j+1}
=F^K\delta x_j+d_j,
}
\]

where \(d_j\) collects new errors generated during the block. A more structured decomposition is

\[
d_j=M_ub_j+R_\nu r_j+a_j.
\]

Here:

- \(b_j\): actuator-consistent scalar control error introduced by bootstrap/arithmetic,
- \(r_j\): replicated-control representation inconsistency,
- \(R_\nu r_j\): state perturbation induced by that inconsistency during reconstruction,
- \(a_j\): remaining CKKS arithmetic error from multiplication, rescaling, alignment, and reconstruction.

The explicit solution is

\[
\delta x_j
=F^{Kj}\delta x_0+
\sum_{\ell=0}^{j-1}F^{K(j-1-\ell)}d_\ell.
\]

Thus even bounded per-block injection can be amplified in an unstable eigendirection of \(F^K\).

## 5. Control and plant errors remain the primary performance metrics

Controller-state numerical error is diagnostic. The control error inside a block can be written schematically as

\[
e_{u,j,i}=HF^i\delta x_j+\eta_{u,j,i},
\qquad i=0,\ldots,K-1,
\]

where \(\eta_{u,j,i}\) denotes newly introduced encrypted-arithmetic error in the control calculation.

The plant error then satisfies

\[
e_p(k+1)=Ae_p(k)+Be_u(k).
\]

Therefore the actual encrypted-control performance should be judged by quantities such as

\[
|e_u(k)|,
\qquad
\|e_p(k)\|,
\]

while \(\delta x_j\) is used to diagnose the internal mechanism that produces them.

## 6. What K and nu control

The two block parameters have different roles:

- \(K\): inherited error propagation through \(F^K\), and refresh frequency.
- \(\nu\): observability/reconstruction conditioning through \(\mathcal O_\nu\), \(L_\nu\), and \(M_u\).

The average number of bootstraps per sample is approximately

\[
\frac{\nu}{K}
\]

for the current reconstruction architecture.

This suggests a joint design problem involving

\[
\rho(F)^K,
\qquad
\|M_u\|,
\qquad
\Gamma_B,
\qquad
\Gamma_R,
\qquad
\frac{\nu}{K}.
\]

A candidate short-horizon bootstrap-error gain is

\[
\Gamma_B=
\max_{0\le i<K}\|HF^iM_u\|_1.
\]

A corresponding representation-inconsistency gain can be defined from the reconstruction columns \(m_\ell=M_u[:,\ell]\) as

\[
\Gamma_R=
\max_{0\le i<K,\ 0\le \ell<\nu}
\left\|HF^i\operatorname{diag}(m_\ell)\right\|.
\]

These gains describe injection sensitivity; they do not replace the homogeneous \(F^K\) propagation term.

## 7. Open HE-side question

The 149-block long-horizon run shows that the production bootstrap error is not well modeled as a small time-independent disturbance floor. The nominal plaintext/control magnitude remains inside the designed range, while the production ciphertext becomes progressively harder to bootstrap/decode.

The unresolved HE-side question is therefore to identify a useful state variable or quality metric \(Q(c)\) such that the bootstrap/arithmetic error can be related to ciphertext history:

\[
\|\mathcal E_B(c)\|\le \phi(Q(c)).
\]

This is separate from the plant closed-loop stability question and is one of the main targets of the next experiments.
