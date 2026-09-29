# Terminology

This repository avoids `common-mode` as the primary description because the intended interpretation is control-oriented rather than circuit-oriented.

For one replicated scalar control ciphertext, let

$$
\hat{\mathbf u}=u\mathbf1_4+e,
\qquad
\Pi=\frac14\mathbf1_4\mathbf1_4^T.
$$

We use:

## Equivalent scalar control-input error

$$
e_{scalar}=\Pi e.
$$

All slots contain the same error component. This is equivalent to the scalar perturbation

$$
u\mapsto u+\Delta u.
$$

It can therefore be interpreted in the usual physical control-input channel.

## Replication-inconsistency error

$$
e_{repl}=(I-\Pi)e.
$$

This is the component that violates the intended replicated representation `[u,u,u,u]`. It cannot be represented by one scalar plant-input disturbance and is therefore treated as an implementation-internal numerical state.

## Controller-internal numerical error

For the same measured output $y_{physical}$, propagate

$$
x_c^{ref,+}=Fx_c^{ref}+Gy_{physical}
$$

and define

$$
\delta x_c=Dec(c_x)-x_c^{ref}.
$$

This diagnostic isolates controller-side numerical propagation. It deliberately removes the correction that would arise from comparing two full physical closed loops, so it should not be called the physical closed-loop state error.

## Physical control error

Primary performance metric:

$$
e_u=u_{CKKS}-u_{semantic}.
$$

The repository also records `u_nominal`; `u_semantic` is the plaintext quantity corresponding to the actual physical trajectory being driven by the decoded encrypted control.

## Plant-state error

Primary closed-loop trajectory metric:

$$
e_p=x_{p,CKKS}-x_{p,nominal}.
$$

