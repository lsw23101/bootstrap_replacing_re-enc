# 03 — OpenFHE CKKS: PT-CT u/y reconstruction

이 단계에서는 controller state와 $u,y$ 신호를 ciphertext로 처리하고, controller/reconstruction coefficient는 plaintext로 둡니다.

## Ciphertext layout

Controller state:

```math
c_x\leftrightarrow[x_1,x_2,x_3,x_4].
```

각 measurement는 서로 다른 ciphertext입니다.

```math
c_{y_i}\leftrightarrow[y(t+i),0,0,0].
```

각 control도 서로 다른 ciphertext입니다.

```math
c_{u_i}\leftrightarrow[u(t+i),0,0,0].
```

즉 `[u(t),u(t+1),u(t+2),u(t+3)]`를 한 ciphertext에 packing하지 않으며, `[u,u,u,u]` 형태도 사용하지 않습니다.

## Control 계산

```math
u(t+i)=P_i x(t)+\sum_{j<i}q_{ij}y(t+j).
```

$P_i x$는 state rotation과 plaintext one-hot mask를 사용하여 slot 0에만 결과가 남도록 계산합니다.

Reconstruction도 각 sparse scalar $u_i,y_i$를 state slot으로 scatter한 뒤 $M_u,M_y$ coefficient를 곱합니다.

자세한 layout은 [`../docs/ENCRYPTED_REPRESENTATION.md`](../docs/ENCRYPTED_REPRESENTATION.md)를 참고하세요.

## 실행

```bash
export OPENFHE_ROOT="$HOME/openfhe-development"
./scripts/build.sh
./scripts/run.sh stable 10
./scripts/run.sh unstable_low 10
./scripts/run.sh unstable_high 10
```

현재 OpenFHE 구현은 `K=nu=4`로 고정되어 있습니다.

결과 폴더에는 다음 기본 plot이 생성됩니다.

```text
u.png
plant_state.png
u_error.png
plant_state_error.png
```

`sample_trace.csv`에서는 `u_active_error`, `u_inactive_leakage`, `y_active_error`, `y_inactive_leakage`도 확인할 수 있습니다.
