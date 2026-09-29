# 04 — OpenFHE CKKS: CT-CT u/y reconstruction

이 단계에서는 state/signals뿐 아니라 controller 및 reconstruction coefficient도 ciphertext로 암호화합니다.

## Signal layout

```math
c_x\leftrightarrow[x_1,x_2,x_3,x_4],
```

```math
c_{y_i}\leftrightarrow[y(t+i),0,0,0],
```

```math
c_{u_i}\leftrightarrow[u(t+i),0,0,0].
```

각 $u(t+i)$와 $y(t+i)$는 **서로 다른 ciphertext**입니다.

## Encrypted coefficients

$P_i x$ 계산에서는 각 $P_i[j]$를

```math
[P_i[j],0,0,0]
```

형태의 별도 encrypted mask로 저장합니다. State를 rotation하여 $x_j$를 slot 0으로 가져온 뒤 ct-ct multiplication을 수행합니다.

$q_{ij}$도

```math
[q_{ij},0,0,0]
```

형태로 암호화합니다.

Reconstruction coefficient $M_u,M_y$는 output state slot별 one-hot encrypted mask로 저장합니다. Sparse scalar $u_i,y_i$를 각 state slot으로 rotation한 뒤 해당 mask와 ct-ct multiplication합니다.

이 구현은 replicated scalar `[u,u,u,u]` 또는 `[y,y,y,y]`를 사용하지 않습니다.

자세한 식은 [`../docs/ENCRYPTED_REPRESENTATION.md`](../docs/ENCRYPTED_REPRESENTATION.md)를 참고하세요.

## 실행

```bash
export OPENFHE_ROOT="$HOME/openfhe-development"
./scripts/build.sh
./scripts/run.sh stable 10
./scripts/run.sh unstable_low 10
./scripts/run.sh unstable_high 10
```

처음에는 3~10 blocks로 smoke test를 권장합니다.

성능 판단은 controller-state error 하나가 아니라 실제 $u$와 plant trajectory를 기준으로 합니다.
