# OpenFHE 신호 표현 및 연산 규칙

03/04 OpenFHE 구현에서 사용하는 ciphertext layout을 정리한 문서입니다.

## 1. Controller state

Controller state는 하나의 ciphertext에 4개 state를 packing합니다.

```math
c_x \leftrightarrow [x_1,x_2,x_3,x_4].
```

## 2. Measurement와 control

각 시점의 measurement와 control은 **각각 별도의 ciphertext**입니다. 시간축으로 4개 값을 한 ciphertext에 packing하지 않습니다.

```math
c_{y_i}\leftrightarrow [y(t+i),0,0,0],
```

```math
c_{u_i}\leftrightarrow [u(t+i),0,0,0].
```

따라서 한 block에서 `K=4`이면

```text
c_y0, c_y1, c_y2, c_y3   # 서로 다른 4개 ciphertext
c_u0, c_u1, c_u2, c_u3   # 서로 다른 4개 ciphertext
```

를 사용합니다.

## 3. Control 계산

Block 시작 state를 $x(t)$라고 할 때

```math
u(t+i)=P_i x(t)+\sum_{j=0}^{i-1}q_{ij}y(t+j).
```

$P_i=[p_{i1},p_{i2},p_{i3},p_{i4}]$라 두면 state를 rotation한 뒤 slot 0의 one-hot coefficient와 곱합니다.

```math
[p_{i1},0,0,0]\odot[x_1,x_2,x_3,x_4]
```

```math
+[p_{i2},0,0,0]\odot[x_2,x_3,x_4,x_1]
```

```math
+[p_{i3},0,0,0]\odot[x_3,x_4,x_1,x_2]
```

```math
+[p_{i4},0,0,0]\odot[x_4,x_1,x_2,x_3]
```

결과는

```math
[P_i x(t),0,0,0]
```

입니다.

Measurement term도

```math
[q_{ij},0,0,0]\odot[y(t+j),0,0,0]
=
[q_{ij}y(t+j),0,0,0]
```

이므로 최종 control ciphertext는

```math
c_{u_i}\leftrightarrow[u(t+i),0,0,0]
```

형태를 유지합니다.

### PT-CT

`P_i`, `q_ij`의 one-hot mask는 plaintext coefficient입니다.

### CT-CT

각 one-hot mask도 setup 시 ciphertext로 암호화하여 재사용합니다.

## 4. Bootstrap

각 control ciphertext를 개별적으로 bootstrap합니다.

```math
c_{u_i}^{B}=Boot(c_{u_i}),\qquad i=0,1,2,3.
```

이상적인 semantic layout은 여전히

```math
[u(t+i),0,0,0]
```

입니다.

따라서 diagnostic은 기존의 replicated-slot spread가 아니라 다음 두 값을 봅니다.

```math
e_{u,\mathrm{active}}=|\hat u_0-u|,
```

```math
e_{u,\mathrm{inactive}}=\max_{r=1,2,3}|\hat u_r|.
```

Measurement ciphertext도 같은 방식으로 active-slot error와 inactive-slot leakage를 확인할 수 있습니다.

## 5. Reconstruction

Reconstruction은

```math
x_c(t+4)=M_uU+M_yY
```

이며

```math
U=[u(t),u(t+1),u(t+2),u(t+3)]^T,
```

```math
Y=[y(t),y(t+1),y(t+2),y(t+3)]^T.
```

각 $u(t+j)$와 $y(t+j)$는 sparse scalar ciphertext이므로, 단순히 full column과 Hadamard product하지 않습니다.

예를 들어 $M_u$의 $j$번째 column을

```math
m_j=[m_{1j},m_{2j},m_{3j},m_{4j}]^T
```

라 하면 $c_{u_j}=[u_j,0,0,0]$를 각 state slot으로 rotation한 뒤 one-hot coefficient와 곱합니다.

```math
[m_{1j}u_j,0,0,0]
+
[0,m_{2j}u_j,0,0]
+
[0,0,m_{3j}u_j,0]
+
[0,0,0,m_{4j}u_j].
```

따라서

```math
m_j u_j
```

가 packed controller-state contribution으로 만들어집니다.

$M_yY$도 동일한 방식으로 계산합니다.

이 scatter 방식은 inactive slot에 생긴 CKKS leakage를 다른 state component로 의도적으로 합산하지 않습니다.

## 6. 핵심 주의사항

- `[u,u,u,u]`는 03/04의 control representation이 아닙니다.
- `[y,y,y,y]`도 03/04의 measurement representation이 아닙니다.
- 각 $u(t+i)$와 $y(t+i)$는 별도의 sparse scalar ciphertext입니다.
- Controller state만 `[x1,x2,x3,x4]` 형태로 packing됩니다.
- 이전 replicated-scalar 구현에서 생성된 long-horizon 결과는 현재 architecture의 결과로 사용하면 안 됩니다.
