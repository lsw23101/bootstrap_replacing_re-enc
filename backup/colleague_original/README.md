# CKKS 기반 동적 제어기와 u/y Reconstruction 실험

이 저장소는 **동적 제어기의 내부 행렬 $F$가 stable 또는 unstable한 경우**, $u,y$ history를 이용한 controller-state reconstruction과 CKKS 기반 암호화 구현을 단계적으로 확인하기 위한 코드입니다.

처음 보는 사람도 구조를 따라갈 수 있도록 다음 네 단계로 구성했습니다.

1. `01_plain_tac` — Plain TAC formulation
2. `02_plain_uy_reconstruction` — Plain $u,y$-based reconstruction
3. `03_openfhe_ptct_uy_reconstruction` — OpenFHE CKKS, PT-CT coefficients
4. `04_openfhe_ctct_uy_reconstruction` — OpenFHE CKKS, CT-CT coefficients

> **중요:** 03/04에서 각 $u(t+i)$와 $y(t+i)$는 서로 다른 ciphertext입니다. Scalar signal은 `[value,0,0,0]` 형태로 표현하며 `[u,u,u,u]` 또는 `[y,y,y,y]` replicated packing을 사용하지 않습니다.

자세한 ciphertext layout은 [`docs/ENCRYPTED_REPRESENTATION.md`](docs/ENCRYPTED_REPRESENTATION.md)를 먼저 확인하세요.

---

## 1. 문제 설정

Plant:

```math
x_p(k+1)=Ax_p(k)+Bu(k)
```

```math
y(k)=Cx_p(k)
```

Dynamic controller:

```math
x_c(k+1)=Fx_c(k)+Gy(k)
```

```math
u(k)=Hx_c(k)
```

전체 closed-loop matrix는

```math
A_{cl}=
\begin{bmatrix}
A & BH\\
GC & F
\end{bmatrix}.
```

이 저장소에서는 특히

```math
\rho(F)>1,\qquad \rho(A_{cl})<1
```

인 controller realization에서도 CKKS implementation이 어떻게 동작하는지 비교합니다.

---

## 2. Controller profile

모든 단계에서 다음 세 가지 profile을 공통으로 사용합니다.

| profile | $\rho(F)$ | $\rho(A_{cl})$ | 특징 |
|---|---:|---:|---|
| `stable` | 0.6000 | 0.6893 | Schur-stable $F$ baseline |
| `unstable_low` | 1.0460 | 0.7852 | unstable $F$, 낮은 internal transient |
| `unstable_high` | 1.0460 | 0.7852 | unstable $F$, 큰 internal transient |

`unstable_high`는 `unstable_low`의 similarity-transformed realization입니다.

```math
F_h=TF_lT^{-1},\qquad G_h=TG_l,\qquad H_h=H_lT^{-1}.
```

초기 controller state도 함께 변환하므로 exact arithmetic에서는 두 profile의 ideal $u$와 plant trajectory가 동일합니다. 내부 realization의 numerical sensitivity만 크게 달라지도록 만든 비교 case입니다.

자세한 내용은 [`docs/CONTROLLER_PROFILES.md`](docs/CONTROLLER_PROFILES.md)를 참고하세요.

---

## 3. Repository 구조

```text
.
├── 01_plain_tac/
├── 02_plain_uy_reconstruction/
├── 03_openfhe_ptct_uy_reconstruction/
├── 04_openfhe_ctct_uy_reconstruction/
├── config/
├── include/
├── tools/
├── docs/
└── experiments/
```

모든 단계의 기본 결과는 가능한 한 다음 형식으로 통일합니다.

```text
u.png
plant_state.png
u_error.png
plant_state_error.png
sample_trace.csv
```

성공 여부는 controller-state error 하나가 아니라 **실제 control $u$와 plant trajectory $x_p$**를 우선해서 판단합니다. Controller-state error는 numerical mechanism을 보기 위한 diagnostic입니다.

---

## 4. Part 1 — Plain TAC

```bash
python3 01_plain_tac/tac_plain.py --profile stable
python3 01_plain_tac/tac_plain.py --profile unstable_low
python3 01_plain_tac/tac_plain.py --profile unstable_high
```

기본 controller

```math
x_c^+=Fx_c+Gy,\qquad u=Hx_c
```

와 TAC identity

```math
x_c^+=(F-RH)x_c+Gy+Ru
```

를 plain arithmetic에서 비교합니다.

---

## 5. Part 2 — Plain u/y reconstruction

기본값은 refresh period $K=4$, reconstruction horizon $\nu=4$입니다.

```bash
python3 02_plain_uy_reconstruction/uy_reconstruction_plain.py \
  --profile unstable_low --blocks 50 --period 4 --nu 4
```

Plain 단계에서는 $K$와 $\nu$를 분리할 수 있습니다.

```bash
python3 02_plain_uy_reconstruction/uy_reconstruction_plain.py \
  --profile unstable_low --blocks 20 --period 12 --nu 4
```

Reconstruction은

```math
x_c(t+K)=M_uU+M_yY
```

형태로 구현합니다.

---

## 6. Part 3 — OpenFHE PT-CT

여기부터 OpenFHE가 필요합니다.

```bash
export OPENFHE_ROOT="$HOME/openfhe-development"
./03_openfhe_ptct_uy_reconstruction/scripts/build.sh
./03_openfhe_ptct_uy_reconstruction/scripts/run.sh stable 10
./03_openfhe_ptct_uy_reconstruction/scripts/run.sh unstable_low 10
./03_openfhe_ptct_uy_reconstruction/scripts/run.sh unstable_high 10
```

### Signal representation

Controller state만 vector packing합니다.

```math
c_x\leftrightarrow[x_1,x_2,x_3,x_4].
```

각 measurement는 개별 sparse scalar ciphertext입니다.

```math
c_{y_i}\leftrightarrow[y(t+i),0,0,0].
```

각 control도 개별 sparse scalar ciphertext입니다.

```math
c_{u_i}\leftrightarrow[u(t+i),0,0,0].
```

PT-CT에서는 $P_i$, $q_{ij}$, $M_u$, $M_y$ coefficient를 plaintext one-hot mask로 사용합니다.

---

## 7. Part 4 — OpenFHE CT-CT

```bash
./04_openfhe_ctct_uy_reconstruction/scripts/build.sh
./04_openfhe_ctct_uy_reconstruction/scripts/run.sh stable 10
./04_openfhe_ctct_uy_reconstruction/scripts/run.sh unstable_low 10
./04_openfhe_ctct_uy_reconstruction/scripts/run.sh unstable_high 10
```

Signal layout은 Part 3과 동일합니다.

다만 controller/reconstruction coefficient도 setup에서 ciphertext로 암호화합니다.

Control 계산의 $P_i x$는 다음 구조입니다.

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

따라서 결과는

```math
[P_i x,0,0,0]
```

입니다.

$q_{ij}y_j$도 sparse scalar layout을 그대로 유지합니다.

Reconstruction에서는 각 sparse scalar $u_i,y_i$를 state slot으로 rotation하여 one-hot $M_u,M_y$ coefficient와 곱합니다. 단순 replicated broadcast-sum은 사용하지 않습니다.

---

## 8. OpenFHE diagnostic

03/04에서는 ideal scalar layout이

```math
[s,0,0,0]
```

이므로 다음을 기록합니다.

Active-slot error:

```math
e_{active}=|\hat s_0-s|
```

Inactive-slot leakage:

```math
e_{inactive}=\max_{r=1,2,3}|\hat s_r|.
```

`sample_trace.csv`와 `bootstrap_diagnostics.csv`에서 $u$와 $y$에 대해 이 값을 확인할 수 있습니다.

---

## 9. Python 환경

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Profile 검증:

```bash
python3 tools/validate_profiles.py
```

Sparse layout algebra 확인:

```bash
python3 tools/verify_sparse_layout.py
```

---

## 10. 현재 advanced experiment 상태

이전 replicated-scalar OpenFHE 구현에서 생성된 long-horizon 결과는 현재 architecture의 결과로 사용하지 않습니다.

현재 v4에서는 corrected sparse-scalar implementation을 먼저 03/04에서 재검증한 뒤 다음 실험을 다시 수행할 예정입니다.

- long-horizon PT-CT / CT-CT
- same-measurement controller-error dynamics
- post-bootstrap level/tower sweep
- active-slot error / inactive-slot leakage
- production-vs-fresh bootstrap comparison

수정 내역은 [`docs/CORRECTION_FROM_V3.md`](docs/CORRECTION_FROM_V3.md)를 참고하세요.

---

## 11. 권장 실행 순서

```text
01_plain_tac
    ↓
02_plain_uy_reconstruction
    ↓
03_openfhe_ptct_uy_reconstruction
    ↓
04_openfhe_ctct_uy_reconstruction
```

OpenFHE는 먼저 `stable` profile로 3~10 blocks를 확인한 뒤 `unstable_low`, `unstable_high` 순서로 진행하는 것을 권장합니다.
