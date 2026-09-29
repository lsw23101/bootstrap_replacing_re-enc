# 05 — OpenFHE CKKS CT-CT: 지연 입력창 ARX 제어

상태 암호문 $x$를 유지하지 않고, **bootstrapping된 과거 입력창**과 **새로 암호화된 출력**만으로 매 step $u(t)$를 계산합니다.
같은 입력창을 $K$ step 동안 재사용하고, 그 사이에 다음 입력창을 bootstrapping합니다.
유도와 결과는 [`../tex/formulation.tex`](../tex/formulation.tex)에 있습니다.

## 기준 설정

128-bit ($N=65536$), $K=64$, $d=61$, 출력 패킹과 입력창 패킹, 입력 레벨 축소를 씁니다.

```bash
./arx_window_ctct --profile unstable_low --K 64 --steps 640 --results-dir OUT \
  --pack-y --pack-u --compact-inputs --secure128
```

`results/`에는 이 설정의 최종 검증 run이 있습니다 (`unstable_low`, `unstable_high`, 각 640 step).

| run | 최대 $\lvert u_{\mathrm{ckks}}-u_{\mathrm{plain}}\rvert$ | bootstrapping 1회 | $u$ 계산 (평균) | 최소 $T_s$ (모델) | 메모리 |
|---|---:|---:|---:|---:|---:|
| `unstable_high_K64_d61_s128` | $1.1\times10^{-7}$ | 6.1 s | 191 ms | 307 ms | 3.7 GB |
| `unstable_low_K64_d61_s128` | $3.3\times10^{-8}$ | 6.7 s | 203 ms | 303 ms | 3.8 GB |

## 계산식

블록 $k$는 $t=kK,\dots,kK+K-1$이고, 입력창 지연은 $\tau_k=kK-d$, $m=t-\tau_k-1$입니다.

```math
u(t)=\sum_{j=0}^{3}\alpha_{m,j}\,u_{\mathrm{bt}}(\tau_k-j)+\sum_{j=0}^{3}\beta_{m,j}\,y(\tau_k-j)+\sum_{l=0}^{m-1}\gamma_l\,y(t-1-l)
```

```math
\alpha_{m,j}=HF^m\bar F^jR,\quad \beta_{m,j}=HF^m\bar F^jG,\quad \gamma_l=HF^lG,\quad \bar F=F-RH,\ \bar F^4=0.
```

- **R:** deadbeat 이득 $R=F^4\mathcal O^{-1}e_4$입니다. 이 값은 $M_u$의 마지막 열과 같습니다.
- **CT-CT:** 모든 계수를 암호화합니다. 각 $u(t)$는 ct-ct 곱셈 한 단(depth 1)으로 계산됩니다.
- **지연 $d$ (중요):** 기본값은 $d=K-n+1$ ($n=4$)입니다.
  - $d\le K-n+1$이면 입력창이 한 블록 구간 안에 있습니다. 그러면 블록 경계의 리셋 점프가 다음 창에 들어가지 않아 오차가 유계입니다.
  - $d>K-n+1$이면 점프가 $J_k=-F^{d-1}\bar F^{K-d+1}J_{k-1}$로 이어집니다. $\rho(F)>1$이고 $K$가 크면 발산합니다. 예를 들어 $K=d=64$에서 block당 6.8배입니다.
  - 이 경우 프로그램이 경고를 출력합니다.
- **허수부 제거:** bootstrapping 직후 $c+\overline{c}$로 허수부를 제거합니다 (기본값).
  - 필요한 conjugation key는 bootstrapping key에 포함되어 있습니다.
  - $\tfrac12$은 암호화된 $\alpha$에 흡수되므로 추가 레벨이 들지 않습니다.
- **초기 구간:** 첫 입력창이 채워지기 전에는 암호화된 $x_c(0)$로 $u(t)=HF^tx_c(0)+\sum_l\gamma_l y(t-1-l)$를 계산합니다. 이 구간의 계수는 쓰기 직전에 암호화합니다 (클라이언트 측, 서버 시간에서 제외).
- **Bootstrapping 정밀도:** 04와 같은 iterative bootstrapping (2회)을 사용합니다.
- **작은 값 처리:** $|v|<10^{-15}$인 값은 암호화 전에 0으로 둡니다. OpenFHE encoder의 `Scaling factor too small` 예외를 피하기 위해서입니다.

## 옵션

| 옵션 | 의미 |
|---|---|
| `--profile` | `stable`, `unstable_low`, `unstable_high` |
| `--K`, `--d` | 입력창 재사용 길이와 지연 (기본 $d=K-3$) |
| `--steps`, `--results-dir` | 실행 step 수, 출력 폴더 |
| `--pack-y` | 출력 패킹: $S$-slot zero-padded $y$, 블록마다 새로 시작하는 shift register $Y\leftarrow\mathrm{Rot}_{-1}(Y)+c_y$, 원소별 곱셈 1회 + $\log_2S$회 rotate-and-add, 공개 주기-4 mask로 $[s,0,0,0]$ 복원 |
| `--pack-u` | 입력창 패킹: sparse $u$ 4개를 rotation으로 $[u(\tau),\dots,u(\tau-3)]$로 모아 **창당 1회** bootstrapping. 입력창 항은 $\sum_j\mathrm{Rot}_j(W)\odot[\alpha_{m,j},0,0,0]$ |
| `--compact-inputs` | $y$와 패킹 계수는 depth−3, $\alpha$, $HF^t$, $x_c(0)$는 depth−2에서 암호화. 128-bit에서 메모리 때문에 필수 |
| `--secure128` | `HEStd_128_classic`. OpenFHE가 $N=65536$을 선택합니다. 없으면 $N=4096$ 연구용 파라미터 |
| `--keep-imag` | 허수부 제거를 끔 (발산 재현용) |
| `--boot-iters`, `--levels-after-bootstrap` | iterative bootstrapping 횟수 (기본 2), bootstrapping 후 레벨 (기본 2) |

## 빌드와 실행

```bash
eval "$(/opt/homebrew/bin/brew shellenv)"; source .venv/bin/activate; export CXX=g++-16
./05_openfhe_ctct_arx_window/scripts/build.sh
./05_openfhe_ctct_arx_window/scripts/run.sh unstable_low 16 200 OUT --pack-y --pack-u --compact-inputs   # N=4096 빠른 확인
```

`run.sh`의 인자는 `PROFILE K STEPS [OUT] [추가 옵션...]`이고, 끝나면 `tools/plot_basic_results.py`로 plot을 만듭니다.

## 출력

| 파일 | 내용 |
|---|---|
| `sample_trace.csv` | 평문 제어기 대비 $u$와 plant 오차, 같은 식을 평문으로 계산한 값과의 차이(step별 CKKS 오차), $u$의 inactive-slot leakage |
| `bootstrap.csv` | bootstrapping 작업별 deadline, 허용 여유(slack), 전후 오차, 허수부 오차 추정치 |
| `timing.csv` | step별 $y$ 암호화, $u$ 계산, register 갱신, 복호화, bootstrapping 시간과 ct-ct 곱셈/rotation 수 |
| `coefficients.csv` | 평문 계수 $\alpha,\beta,\gamma,\eta$ |
| `run.log`, `time.log` | 요약과 `/usr/bin/time -l` 결과 (메모리) |

`min_feasible_sampling_period_ms`는 실시간 스케줄 모델로 계산한 값입니다.
- **가정:** 서버는 $y(t)$를 받자마자 $u(t)$를 계산하고, 별도의 bootstrapping 작업자가 작업을 FIFO 순서로 처리합니다.
- **계산:** 이 조건에서 모든 deadline을 지키는 최소 sampling period를 측정된 연산 시간으로 구합니다.
- **한계:** 실제 병렬 실행이 아니라 측정 시간 기반 모델입니다.

비교 실험(N=4096 스칼라/패킹 비교, $K$ sweep, $d=K$ 발산 run 등)의 결과는 `../backup/05_results_archive/`에 있습니다.
