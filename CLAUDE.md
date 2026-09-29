# CLAUDE.md

CKKS 동형암호로 선형 동적 제어기(내부 행렬 $F$가 불안정할 수 있음)를 무한 시간 동안 운용하는 연구 코드입니다.
동료 저장소(01–04, 지금은 `backup/colleague_original/`)에서 출발했고, 현재 작업은 `05_openfhe_ctct_arx_window/`와 `tex/`입니다.
전체 유도와 결과, 작업 경과는 `tex/formulation.tex`(`formulation.pdf`)에 있습니다.

## 목표

- **re-encryption 대신 bootstrapping:** 재사용할 $u$ 암호문의 곱셈 레벨을 bootstrapping으로 회복합니다.
- **TAC + ARX:** 상태 암호문 없이 입출력만으로 $u(t)$를 계산합니다. TAC는 Kim, arXiv 2307.16445입니다 (로컬에 `24_Kim_TAC_Further.pdf`가 있을 수 있으나 저장소에는 올리지 않음).
- **K스텝 재사용:** bootstrapping된 입력창을 $K$ step 재사용해서 bootstrapping 시간을 확보합니다.
- **CT-CT:** 제어기 계수까지 암호화합니다.
- **판단 기준:** 오차 크기보다 **발산 여부**입니다. $10^{-5}$ 수준의 오차는 허용합니다 (사용자 기준).

## 디렉터리

| 경로 | 내용 |
|---|---|
| `05_openfhe_ctct_arx_window/` | **현재 구현.** `src/arx_window_ctct.cpp`, `scripts/`, `README.md`, `results/` (최종 설정 검증 run 2개) |
| `tex/formulation.tex` | 전체 정리: 작업 경과 → 문제 설정 → TAC/ARX → 지연 입력창 → $d$ 안정 조건 → CKKS 구현 → 허수부 문제 → 실험 결과 |
| `tex/issue.tex` | 허수부 누적 문제 단독 보고서 (03 기준) |
| `tex/figures/` | 문서 그림. 생성 스크립트는 `backup/analysis_scripts/` |
| `config/controller_profiles.py` | plant $(A,B,C)$와 제어기 profile 원본. 수정하면 `python tools/generate_controller_profiles.py`로 `include/controller_profiles.h`를 다시 생성 |
| `include/controller_profiles.h` | 05가 include하는 생성 헤더 (A,B,C,F,G,H,$M_u$,…) |
| `tools/` | `generate_controller_profiles.py`, `validate_profiles.py`, `plot_basic_results.py` (05 `run.sh`가 사용) |
| `backup/` | 동료 원본 01–04, 05 비교 실험 결과, 분석 스크립트. 자세한 내용은 `backup/README.md` |

## 환경과 빌드 (macOS, Apple Silicon, 16 GB)

```bash
eval "$(/opt/homebrew/bin/brew shellenv)"
source .venv/bin/activate          # run.sh가 plot에 python3를 쓰므로 필요
export CXX=g++-16                  # Apple clang은 -fopenmp 미지원
./05_openfhe_ctct_arx_window/scripts/build.sh
```

- **Python:** `.venv`에 `requirements.txt`가 설치되어 있습니다 (시스템 Python 3.9).
- **OpenFHE:** v1.5.1 (commit `1306d14f`)을 `~/openfhe-development`에 받아 in-tree `build/`로 빌드했습니다 (gcc-16, OpenMP on).
- **LaTeX:** `tex/`에서 `xelatex formulation.tex`를 두 번 실행합니다 (kotex, AppleMyungjo). 빌드 부산물(aux/log/out/toc)은 지웁니다.
- **실행 중 빌드 금지:** 실행 중인 `arx_window_ctct`가 있을 때 다시 빌드하지 마세요 (바이너리 덮어쓰기).

## 05 실행

```bash
# 기준 설정: 128-bit, K=64, d=61(기본값 K-3), 출력/입력창 패킹, 입력 레벨 축소
05_openfhe_ctct_arx_window/arx_window_ctct --profile unstable_low --K 64 --steps 640 \
  --results-dir OUT --pack-y --pack-u --compact-inputs --secure128      # 약 4분, 3.8 GB

# 빠른 확인 (N=4096 연구용, 수십 초)
./05_openfhe_ctct_arx_window/scripts/run.sh unstable_low 16 200 OUT --pack-y --pack-u --compact-inputs
```

옵션 표와 출력 파일 설명은 `05_openfhe_ctct_arx_window/README.md`에 있습니다.
`run.log`의 `min_feasible_sampling_period_ms`는 측정 시간 기반 스케줄 모델 값입니다. bootstrapping 작업자가 별도 자원에서 병렬로 돈다고 가정합니다.

## 핵심 사실 (코드를 바꿀 때 지켜야 할 것)

1. **허수부 제거는 필수입니다.**
   - bootstrapping은 CKKS slot의 허수부를 보존합니다. 그런데 plant는 복호화된 실수부만 받습니다.
   - 허수부를 남기면 그 성분이 plant를 거치지 않고 $F$로 개루프 전파되어, $\rho(F)>1$에서 발산합니다.
   - bootstrapping 직후 $c+\overline{c}$를 계산합니다. conjugation key(index $M-1$)는 bootstrapping key에 포함되어 있습니다.
   - $\tfrac12$은 $\alpha$에 흡수합니다. 이 때문에 곱셈 레벨을 쓰지 않습니다.
2. **지연 $d\le K-n+1$ ($n=4$) 필수입니다.**
   - 이 조건을 넘으면 입력창이 이전 블록의 리셋 시점을 포함합니다. 그러면 리셋 점프가 $J_k=-F^{d-1}\bar F^{K-d+1}J_{k-1}$로 누적됩니다.
   - 예: $K=d=64$에서 block당 6.8배로 발산합니다. 평문에서도 같은 값이 나오므로 CKKS와 무관한 구조 문제입니다.
   - 기본값은 $d=K-3$이고, 조건을 넘는 $d$를 주면 경고를 출력합니다.
3. **bootstrapping 대상 $u$는 4-slot sparse `[u,0,0,0]`를 유지합니다.** 사용자 요청이고, 입력창 패킹의 전제입니다.
   - $y$는 $S$-slot 제로패딩으로 패킹합니다.
   - 주기-4 mask로 결과를 sparse로 되돌립니다. $y$ 항은 fresh 암호문에서 출발하므로 이 mask 레벨은 bootstrapping 예산과 무관합니다.
4. **오래 사는 rotation 누적 암호문을 만들지 않습니다.** $y$ register는 block마다 새로 시작합니다. 그래야 key-switching 노이즈가 유계입니다.
5. **레벨 예산:** bootstrapping 후 2 레벨(ct-ct 곱셈 1 + 여유 1)입니다. `--compact-inputs`에서는 $y$와 패킹 계수를 depth−3, $\alpha$, $HF^t$, $x_c(0)$를 depth−2에서 암호화합니다. 128-bit에서는 메모리 때문에 필수입니다.
6. **작은 값 인코딩:** $0<|v|\cdot2^{59}<1$이면 OpenFHE가 `Scaling factor too small` 예외를 던집니다. 05는 $|v|<10^{-15}$를 0으로 둡니다.
7. **기존 reconstruction의 정체:** $M_u, M_y$ reconstruction은 deadbeat TAC ($R=F^n\mathcal O^{-1}e_n$, $R$ = $M_u$의 마지막 열)의 ARX 형태와 같습니다.
8. **realization 불변:** ARX 계수는 좌표변환에 불변입니다. `unstable_low`와 `unstable_high`의 계수가 동일합니다.
9. **계수 증폭:** 계수는 $\rho(F)^{2K-2}$로 커집니다 ($K=64$: $\max|\alpha|=187$, $\max|\beta|=319$). 128-bit의 낮은 bootstrapping 정밀도와 곱해져 오차가 $10^{-7}$ 수준이 됩니다.

## 현재 설정과 결과 (2026-09-29)

- **기준 설정: 128-bit ($N=65536$), $K=64$, $d=61$, `--pack-y --pack-u --compact-inputs`.**
  - 640 step 최대 $u$ 오차: `unstable_high` $1.1\times10^{-7}$, `unstable_low` $3.3\times10^{-8}$. 모두 유계입니다.
  - bootstrapping 1회 약 6–7 s, $u$ 계산 약 0.2 s/step, 최소 $T_s$ 약 0.3 s (모델)입니다.
  - bootstrapping 여유는 약 3배이고, 메모리는 약 3.8 GB입니다.
- **$K$ 선택 기준:** $K\gtrsim T_{bt}/T_u$ ($\approx34$)이면 bootstrapping이 병목에서 빠집니다. $K=32$는 여유가 1.2배라 bootstrapping 시간 변동에 취약해서 $K=64$로 정했습니다.
- **N=4096** (연구용, 128-bit 아님): 오차 $\le10^{-12}$입니다. 비교 결과는 `backup/05_results_archive/`에 있습니다.

## 남은 과제

- 실제 병렬 실행(bootstrapping 작업자 분리)으로 스케줄 모델을 검증합니다.
- 128-bit 최종 설정을 수천 step으로 확인합니다.
- 초기 구간의 계산 최댓값을 줄입니다.
- 복호화 noise flooding 모드(`NOISE_FLOODING_DECRYPT`)와 허수부 제거의 상호작용을 정리합니다. 허수부를 지우면 기본 모드의 오차 추정이 작아집니다.
- 4-tank plant로 확장합니다 (SISO/MIMO, 선형/비선형 미정, 사용자가 보류함).

## 작업 규칙

- 새 기능은 05(또는 새 단계)에 추가합니다. `backup/`은 참고용입니다.
- 비교나 테스트용으로 돌린 결과는 `backup/05_results_archive/`로 옮기고, `05_*/results/`에는 최신 기준 결과만 남깁니다.
- `05_*/results/`(최신 기준 run)는 git으로 추적합니다. 비교 결과(`backup/05_results_archive/*/`)와 바이너리는 gitignore 대상입니다.
- 작업 대화 흐름은 `CONTEXT.md`에 있습니다. 큰 결정이 생기면 이어서 기록합니다.
- 커밋과 푸시는 사용자가 요청할 때만 합니다. 원격: `origin` = 동료 저장소, `bootstrap` = [lsw23101/bootstrap_replacing_re-enc](https://github.com/lsw23101/bootstrap_replacing_re-enc).
