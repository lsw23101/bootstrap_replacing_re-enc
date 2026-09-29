# Bootstrapping으로 re-encryption 대체하기: CKKS 동적 제어기

불안정한 내부 행렬 $F$를 가질 수 있는 선형 동적 제어기를 CKKS(OpenFHE)로 무한 시간 운용하는 연구 코드입니다.
TAC(Kim, [arXiv:2307.16445](https://arxiv.org/abs/2307.16445))의 re-encryption을 **bootstrapping**으로 대체합니다.
ARX 표현으로 상태 암호문 없이 입출력만으로 제어 입력 $u(t)$를 계산합니다.

## 먼저 읽을 것

1. **[`tex/formulation.pdf`](tex/formulation.pdf):** 전체 정리입니다. 작업 경과 → 문제 설정 → TAC/ARX 유도 → 지연 입력창 → 블록 리셋 점프와 안정 조건 → CKKS 구현 → 허수부 문제 → 실험 결과 순서입니다.
2. **[`CONTEXT.md`](CONTEXT.md):** 작업 대화 흐름입니다. 어떤 질문에서 무엇을 발견했고, 무엇을 결정했는지 시간 순서로 정리했습니다.
3. [`tex/issue.pdf`](tex/issue.pdf): 허수부 누적 문제만 다룬 단독 보고서입니다.
4. [`CLAUDE.md`](CLAUDE.md): Claude Code로 이어서 작업할 때 쓰는 규칙과 핵심 사실입니다 (자동으로 읽힙니다).

## 작업 흐름 한눈에 보기

```
동료 저장소 (01 평문 TAC → 02 평문 u/y reconstruction → 03 PT-CT → 04 CT-CT)
   │  500 step으로 늘리니 unstable F에서 u 오차가 ρ(F)^4 배/block으로 발산
   ▼
원인: bootstrapping이 CKKS 허수부를 보존 → plant는 실수부만 받음 → 허수부가 F로 개루프 전파
   │  해결: bootstrapping 직후 c + conj(c) (레벨 소모 없음)
   ▼
구조 재정립: u/y reconstruction = deadbeat TAC의 ARX 형태
   │  목표: bootstrapping된 u 입력창을 K step 재사용해서 bootstrapping 시간 확보
   ▼
05: 지연 입력창 ARX, CT-CT
   │  y 패킹 (zero-padded shift register + rotate-sum)  → u 계산이 K에 대해 로그로만 증가
   │  입력창 패킹 (창당 bootstrapping 1회)            → bootstrapping 1/K 회/step
   ▼
128-bit (N=65536): bootstrapping ≈ 6 s, u 계산 ≈ 0.2 s/step → K=64 선택
   │  K=d=64에서 블록 리셋 점프가 누적되어 발산 발견 (평문에서도 동일)
   ▼
안정 조건 d ≤ K − n + 1  →  최종 설정 K=64, d=61: 오차 ≤ 1.1e-7 (유계), 여유 ≈ 18 s vs 필요 ≈ 6 s
```

## 저장소 구성

| 경로 | 내용 |
|---|---|
| `05_openfhe_ctct_arx_window/` | **현재 구현** (C++/OpenFHE). `results/`에 최종 설정 검증 run 2개 (128-bit, K=64, d=61, 640 step) |
| `tex/` | `formulation.tex/.pdf` (전체 정리), `issue.tex/.pdf` (허수부 보고서), `figures/` |
| `config/`, `include/` | plant/제어기 profile 정의(Python)와 생성된 C++ 헤더 |
| `tools/` | profile 생성과 검증, 결과 plot |
| `backup/colleague_original/` | 동료 원본 01–04와 문서. 03에는 `--clear-imag` 수정이 들어 있음 |
| `backup/analysis_scripts/` | 문서 그림 생성과 진단 스크립트 |
| `backup/05_results_archive/` | 비교 실험 결과. 저장소에는 `summary.csv`만 올라감 |

## 재현

```bash
# macOS (Apple Silicon) 기준
brew install gcc cmake
git clone https://github.com/openfheorg/openfhe-development ~/openfhe-development
cd ~/openfhe-development && git checkout 1306d14f   # v1.5.1
mkdir build && cd build && CC=gcc-16 CXX=g++-16 cmake .. -DBUILD_UNITTESTS=OFF -DBUILD_EXAMPLES=OFF -DBUILD_BENCHMARKS=OFF && make -j
cd -   # 이 저장소로

python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt
export CXX=g++-16
./05_openfhe_ctct_arx_window/scripts/build.sh

# 빠른 확인 (N=4096 연구용 파라미터, 수십 초)
./05_openfhe_ctct_arx_window/scripts/run.sh unstable_low 16 200 OUT --pack-y --pack-u --compact-inputs

# 최종 설정 (128-bit, 약 4분, 메모리 약 3.8 GB)
./05_openfhe_ctct_arx_window/arx_window_ctct --profile unstable_low --K 64 --steps 640 \
  --results-dir OUT --pack-y --pack-u --compact-inputs --secure128
```

옵션과 출력 파일 설명은 [`05_openfhe_ctct_arx_window/README.md`](05_openfhe_ctct_arx_window/README.md)에 있습니다.

## 주요 결과

| 항목 | 값 |
|---|---|
| 허수부 제거 전 (03, 500 step) | block당 $\rho(F)^4\approx1.197$배로 발산, 360 step대에서 복호화 실패 |
| 허수부 제거 후 (03) | $u$ 오차 약 $10^{-8}$에서 유지 |
| 05, $N=4096$, 패킹 | $u$ 오차 $\le10^{-12}$, 최소 $T_s$ 약 60–100 ms (모델) |
| 05, 128-bit, $K=d=64$ | 블록 리셋 점프 누적으로 발산 (block당 6.8배) |
| **05, 128-bit, $K=64$, $d=61$ (최종)** | $u$ 오차 $\le1.1\times10^{-7}$ 유계, bootstrapping 약 6 s, $u$ 계산 약 0.2 s/step, 최소 $T_s$ 약 0.3 s (모델) |

최소 $T_s$는 측정 시간 기반 스케줄 모델 값입니다. bootstrapping 작업자가 별도 자원에서 병렬로 돈다고 가정합니다.
실제 병렬 실행으로 검증하는 것은 남은 과제입니다.
