# backup

현재 작업(`05_openfhe_ctct_arx_window/`, `tex/`)에 직접 필요하지 않은 코드와 결과를 보관합니다.
지운 것은 없고, git으로 추적되던 파일은 `git mv`로 옮겼으므로 이력이 남아 있습니다.

| 경로 | 내용 |
|---|---|
| `colleague_original/` | 동료에게 받은 원래 저장소: 01 평문 TAC, 02 평문 u/y reconstruction, 03 PT-CT, 04 CT-CT, `docs/`, `experiments/` (v3 archive), `scripts/`, 원래 `README.md`/`README_EN.md`, `tools/verify_sparse_layout.py` |
| `05_results_archive/` | 05의 비교 실험 결과 (gitignore 대상) |
| `analysis_scripts/` | 분석 및 그림 생성 스크립트, 진단용 코드 |
| `misc/less_help_accidental.txt` | 원래 저장소 루트에 실수로 커밋되어 있던 `less` 도움말 파일 (`truction \`) |

## colleague_original에서 우리가 바꾼 것

- `03_openfhe_ptct_uy_reconstruction/src/uy_reconstruction_ptct.cpp`: `--clear-imag` 옵션 추가 (bootstrapping 직후 허수부 제거, $\tfrac12$은 $M_u$ mask에 흡수).
- `03_openfhe_ptct_uy_reconstruction/scripts/run.sh`: 4번째 이후 인자를 실행 파일에 전달.
- 결과 폴더: `03_*/results/`에 허수부 분석 run이 있습니다.
  - `<profile>_500`: 500 step, bootstrapping만
  - `<profile>_200_clearimag`: 허수부 제거
- 결과 폴더: `04_*/results/`에도 500 step run이 있습니다.

**주의:** 옮긴 뒤로 01–04의 build/run 스크립트는 저장소 루트 기준 상대 경로(`include/`, `tools/`, `config/`)가 맞지 않습니다. 다시 돌리려면 저장소 루트로 복사하거나 경로를 고쳐야 합니다.

## 05_results_archive

| 이름 패턴 | 설정 |
|---|---|
| `<profile>_K<K>` | $N=4096$, 스칼라 방식, $d=K$ |
| `<profile>_K<K>_pack` | 출력 패킹 |
| `<profile>_K<K>_packyu` | 출력 + 입력창 패킹 |
| `unstable_low_K4_d1` | 기존 구조와 같은 지연 ($d=1$) |
| `unstable_low_K4_keepimag` | 허수부 제거 없음 (오차 증가 재현) |
| `*_500` | 500 step 유계성 확인 |
| `<profile>_K<K>_s128` | 128-bit, $d=K$ ($K=16,32,64$ sweep) |
| `unstable_high_K64_s128_1000` | 128-bit, $K=d=64$, 1000 step: **리셋 점프 누적으로 발산** (block당 6.8배) |
| `summary.csv` | 위 run들의 요약 (`analysis_scripts/analyze05.py` 출력) |

## analysis_scripts

| 파일 | 용도 |
|---|---|
| `analyze05.py` | 05 비교 run 요약표와 `tex/figures/arx_packyu_u_error.png`, `arx_min_Ts.png` 생성 |
| `final_d_stability_plot.py` | `tex/figures/final_d_stability.png` 생성 ($d=64$ vs $d=61$) |
| `arx_window_check.py` | 지연 입력창 ARX 식이 원래 제어기와 같은지 평문 확인 |
| `arx_noise_check.py`, `arx_noise_check_d.py` | 평문 + 노이즈로 $K$, $d$별 오차 안정성 확인 (발산 원인 규명) |
| `u_error_500.py` | 03/04 500 step $u$ 오차 그림 |
| `ptct_imag.cpp`, `build_imag.sh` | 03에 허수부 오차 추정치(`GetLogError`) 로깅을 넣은 진단용 복사본 |
| `secparams.cpp` | 현재 파라미터의 $\log QP$와 128-bit에 필요한 ring dimension 계산 |

`ptct_imag.cpp`, `build_imag.sh`, `u_error_500.py`는 작성 당시 경로(스크래치 폴더, 03/04 원래 위치)를 기준으로 합니다. 다시 쓰려면 경로를 맞춰야 합니다.
