# Advanced experiments

현재 공유용 v4에서는 이전 long-horizon 결과를 포함하지 않습니다.

이유는 이전 결과가 `[u,u,u,u]`, `[y,y,y,y]` replicated-scalar OpenFHE 구현을 기반으로 생성되었기 때문입니다. 현재 03/04의 의도한 architecture는

```math
c_{u_i}\leftrightarrow[u(t+i),0,0,0],
\qquad
c_{y_i}\leftrightarrow[y(t+i),0,0,0]
```

이며 각 시점의 $u,y$는 별도 ciphertext입니다.

따라서 다음 실험은 v4 corrected implementation으로 다시 수행한 뒤 이 디렉터리에 추가해야 합니다.

- long-horizon PT-CT / CT-CT,
- same-measurement controller-error dynamics,
- post-bootstrap level/tower sweep,
- active-slot error / inactive-slot leakage,
- bootstrap production-vs-fresh comparison.

이전 구현과 수정 내용은 [`../docs/CORRECTION_FROM_V3.md`](../docs/CORRECTION_FROM_V3.md)를 참고하세요.
