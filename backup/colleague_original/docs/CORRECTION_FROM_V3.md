# v3 -> v4 구현 수정 사항

이전 v3 OpenFHE 코드는 control/measurement scalar를 `[s,s,s,s]` 형태로 replicated packing하여 사용했습니다.

현재 의도한 구조는 다음과 같습니다.

```math
c_{y_i}\leftrightarrow[y(t+i),0,0,0],
\qquad
c_{u_i}\leftrightarrow[u(t+i),0,0,0].
```

또한 각 시점의 $u$와 $y$는 서로 다른 ciphertext입니다.

따라서 v4에서는 다음을 수정했습니다.

1. PT-CT의 $P_i x$ 계산을 rotation + slot-0 one-hot coefficient 방식으로 변경.
2. CT-CT의 $P_i$ coefficient를 16개의 encrypted one-hot mask로 변경.
3. 모든 $y(t+i)$를 별도의 sparse scalar ciphertext로 암호화.
4. 모든 $u(t+i)$를 별도의 sparse scalar ciphertext로 계산 및 bootstrap.
5. $M_uU$와 $M_yY$ reconstruction을 sparse-scalar scatter 방식으로 변경.
6. `slot_spread` diagnostic을 제거하고 `active-slot error`와 `inactive-slot leakage`로 변경.
7. replicated-scalar architecture를 사용했던 이전 advanced result archive를 현재 공유용 repo에서 제거.

이 수정 이후 03/04의 장기 실험은 다시 수행해야 합니다.
