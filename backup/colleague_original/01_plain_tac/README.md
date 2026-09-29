# 01 — Plain TAC reformulation

Purpose: verify the TAC algebra in ordinary floating-point arithmetic before introducing CKKS.

Run:

```bash
python3 tac_plain.py --profile stable
python3 tac_plain.py --profile unstable_low
python3 tac_plain.py --profile unstable_high
```

The script chooses $R$ by dual pole placement so that $F-RH$ has poles `0.2, 0.3, 0.4, 0.5` by default. Override them with `--tac-poles`.

Expected result: `u_error` and plant-state error should remain at floating-point roundoff level because

```math
(F-RH)x+Gy+Ru = Fx+Gy
```

whenever $u=Hx$.
