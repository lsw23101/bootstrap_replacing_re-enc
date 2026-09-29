# mode_1000 long-horizon run

Requested: 1000 blocks.

Observed:

- 149 full blocks completed.
- 598 samples completed before failure.
- failure during block 149, sample 2.
- OpenFHE reported excessive CKKS approximation error during decode.
- measured same-measurement internal-error growth over blocks 20+: approximately 1.1973239/block.
- theoretical factor: rho(F)^4 = 1.1972807715668585.
- maximum recorded nominal-control error: approximately 1.66096e-1.
- maximum recorded plant-state infinity-norm error: approximately 1.74812e-2.
- maximum recorded production bootstrap-added error: approximately 3.80062e-1.

See `../../docs/LONG_HORIZON_1000.md` for interpretation.
