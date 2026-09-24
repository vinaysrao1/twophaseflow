# Preliminary checks of candidate packages (2026-09-24)

Run on this repository's development container (4 CPU cores, no GPU), Julia 1.13,
WaterLily.jl v1.8.0 (commit of 2026-09-21), InterfaceAdvection.jl 1.0.0-DEV (2026-09-14).

```bash
julia --project=<env with WaterLily> -t 4 pipe_laminar.jl
julia --project=<env with WaterLily> -t 4 throughput.jl
```

## 1. Laminar pipe flow with WaterLily's immersed boundary (BDIM) — `pipe_laminar.jl`

Periodic pipe, body-force driven, steady Poiseuille flow; exact u_max = 1.

| R [cells] | Q / Q_exact | u_max | effective radius − R |
|---|---|---|---|
| 8  | 1.161 | 1.068 | +0.30 cells |
| 16 | 1.082 | 1.038 | +0.32 cells |
| 32 | 1.041 | 1.020 | +0.32 cells |

The wall sits ~0.32 cells outside the geometric radius: first-order, systematic, and so
correctable by shifting the signed-distance function. Uncorrected at D/100 (R = 50) the
laminar flow-rate error would be ≈ 2.6 %.

## 2. Throughput — `throughput.jl` (single phase, Float32, pipe body, 4 threads)

| grid | cells | cell-steps/s |
|---|---|---|
| 96×48×48 | 2.2e5 | 3.1e6 |
| 192×96×96 | 1.8e6 | 7.6e6 |

Grid sizes must be a·2ⁿ with n > 2 (multigrid requirement).

## 3. InterfaceAdvection.jl test suite

11 test sets. One Float32 test in PLIC.jl uses exact `==` on results that differ in the
last bit (0.20147918 vs 0.20147927) and aborts the run; with that comparison changed to `≈`
all 136 tests pass. Coverage is unit-level (6 flow tests); no tests involve bodies or
walls, since bodies are not supported.
