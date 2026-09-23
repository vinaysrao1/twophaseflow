# twophaseflow

Numerical models of incompressible oil–water flow in a horizontal pipe with an optional
ISO 5167-4 venturi. Reference case: 10 cm pipe, very light crude + water.

| Sub-project | Status | |
|---|---|---|
| [`onedim/`](onedim/README.md) — PipeFlow1D | implemented, tested | steady 1D mechanistic model: pressure profile, holdup, flow pattern, venturi Δp / C / permanent loss, sweeps over rate and water cut |
| `threed/` | planned | 3D LES + VOF on GPU (Julia), for flow visualization; runs on Modal |

The plan and literature survey are in [PLAN.md](PLAN.md).
