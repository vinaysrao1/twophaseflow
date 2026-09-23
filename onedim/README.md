# PipeFlow1D — steady 1D oil–water pipe flow

Sub-project 1 of `twophaseflow`: a fast, steady, one-dimensional model of incompressible
single-phase and oil–water flow in a horizontal pipe, including a classical (ISO 5167-4)
venturi. It is the engineering companion and cross-check for the 3D solver (sub-project 2).

Pure Julia, standard library only (no registry packages needed). Julia ≥ 1.10.

```bash
cd onedim
julia --project=. test/runtests.jl          # 842 tests, ~20 s
julia --project=. examples/single_phase.jl
julia --project=. examples/oilwater_sweep.jl
julia --project=. examples/oilwater_venturi.jl
```

```julia
using PipeFlow1D
sys = OilWater()                                   # light crude + water, σ = 0.025 N/m
v   = iso_venturi(0.10; β=0.5)                     # 10 cm pipe, ISO classical venturi
sol = solve(v, sys, Rates(Qtotal=m3h(40), wc=0.3)) # 40 m³/h, 30 % water cut
venturi_report(sol)                                # Δp, C, permanent loss
sol.p, sol.Hw, sol.pattern                         # profiles along x
write_csv("profile.csv", sol)
```

## What it computes

Inputs: pipe geometry (any sequence of cylinders and cones), wall roughness, fluid
properties, and either `(Qo, Qw)` or `(Qtotal, water cut)`. Rates helpers: `m3h`, `bpd`.

Outputs along the pipe: pressure `p(x)`, in-situ water holdup `Hw(x)` (differs from the
water cut when the phases slip), phase velocities, frictional gradient, flow pattern.
For venturis: upstream–throat Δp, effective discharge coefficient, permanent loss.

## Model

### Momentum balance (march along x)

Steady mixture momentum for a control volume of the pipe:

    dp/dx = −(1/A) d(A J)/dx + (dp/dx)_friction − (dp/dx)_local
    J = ρ_w α_w u_w² + ρ_o α_o u_o²          (momentum flux)

The acceleration term is discretised with the harmonic-mean area, which reproduces
Bernoulli exactly for single-phase and homogeneous flow (tested to round-off).
Local states are evaluated as fully developed at the local diameter (quasi-1D).

### Wall friction (`friction=`)

| Model | Notes |
|---|---|
| `Churchill()` (default) | Churchill (1977); one expression for laminar, transitional and turbulent flow; continuous in Re |
| `Colebrook()` | Colebrook–White, 64/Re below Re = 2300 |
| `Haaland()` | explicit approximation to Colebrook |
| `Laminar()`, `NoFriction()` | verification |

**Laminar–turbulent transition in 1D.** Churchill blends smoothly across Re ≈ 2000–4000.
The real transition is intermittent (puffs/slugs, Avila et al. 2011), so any 1D friction
factor in that band is uncertain by tens of per cent. For a 10 cm pipe this band means
velocities below ~0.04 m/s (water) or ~0.15 m/s (light crude), outside normal production.

### Venturi losses (`losses=`)

Convergent and throat: wall friction only (a venturi-quality contraction loses little;
the resulting C ≈ 0.99 matches ISO 5167-4). Divergent: Crane TP-410 loss coefficient
`K = 2.6 sin(θ/2)(1 − β²)²` on the throat velocity, spread over the diffuser length,
replacing wall friction there. For β = 0.5 and a 15° diffuser this gives ~22 % permanent
loss, at the high end of the 5–20 % typical of classical venturis; `CraneDiffuser(factor=…)`
calibrates it.

### Oil–water closures (`model=`)

- **`Homogeneous()`** — no slip (holdup = water cut), volume-averaged density, emulsion
  viscosity with the continuous phase chosen by the inversion model:
  - viscosity: `Brinkman()` (default), `KriegerDougherty()`, `Taylor()` (dilute);
  - inversion: `Arirachakaran()` (default: `wc_I = 0.5 − 0.1108 log10(μ_o/1 mPa·s)`, 0.447 for 3 mPa·s),
    `FixedInversion(wc)`.
  The emulsion viscosity peaks at inversion, producing the pressure-gradient peak seen
  experimentally.
- **`StratifiedTwoFluid()`** — two-fluid momentum balance for water under oil (Taitel &
  Dukler 1976; Brauner & Moalem Maron 1989). Solves one nonlinear equation for the
  interface position; all roots are found by scanning + Brent's method. Multiple roots
  (Ullmann et al. 2003) are reported in `nroots`; `root=:lowest|:highest` chooses. The
  faster phase sees the interface as a wall (hydraulic-diameter convention); where that
  switch brackets the solution the phases lock at equal velocity (`locked`).
- **`Mechanistic()`** (default) — predicts the pattern at each section and applies the
  stratified or homogeneous model. Dispersion criterion `HinzeBarnea()`: drops from
  turbulent break-up (Hinze 1955) must be smaller than the critical size for
  buoyancy-driven separation or deformation (Barnea 1986/87, as applied by Brauner 2001).
  Inside the venturi (`homogeneous_in_venturi=true`) the homogeneous model is used, as in
  multiphase venturi metering practice.

## Results for the reference case (D = 10 cm, light crude 820 kg/m³ / 3 mPa·s, water)

From `examples/`:

- Water, 1 m/s (Re 10⁵): f_D = 0.0179, venturi β = 0.5: Δp = 7.65 kPa, C = 0.990.
- Pattern: stratified below ~1.5 m/s mixture velocity at all water cuts; dispersed above
  (water-in-oil for wc < 0.447, oil-in-water above).
- Dispersed flow near inversion has up to ~1.5× the single-phase pressure gradient
  (Brinkman viscosity ~3.6× the oil viscosity at wc = 0.4).
- Light crude and water have similar viscosity, so slip is small in stratified flow
  (holdup within ~0.05 of the water cut).

## Validation (test/runtests.jl)

| Check | Result |
|---|---|
| Hagen–Poiseuille Δp and linear profile | exact (1e-10) |
| Colebrook implicit equation; Moody values at Re 10⁵, 10⁶ | 1e-10; 0.2 % |
| Churchill / Haaland vs Colebrook, Re 10⁴–10⁷, ε/D 0–10⁻³ | within 3 % |
| Frictionless venturi vs Bernoulli at every node, full pressure recovery | round-off |
| Venturi C with friction, Re_D 10⁵–5×10⁵ | 0.975–1.0 (ISO 0.984–0.995) |
| Venturi grid convergence (dx halved) | < 10⁻⁴ relative |
| Stratified: mass conservation, equal dp/dx from both phase equations | 10⁻¹⁰, 10⁻⁶ |
| Stratified: viscous oil accumulates; equal fluids → holdup ≈ water cut | ✓ |
| Vanishing phase → single-phase dp/dx | within 2 % |
| Emulsion viscosity peak at inversion; Taylor/Krieger → Einstein 2.5 | ✓ |
| Pattern monotone in velocity, transition 0.5–3 m/s | ✓ |

Not yet done: comparison with published oil–water datasets (Trallero 1995, Elseth 2001,
Rodriguez & Oliemans 2006, Angeli & Hewitt 2000). Those data are mostly in figures and
need digitising; this is the next validation step.

## Limitations

- Accuracy of any 1D oil–water model is ±20–40 % on pressure gradient, worse near
  phase inversion (Chakrabarti et al. 2005; Rodriguez & Oliemans 2006).
- The pattern map resolves stratified vs dispersed only. Intermediate patterns (stratified
  with mixing, dual dispersions) are assigned to the nearer limit, so dp/dx jumps at the
  transition, whereas real transitions are gradual.
- The hydraulic-diameter stratified model can predict slightly lower dp/dx than single
  phase water at high water cut (the water "sees" the interface as a free surface): a
  known artefact of the approach.
- The Hinze drop size is the dilute-limit value; `HinzeBarnea(dense_factor=…)` adds a
  simple concentration correction.
- Phase-inversion correlations are uncertain and show hysteresis in experiments.
- Isothermal, Newtonian, no gas, hydrostatic level-gradient terms neglected (so no
  hydraulic-jump behaviour of stratified flow through the venturi).
- Steady state only.

## References

- Agrawal, Gregory & Govier (1973), Can. J. Chem. Eng. 51:280.
- Angeli & Hewitt (2000), Int. J. Multiphase Flow 26:1117.
- Arirachakaran et al. (1989), SPE 18836.
- Avila et al. (2011), Science 333:192.
- Barnea (1986), Int. J. Multiphase Flow 12:733; Barnea (1987), Int. J. Multiphase Flow 13:1.
- Brauner (2001), Int. J. Multiphase Flow 27:885; Brauner (2003), *Liquid–liquid two-phase flow systems*, CISM.
- Brauner & Moalem Maron (1989), Int. J. Multiphase Flow 15; (1992) 18:103.
- Brinkman (1952), J. Chem. Phys. 20:571.
- Churchill (1977), Chem. Eng. 84(24):91.
- Colebrook (1939), J. ICE 11:133; Haaland (1983), J. Fluids Eng. 105:89.
- Crane Co., *Flow of Fluids*, Technical Paper 410.
- Hinze (1955), AIChE J. 1:289.
- ISO 5167-4, Measurement of fluid flow by pressure differential devices — Venturi tubes.
- Krieger & Dougherty (1959), Trans. Soc. Rheol. 3:137.
- Rodriguez & Oliemans (2006), Int. J. Multiphase Flow 32:323.
- Taitel & Dukler (1976), AIChE J. 22:47.
- Trallero, Sarica & Brill (1997), SPE Prod. & Facilities 12:165.
- Ullmann, Zamir, Gat & Brauner (2003), Int. J. Multiphase Flow 29:1583.
