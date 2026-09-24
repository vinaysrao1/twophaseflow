# twophaseflow — Plan (draft for review)

Status: **v3**. Decided: Julia; D = 10 cm; very light crude + water; ISO 5167-4 classical venturi; full 3D flow visualization as the end goal.

Two sub-projects:
1. **`onedim/` — PipeFlow1D** (steady 1D mechanistic model). **Implemented** — see `onedim/README.md`.
2. **`threed/` — 3D LES/VOF solver** (sections 3–9 below). Not started. Production runs on Modal (GPU); deployment details to be worked out.

## 1. Scope

Numerical model of incompressible, fully liquid-filled ("full pipe", no gas) flow in a **horizontal** pipe.

| Stage | Physics | Geometry |
|---|---|---|
| 1 | Single phase, incompressible, isothermal | Straight circular pipe |
| 2 | Single phase | Pipe with venturi (converging → throat → diverging) |
| 3 | Two immiscible liquids (crude oil + water), any water cut, any total rate | Straight pipe |
| 4 | Oil + water | Pipe with venturi |
| 5 (optional) | Transient (time-varying production rates) | Any of the above |

Inputs throughout: pipe diameter `D`, length `L`, wall roughness `ε`, area profile `A(x)`, fluid properties (ρ, μ, interfacial tension σ), and either (`Q_oil`, `Q_water`) or (`Q_total`, water cut `WC = Q_w/Q_total`).
Outputs: pressure profile `p(x)`, pressure gradient, in-situ water holdup `H_w(x)` (≠ input water cut when phases slip), phase velocities, flow pattern, venturi differential pressure and permanent loss.

## 2. Research survey

### 2.1 Single-phase pipe flow (well established)
- **Laminar** (Re < ~2300): Hagen–Poiseuille, `f = 64/Re` (Darcy). Exact.
- **Turbulent**: Darcy–Weisbach `Δp = f (L/D) ρV²/2` with Colebrook–White (implicit, the standard), explicit approximations by Haaland (1983) and Swamee–Jain (1976). Churchill (1977) gives one continuous expression over laminar/transition/turbulent — useful for avoiding discontinuities in a solver. Validation: Moody chart, Nikuradse data.
- **Developing flow**: entrance length ≈ 0.05·Re·D (laminar), ≈ 10–60 D (turbulent); Shah & London (1978) apparent friction factors. Usually negligible for long pipes; relevant immediately downstream of a venturi.
- **Full CFD** (resolving the velocity field): finite-volume incompressible Navier–Stokes with pressure–velocity coupling — SIMPLE (Patankar & Spalding 1972), PISO (Issa 1986), projection (Chorin 1968). Turbulence needs RANS (k-ε, k-ω SST) or LES. Mature open-source: OpenFOAM (`simpleFoam`/`pimpleFoam`).
- **Transients** in incompressible flow: rigid-column model, `ρ L_eff dV/dt = Δp_driving − losses`. (True water hammer requires compressibility / method of characteristics, Wylie & Streeter — out of scope given the incompressible assumption.)

### 2.2 Venturi (single phase)
- 1D Bernoulli + continuity: `Q = C_d A_t √(2Δp / (ρ(1−β⁴)))`, `β = d_t/D`.
- Discharge coefficient and geometry standardized in **ISO 5167-4** (classical venturi: convergent 21°±1°, divergent 7–15°, `C_d ≈ 0.984–0.995` for Re_D > 2×10⁵; lower and Re-dependent at low Re — relevant for viscous crude).
- Permanent (unrecovered) pressure loss ≈ 5–20 % of Δp; diffuser losses from Idelchik (Handbook of Hydraulic Resistance) / Crane TP-410 loss coefficients, which depend on divergence angle and area ratio. Diffusers above ~7–8° half-angle risk separation.

### 2.3 Oil–water flow in horizontal pipes
Key reviews: **Brauner (2003)** "Liquid–liquid two-phase flow systems" (the standard reference); Angeli & Hewitt (2000); Ismail et al. (2015); Shi & Yeung and later reviews. Unlike gas–liquid, density ratio is ~1 and viscosity ratio can be 1–10⁴, so interfacial tension, wettability and oil viscosity dominate.

**Flow patterns** (Trallero, Sarica & Brill 1997 classification, widely adopted):
- ST — stratified (smooth)
- ST&MI — stratified with mixing at interface
- Do/w & w — dispersion of oil in water over a water layer
- Do/w — oil-in-water dispersion (emulsion), water continuous
- Dw/o & Do/w — dual dispersion
- Dw/o — water-in-oil dispersion, oil continuous
- Dual-continuous (Lovick & Angeli 2004); annular / core-annular for heavy oils (Joseph et al. 1997, Bannwart 2001).

Pattern transitions: Brauner & Moalem Maron (1992) (stability of stratified flow, Kelvin–Helmholtz type), Trallero (1997) mechanistic model; droplet-size-based dispersion criteria from Hinze (1955), Brauner (2001).

**Stratified flow — two-fluid model** (Taitel & Dukler 1976 for gas–liquid; adapted to liquid–liquid by Brauner & Moalem Maron 1989, Trallero 1995, Rodriguez & Oliemans 2006). Steady, fully developed, horizontal; momentum balance per phase:

```
−A_w dp/dx − τ_w S_w + τ_i S_i = 0
−A_o dp/dx − τ_o S_o − τ_i S_i = 0
```
Eliminating `dp/dx` gives one nonlinear equation in the interface height `h/D`, solved by root-finding; then `dp/dx` and holdup follow. Closure choices: wall shear via phase hydraulic diameters (Agrawal et al. 1973 convention: the faster phase "sees" the interface as a wall), interfacial friction factor (smooth vs. wavy; Ullmann & Brauner 2006 correlations). Known issue: **multiple holdup solutions** exist in some ranges (Ullmann, Zamir, Gat & Brauner 2003); must be detected and reported, not hidden.

**Dispersed flow — homogeneous (no-slip) model**: `H_w = WC`, mixture density `ρ_m = WC ρ_w + (1−WC) ρ_o`, mixture viscosity from an emulsion law, then single-phase friction correlation. Emulsion viscosity options: Einstein (1906, dilute), Brinkman (1952), Pal & Rhodes (1989), Krieger–Dougherty (1959). Emulsion viscosity peaks near phase inversion — a major pressure-drop effect.

**Phase inversion** (switch water-continuous ↔ oil-continuous): Arirachakaran et al. (1989) empirical correlation (inversion water cut falls as oil viscosity rises); Brauner & Ullmann (2002) energy-minimization model; Yeh et al. (1964). Accuracy is limited; hysteresis is observed experimentally.

**Slip in intermediate patterns**: drift-flux (Zuber & Findlay 1965) or two-fluid with entrainment; mechanistic models by Trallero (1995), Al-Wahaibi & Angeli.

**Expected accuracy** of published 1D models vs. data: typically ±20–40 % on pressure gradient, larger near inversion (Chakrabarti et al. 2005 reported 40–200 % spread; Rodriguez & Oliemans ~35 %; Ullmann & Brauner wavy-stratified ~±20 %). This sets realistic expectations for the model.

### 2.4 Oil–water through a venturi
- Multiphase flow metering literature: Falcone, Hewitt & Alimonti, *Multiphase Flow Metering* (2009); Pal (1993) emulsions through orifice/venturi meters; Meng et al. (2010) venturi + conductance sensor.
- Common and reasonably accurate practice: **homogeneous model** through the venturi (oil–water mixes readily, density ratio ~1). Errors grow when phases are poorly mixed (stratified upstream).
- Stratified flow through a contraction is non-trivial: the level-gradient (hydrostatic) term `ρ g dh/dx` becomes important, analogous to shallow-water flow; critical-flow / hydraulic-jump-like behavior is possible. This is a research-level extension, not a standard closure.

### 2.5 High-fidelity multiphase CFD (for reference / optional validation)
- Interface capturing: VOF (Hirt & Nichols 1981), level set (Osher & Sethian), phase-field. OpenFOAM `interFoam`, `multiphaseEulerFoam`.
- Eulerian–Eulerian two-fluid CFD with population balance for droplets.
- Cost: 3D, turbulent, long pipe → hours to days per case. Useful to study the venturi region, not for parametric sweeps.

### 2.6 1D transient multiphase (industry)
- Commercial: OLGA, LedaFlow; open: TRACE-style. Based on the 1D transient two-fluid model, which is conditionally **ill-posed** (complex characteristics) beyond the Kelvin–Helmholtz limit — needs care (Barnea & Taitel 1994; Issa & Kempf 2003 slug capturing). Only needed if Stage 5 is required for two-phase.


## 3. Operating conditions (decided)

| Quantity | Value |
|---|---|
| Pipe diameter D | 0.10 m (A = 7.85×10⁻³ m²), horizontal, gravity in −y |
| Water | ρ = 1000 kg/m³, μ = 1.0 mPa·s |
| Light crude (≈ 40° API) | ρ ≈ 820 kg/m³, μ ≈ 3 mPa·s (configurable 2–5) |
| Oil–water interfacial tension | σ ≈ 0.025 N/m (configurable) |
| Venturi | ISO 5167-4 classical: β = d/D = 0.5 (configurable 0.3–0.75), convergent 21°, throat length d, divergent 15° included angle |
| Mixture velocity (proposed) | 0.1–3 m/s ≈ 3–85 m³/h ≈ 430–13 000 bbl/d |
| Water cut | 0–100 % |

Reynolds numbers, `Re = ρVD/μ`:

| V (m/s) | Water | Light crude |
|---|---|---|
| 0.1 | 1.0×10⁴ | 2.7×10³ |
| 0.5 | 5.0×10⁴ | 1.4×10⁴ |
| 1.0 | 1.0×10⁵ | 2.7×10⁴ |
| 3.0 | 3.0×10⁵ | 8.2×10⁴ |

Flow is turbulent over essentially the entire operating range; in the venturi throat (β = 0.5) velocity and Re are 4× higher. Expected patterns for light oil in a 10 cm pipe (from Trallero 1997 / Elseth 2001 type maps): stratified/wavy at V ≲ 0.5 m/s, stratified with interfacial mixing up to ~1–1.5 m/s, dispersed above.

## 4. Laminar–turbulent transition and turbulence treatment

### 4.1 Physics
- Hagen–Poiseuille flow is **linearly stable at all Re**; transition is subcritical and needs finite-amplitude disturbances. Turbulence becomes sustained above Re ≈ 2040 (Avila et al. 2011, *Science*: crossover of puff decay and splitting); intermittent puffs/slugs in ~2000–3000; fully turbulent above ~4000 in normal (disturbed) conditions. With extreme care, laminar flow has been kept to Re ~10⁵ (Pfenninger 1961) — irrelevant in production.
- For this project Re < 2300 only for V < 0.023 m/s (water) or V < 0.084 m/s (light crude), i.e. < 0.7 m³/h. **Transition is not a design regime here**; it appears only in the lowest-rate cases.
- **Venturi convergent — relaminarization risk.** Strong acceleration can laminarize the wall layer when `K = (ν/U²)(dU/dx) ≳ 3×10⁻⁶` (Launder 1964; Narasimha & Sreenivasan 1979). Estimate for β = 0.5: at 1 m/s, K ≈ 3×10⁻⁶ for water and ≈ 1×10⁻⁵ for light crude; K scales as 1/U. So partial relaminarization in the convergent is plausible at low rates, especially oil-continuous flow — a real effect the model should be able to show.
- **Venturi divergent**: adverse pressure gradient re-triggers turbulence and can cause separation (15° included angle is near the limit).
- **Two-phase**: each layer has its own Re; at low rates a stratified oil layer can be transitional while water is turbulent. The interface damps turbulence normal to it (Egorov 2004 damping; Fabre & Masbernat).

### 4.2 How the solver handles it
Transition is not modeled by a switch or correlation. The 3D solver resolves the unsteady flow and transition emerges from the equations, provided three things:
1. **Disturbances are supplied.** A laminar initial field at Re = 10⁵ stays laminar numerically for a long time. Inflow comes from a periodic precursor pipe that is seeded with perturbations and sustains its own turbulence (or decays to laminar if Re is too low, which is the physically correct result).
2. **The subgrid (LES) model vanishes in laminar flow and at walls.** Use **WALE** (Nicoud & Ducros 1999) or **Vreman** (2004); both give zero eddy viscosity in pure shear/laminar flow. Plain Smagorinsky is rejected: it adds viscosity to laminar shear and suppresses transition. Dynamic Smagorinsky (Germano 1991) is an option but is noisier on immersed boundaries.
3. **The wall treatment does not assume turbulence.** An equilibrium log-law wall model is wrong in laminar or relaminarizing flow. Use an ODE-based equilibrium wall model with a damped mixing length (e.g. Kawai & Larsson 2012 / Bose & Park 2018 review), which reduces to the laminar profile when the local Re is low, and switch it off where the first cell is wall-resolved (y⁺ < ~5).

RANS transition models (γ–Re_θ, Langtry & Menter 2009) are not used: they give only mean fields (not the 3D unsteady visualization wanted) and are calibrated for external boundary layers.

### 4.3 Resolution budget (Re_τ ≈ 0.09 Re^0.88)
| Case | Re_τ | Uniform-grid approach | Cells (15 D domain) |
|---|---|---|---|
| Verification, Re = 5300 | ~180 | Wall-resolved LES / quasi-DNS, D/200 | ~1–2×10⁸ (5 D periodic: ~2×10⁷) |
| Oil, V = 1 m/s, Re = 2.7×10⁴ | ~720 | Wall-modeled LES, D/100 | ~1.5×10⁷ |
| Water, V = 1 m/s, Re = 10⁵ | ~2300 | Wall-modeled LES, D/100 | ~1.5×10⁷ |
| Wall-resolved at Re = 10⁵ | ~2300 | D/2300 uniform | ~10¹¹ — infeasible |

Wall modeling is therefore required for production cases on a Cartesian grid. Time step from CFL in the throat: ~10⁻⁴ s; ~10⁵ steps for statistics → roughly 0.5–4 h per case on one data-centre GPU, ~2 days on 4 CPU cores.

## 5. Numerical method

**Governing equations** (one-fluid formulation for both stages):
```
∇·u = 0
ρ(φ)[∂u/∂t + (u·∇)u] = −∇p + ∇·[(μ(φ)+μ_sgs)(∇u+∇uᵀ)] + ρ(φ)g + σκ n δ_s + f_drive
∂φ/∂t + ∇·(φu) = 0          (φ = water volume fraction; φ ≡ 1 for single phase)
```
- **Grid**: uniform Cartesian, staggered (MAC).
- **Pipe and venturi walls**: immersed boundary via a signed-distance function of the ISO venturi profile (Boundary Data Immersion Method, Weymouth & Yue 2011). No mesh generation; geometry is a function `r_wall(x)`.
- **Time integration**: explicit 2nd-order predictor–corrector projection (Chorin/Kim–Moin type); pressure Poisson solved by geometric multigrid.
- **Convection**: low-dissipation central/QUICK-limited scheme + WALE SGS.
- **Two-phase**: conservative geometric VOF (PLIC interface reconstruction, mass conserving to machine precision), consistent mass–momentum transport (Rudman 1997 type, required for stability), surface tension by CSF (Brackbill 1992) with height-function curvature. Density ratio 0.82 and viscosity ratio ~3 are easy compared with gas–liquid.
- **Driving**: constant flow rate enforced by a uniform body force adjusted each step (precursor) — so volumetric production is a direct input; water cut set by inflow φ profile.
- **Inflow/outflow**: turbulent inflow from a periodic precursor (straight pipe, 5–6 D), convective outflow.

### 5.1 Implementation base (Julia)
Recommended: build on **WaterLily.jl** (Weymouth et al., Comput. Phys. Commun. 2025) — incompressible, Cartesian, BDIM immersed boundaries, GPU/CPU via KernelAbstractions.jl, differentiable, <1000 lines — and the VOF extension **InterfaceAdvection.jl** (Huang & Weymouth), which provides conservative VOF on WaterLily.
What this project adds: pipe/venturi geometry, flow-rate forcing, precursor inflow, gravity + two-phase inflow, WALE SGS, wall model on the immersed boundary, diagnostics (p(x), holdup, ΔP_venturi, C_d, friction factor), output.
Alternative: write the solver from scratch in Julia on KernelAbstractions.jl. More control, ~3–5× more work, same numerics. Other Julia options considered: Oceananigans.jl (ocean-oriented, less suited to pipe walls), Trixi.jl (compressible), Gridap/Ferrite (FEM, slower for LES).

### 5.2 Known limitation: dispersed flow
In dispersed patterns droplets are ~0.1–2 mm (Hinze 1955 scaling), i.e. at or below the grid spacing (1 mm at D/100). VOF cannot resolve them; the interface would be numerically smeared rather than physically dispersed. VOF-LES is valid for stratified, wavy and interfacial-mixing regimes (V ≲ 1–1.5 m/s). For fully dispersed flow the options are a single-fluid mixture model with emulsion viscosity (cheap, no slip) or Euler–Euler with population balance (heavy). Proposed: mixture-model option in Stage 4b.

## 6. Code layout
```
Project.toml
src/TwoPhaseFlow.jl
  fluids.jl        # Fluid, presets (water, light crude)
  geometry.jl      # pipe SDF, ISO 5167-4 venturi SDF
  setup.jl         # case builder: D, Q_total, water cut, grid, domain
  forcing.jl       # flow-rate forcing, precursor/inflow coupling
  sgs.jl           # WALE / Vreman eddy viscosity
  wallmodel.jl     # ODE equilibrium wall model on the immersed boundary
  twophase.jl      # VOF coupling, gravity, surface tension settings
  diagnostics.jl   # p(x), friction factor, holdup, ΔP, C_d, statistics
  io.jl            # VTK (WriteVTK.jl → ParaView), JLD2 checkpoints
  onedim.jl        # small 1D reference model (Colebrook, ISO C_d, 2-fluid stratified) for cross-checks
examples/          # straight_pipe.jl, venturi.jl, oilwater_pipe.jl, oilwater_venturi.jl
test/              # runtests.jl
docs/
```
Visualization: VTK output for ParaView (iso-surfaces of the oil–water interface, Q-criterion vortices, velocity slices); Makie.jl for quick plots and animations.

## 7. Stages and validation

| Stage | Deliverable | Validation |
|---|---|---|
| 0 | Geometry, grid, I/O, 1D reference model | Unit tests; SDF accuracy; Colebrook/ISO values |
| 1a | Laminar pipe | Hagen–Poiseuille profile and `f = 64/Re` at Re = 500–1500 (grid convergence study) |
| 1b | Transition check | Re = 1500 perturbed → decays to laminar; Re = 3000–5300 → sustained turbulence |
| 1c | Turbulent pipe | Re = 5300 vs DNS of El Khoury et al. (2013) (mean profile, rms, friction); Re = 10⁴–10⁵ with wall model vs Colebrook/McKeon (2004) friction factor |
| 2 | Single-phase venturi | ΔP and C_d vs ISO 5167-4; permanent loss 5–20 % of ΔP; relaminarization check in convergent |
| 3 | Oil–water straight pipe | Stratified laminar two-layer analytical solution (exact); holdup and dp/dx vs 1D two-fluid model and published data (Elseth 2001, Rodriguez & Oliemans 2006); pattern vs Trallero map |
| 4 | Oil–water venturi | ΔP vs homogeneous-model prediction and published oil–water venturi data; effect of upstream pattern |

## 8. Assumptions and limitations
- Isothermal, Newtonian, immiscible, no gas, no solids.
- Smooth walls initially (roughness via wall model later).
- Cartesian immersed boundary: walls are resolved to ~1 cell; near-wall accuracy relies on the wall model at production Re.
- Dispersed-regime droplets not resolved (see 5.2).
- Development and small tests run in this environment (4 CPU cores, no GPU); production runs need a GPU.

## 9. Decisions (3D)
| Topic | Decision |
|---|---|
| Hardware | Modal, single **H100 (80 GB)** |
| Software | Reuse existing packages where possible, but **verify each thoroughly** (code review, own test suite, benchmarks) before relying on it |
| Rate range | Mixture velocity 0.1–3 m/s (≈ 3–85 m³/h); 3 m/s is the practical maximum |
| Dispersed flow | **Must be handled** — VOF alone is insufficient (droplets below grid scale). Approach under discussion, see §10 |
| 1D model | Done (`onedim/`); used for inflow holdup and cross-checks |

### 9.1 Findings from checking candidate packages (2026-09-24)
- **WaterLily.jl** v1.8.0 (MIT, active, last commit 2026-09-21): incompressible, Cartesian, BDIM immersed boundaries, CUDA/AMD extensions, VTK/JLD2 output, SGS hook (`sgs!` with a Smagorinsky example; WALE must be added). Single phase only.
- **InterfaceAdvection.jl** v1.0.0-DEV (MIT, unregistered, last commit 2026-09-14): conservative VOF with consistent mass–momentum transport, surface tension, bounded viscosity interpolation, GPU. **Does not support immersed bodies yet** (README goal: "Reintroduce the boundary data immersion method"; `# TODO: support BDIM body` in source) and has an open issue with symmetry BCs under gravity. So pipe walls cannot be represented in the two-phase solver as-is. Options: (a) add BDIM to it (consistent with its momentum-form scheme — non-trivial, must be verified); (b) represent the pipe wall another way (volume penalisation); (c) write our own variable-density step on WaterLily. To be decided after a code review.

## 10. Dispersed flow (under discussion)
See the discussion in the session; summary of options:
| Option | What it resolves | Feasible on one H100? |
|---|---|---|
| Interface-resolved VOF in the pipe | every drop | No — needs Δ ≈ 50 µm, ~10¹¹ cells for 15 D |
| LES + mixture (drift-flux) model + transported drop size | concentration field, creaming layer, emulsion viscosity, pressure drop | Yes — hours per case |
| Euler–Euler two-fluid + population balance | as above plus separate phase velocities | Yes, 2–3× cost; kernels uncertain, LES with E–E less mature |
| Euler–Lagrange point particles | individual drops (dilute only) | No — invalid at 10–50 % volume fraction and d ≈ Δ |
| Hybrid VOF + sub-grid dispersed phase | large interfaces resolved, small drops modelled | Yes, research-grade |
| Interface-resolved "microscope" box (~2–3 cm, periodic, forced turbulence) | drop break-up/coalescence physics, drop sizes, emulsion viscosity calibration | Yes — ~1 day per case |

## 11. Open questions (3D)
1. Dispersed-flow approach (§10).
2. WaterLily + InterfaceAdvection vs own variable-density step, after code review (§9.1).
