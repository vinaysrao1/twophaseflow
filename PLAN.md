# twophaseflow — Plan (draft for review)

Status: **draft, not yet approved**. No simulation code is written until this plan is finalized.

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

## 3. Recommended approach

**Primary engine: 1D steady mechanistic model** marching along `x` over an arbitrary area profile `A(x)`.
Rationale: this is what the research literature and industry use for design with variable water cut and rates; it is fast (milliseconds per case, so full sweeps over WC × Q are cheap); its accuracy is bounded by the closures, not the numerics; and CFD for turbulent oil–water is expensive with no guarantee of better accuracy without calibration.

**Secondary (optional) verification solver:** 2D axisymmetric laminar finite-volume Navier–Stokes (single phase) to verify the 1D venturi treatment against a resolved velocity field at low Re, and to show recirculation in the diffuser. For turbulent/two-phase high-fidelity cases, generate OpenFOAM case files rather than writing our own turbulent/VOF solver.

### 3.1 Architecture (Python)

```
twophaseflow/
  fluids.py          # Fluid dataclass (ρ, μ); presets: water, light/medium/heavy crude
  geometry.py        # Pipe, AreaProfile; Venturi (ISO 5167-4 geometry builder)
  friction.py        # laminar, Colebrook, Haaland, Churchill
  losses.py          # contraction/diffuser loss coefficients (Idelchik/Crane)
  single_phase.py    # 1D marching solver: dp/dx = −f ρV²/(2D) − ρ V dV/dx (+ local losses)
  mixture.py         # emulsion viscosity laws, phase inversion criteria
  patterns.py        # flow-pattern prediction (Trallero / Brauner criteria)
  stratified.py      # two-fluid stratified model, geometry of segments, multiple-root detection
  dispersed.py       # homogeneous model
  two_phase.py       # pattern-dispatching solver, same marching interface as single phase
  transient.py       # (Stage 5) rigid-column transient
  cfd/axisym.py      # (optional) 2D axisymmetric laminar FV solver
  cli.py / plotting.py
tests/               # pytest: analytic + literature validation
examples/            # scripts/notebooks: pressure profile, WC × Q sweeps, flow-pattern maps
```
Dependencies: numpy, scipy, matplotlib, pytest. Closures are pluggable (strategy objects) so correlations can be swapped and compared.

### 3.2 Staged deliverables and validation

| Stage | Deliverable | Validation |
|---|---|---|
| 1 | Single-phase straight pipe | Hagen–Poiseuille (exact); Colebrook vs Moody/Haaland; mass conservation |
| 2 | Venturi (1D) | Ideal Bernoulli limit (f=0, no losses) exact; ISO 5167-4 `C_d` and permanent-loss ranges; optional 2D axisymmetric laminar CFD comparison |
| 3 | Oil–water straight pipe: pattern map, holdup, dp/dx vs WC and Q | Limits WC→0 and WC→1 reduce to single phase; published data (Trallero 1995, Elseth 2001, Angeli & Hewitt 2000, Rodriguez & Oliemans 2006 — digitized from papers where data are not tabulated); report error bands |
| 4 | Oil–water venturi | Homogeneous model vs published venturi oil–water data; sensitivity to mixing assumption |
| 5 | Transient (optional) | Rigid-column analytical solutions (step change in Δp → exponential/tanh approach to steady state) |

Each stage: tests pass, example plots committed, short results note in `docs/`.

## 4. Assumptions and limitations (stated up front)
- Isothermal, Newtonian fluids (waxy/heavy crude can be non-Newtonian — excluded unless requested).
- No gas, no solids, no mass transfer; immiscible liquids.
- Fully developed flow in straight sections; venturi handled quasi-1D with loss coefficients.
- Phase inversion and emulsion viscosity have large empirical uncertainty; the model exposes these as selectable closures rather than one "true" answer.
- Stratified flow through a venturi uses the homogeneous assumption by default; a stratified-with-level-gradient variant is flagged as experimental.

## 5. Open questions (need your answer to finalize)
1. **Purpose**: engineering prediction (pressure drop, holdup, venturi ΔP for metering) — suits the 1D model — or detailed flow-field visualization — needs CFD?
2. **Steady vs transient**: is steady state per operating point enough, or do you need time-varying production rates?
3. **Operating ranges**: pipe diameter, rates (e.g. bbl/d or m³/h), crude API/viscosity. Heavy crude (> ~500 cP) brings core-annular flow and non-Newtonian effects into scope.
4. **2D verification solver**: include the in-house 2D axisymmetric laminar solver, OpenFOAM case generation, or neither?
5. **Language/tooling**: Python (numpy/scipy) OK? Any preference for CLI vs notebooks vs a simple GUI/dashboard?
6. **Venturi geometry**: ISO 5167-4 classical venturi, or a custom profile (you specify β, angles)?
