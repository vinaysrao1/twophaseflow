# Local (cross-sectional) closure models. Each returns a LocalState at a given diameter.

abstract type TwoPhaseModel end

"""
    Homogeneous(; viscosity=Brinkman(), inversion=Arirachakaran())

No-slip dispersion: in-situ holdup equals the input water cut; mixture density
is volume-averaged; viscosity from an emulsion law with the continuous phase set by the
inversion model. Pressure gradient from the single-phase friction law at the mixture
Reynolds number.
"""
Base.@kwdef struct Homogeneous{V<:EmulsionViscosity,I<:InversionModel} <: TwoPhaseModel
    viscosity::V = Brinkman()
    inversion::I = Arirachakaran()
end

"""
    StratifiedTwoFluid(; interfacial=FasterPhase(), root=:lowest)

Two-fluid stratified model (see `stratified_solutions`). When several holdup solutions
exist, `root = :lowest` takes the lowest water level, `:highest` the highest.
"""
Base.@kwdef struct StratifiedTwoFluid{I<:InterfacialFriction} <: TwoPhaseModel
    interfacial::I = FasterPhase()
    root::Symbol = :lowest
end

"""
    Mechanistic(; homogeneous=Homogeneous(), stratified=StratifiedTwoFluid(),
                criterion=HinzeBarnea(), homogeneous_in_venturi=true)

Predicts the flow pattern at each cross-section and applies the stratified or
homogeneous model accordingly. With `homogeneous_in_venturi = true` (default), all
non-`:pipe` segments (entrance, convergent, throat, divergent) use the homogeneous model,
following multiphase venturi-metering practice.
"""
Base.@kwdef struct Mechanistic{H<:Homogeneous,S<:StratifiedTwoFluid,C<:DispersionCriterion} <: TwoPhaseModel
    homogeneous::H = Homogeneous()
    stratified::S = StratifiedTwoFluid()
    criterion::C = HinzeBarnea()
    homogeneous_in_venturi::Bool = true
end

"""
    LocalState

Cross-sectional state: pattern, in-situ water holdup `Hw`, phase velocities, frictional
pressure gradient `dpdx` [Pa/m] (negative for flow in +x), momentum flux `J` [Pa],
no-slip mixture density `ρm`, effective viscosity `μm` (NaN for stratified), and the
number of stratified solutions found `nroots`.
"""
struct LocalState
    pattern::FlowPattern
    Hw::Float64
    uw::Float64
    uo::Float64
    dpdx::Float64
    J::Float64
    ρm::Float64
    μm::Float64
    nroots::Int
    locked::Bool
end

"Momentum flux ρ_w α_w u_w² + ρ_o α_o u_o² for given holdup."
function momentum_flux(sys::OilWater, Qo, Qw, A, Hw)
    J = 0.0
    Qw > 0 && (J += sys.water.ρ * Qw^2 / (Hw * A^2))
    Qo > 0 && (J += sys.oil.ρ * Qo^2 / ((1 - Hw) * A^2))
    return J
end

function single_phase_state(fl::Fluid, pattern, Q, D, ε, fm)
    A = π * D^2 / 4
    U = Q / A
    τ = wall_shear(fm, fl.ρ, fl.μ, U, D, ε)
    Hw = pattern === SINGLE_WATER ? 1.0 : 0.0
    uw, uo = pattern === SINGLE_WATER ? (U, 0.0) : (0.0, U)
    return LocalState(pattern, Hw, uw, uo, -4τ / D, fl.ρ * U^2, fl.ρ, fl.μ, 1, false)
end

"""
    local_state(model, sys, Qo, Qw, D; ε=0, friction=Churchill(), segment_tag=:pipe) -> LocalState
"""
function local_state(m::TwoPhaseModel, sys::OilWater, Qo, Qw, D; ε=0.0,
                     friction::FrictionModel=Churchill(), segment_tag::Symbol=:pipe)
    Qw == 0 && return single_phase_state(sys.oil, SINGLE_OIL, Qo, D, ε, friction)
    Qo == 0 && return single_phase_state(sys.water, SINGLE_WATER, Qw, D, ε, friction)
    return two_phase_state(m, sys, Qo, Qw, D, ε, friction, segment_tag)
end

function two_phase_state(m::Homogeneous, sys, Qo, Qw, D, ε, fm, tag)
    wc = Qw / (Qo + Qw)
    A = π * D^2 / 4
    U = (Qo + Qw) / A
    ρm = mixture_density(sys, wc)
    μm = mixture_viscosity(m.viscosity, m.inversion, sys, wc)
    τ = wall_shear(fm, ρm, μm, U, D, ε)
    pat = continuous_phase(m.inversion, sys, wc) === :water ? DISPERSED_OW : DISPERSED_WO
    return LocalState(pat, wc, U, U, -4τ / D, ρm * U^2, ρm, μm, 1, false)
end

function two_phase_state(m::StratifiedTwoFluid, sys, Qo, Qw, D, ε, fm, tag)
    sols = stratified_solutions(sys, Qo, Qw, D; ε, friction=fm, interfacial=m.interfacial)
    s = m.root === :highest ? sols[end] : sols[1]
    A = π * D^2 / 4
    ρm = mixture_density(sys, Qw / (Qo + Qw))
    return LocalState(STRATIFIED, s.Hw, s.uw, s.uo, s.dpdx, momentum_flux(sys, Qo, Qw, A, s.Hw),
                      ρm, NaN, length(sols), s.locked)
end

function two_phase_state(m::Mechanistic, sys, Qo, Qw, D, ε, fm, tag)
    if m.homogeneous_in_venturi && tag !== :pipe
        return two_phase_state(m.homogeneous, sys, Qo, Qw, D, ε, fm, tag)
    end
    pat, _ = predict_pattern(m.criterion, m.homogeneous.inversion, sys, Qo, Qw, D; ε, friction=fm)
    sub = pat === STRATIFIED ? m.stratified : m.homogeneous
    return two_phase_state(sub, sys, Qo, Qw, D, ε, fm, tag)
end
