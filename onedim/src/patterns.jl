# Flow-pattern prediction for horizontal oil–water flow (simplified mechanistic map).
#
# The full classification (Trallero et al. 1997) has six patterns; this model resolves
# the two limits that need different physics:
#   STRATIFIED   – gravity-separated layers (includes wavy / interfacial mixing: ST, ST&MI)
#   DISPERSED_OW – oil drops in continuous water (Do/w)
#   DISPERSED_WO – water drops in continuous oil (Dw/o)
# Intermediate patterns (Do/w & w, dual dispersions) are assigned to the nearer limit.

@enum FlowPattern SINGLE_WATER SINGLE_OIL STRATIFIED DISPERSED_OW DISPERSED_WO

const GRAVITY = 9.80665

abstract type DispersionCriterion end

"""
    HinzeBarnea(; C_H=0.725, dense_factor=0.0)

A dispersion is sustained when turbulence breaks drops smaller than the critical size
above which they either cream/settle or deform (Barnea 1986, 1987; applied to
liquid–liquid flow by Brauner 2001):

  d_max = C_H (σ/ρ_c)^0.6 ϵ^-0.4 (1 + dense_factor·φ_d)     Hinze (1955), ϵ = 2 f_F U_m³/D
  d_cσ  = [0.4 σ / (Δρ g)]^0.5                               deformation limit
  d_cb  = (3/8) ρ_c f_F U_m² / (Δρ g)                        turbulence vs buoyancy
  dispersed  ⇔  d_max < min(d_cσ, d_cb)

`dense_factor` (default 0) optionally enlarges d_max for concentrated dispersions
(coalescence); the dilute Hinze value is known to under-predict d_max near inversion.
"""
Base.@kwdef struct HinzeBarnea <: DispersionCriterion
    C_H::Float64 = 0.725
    dense_factor::Float64 = 0.0
end

"""
    predict_pattern(crit, inv, sys, Qo, Qw, D; ε=0, friction=Churchill())
        -> (pattern, info)

Predicted pattern at a cross-section of diameter `D`. `info` holds the drop sizes used
by the criterion (NaN for single-phase flow).
"""
function predict_pattern(crit::HinzeBarnea, inv::InversionModel, sys::OilWater, Qo, Qw, D;
                         ε=0.0, friction::FrictionModel=Churchill())
    nan = (d_max=NaN, d_crit_sigma=NaN, d_crit_buoyancy=NaN, continuous=:none)
    Qw == 0 && return SINGLE_OIL, nan
    Qo == 0 && return SINGLE_WATER, nan
    wc = Qw / (Qo + Qw)
    cont = continuous_phase(inv, sys, wc)
    c = cont === :water ? sys.water : sys.oil
    φd = cont === :water ? 1 - wc : wc
    Um = (Qo + Qw) / (π * D^2 / 4)
    f = fanning(friction, c.ρ * Um * D / c.μ, ε / D)
    ϵ = 2f * Um^3 / D
    Δρ = abs(sys.water.ρ - sys.oil.ρ)
    d_max = crit.C_H * (sys.σ / c.ρ)^0.6 * ϵ^-0.4 * (1 + crit.dense_factor * φd)
    d_cσ = Δρ > 0 ? sqrt(0.4sys.σ / (Δρ * GRAVITY)) : Inf
    d_cb = Δρ > 0 ? 3 / 8 * c.ρ * f * Um^2 / (Δρ * GRAVITY) : Inf
    pattern = d_max < min(d_cσ, d_cb) ? (cont === :water ? DISPERSED_OW : DISPERSED_WO) : STRATIFIED
    return pattern, (d_max=d_max, d_crit_sigma=d_cσ, d_crit_buoyancy=d_cb, continuous=cont)
end
