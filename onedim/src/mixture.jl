# Dispersion (emulsion) properties and phase inversion.

abstract type EmulsionViscosity end

"""
    Brinkman()

Brinkman (1952), J. Chem. Phys. 20:571: μ = μ_c (1 − φ)^(−2.5), φ = dispersed fraction.
Widely used for oil–water dispersions in pipes (e.g. Brauner 2003). Default.
"""
struct Brinkman <: EmulsionViscosity end

"""
    KriegerDougherty(; φmax=0.74, intrinsic=2.5)

Krieger & Dougherty (1959), Trans. Soc. Rheol. 3:137:
μ = μ_c (1 − φ/φmax)^(−intrinsic·φmax). φ is capped just below φmax.
"""
Base.@kwdef struct KriegerDougherty <: EmulsionViscosity
    φmax::Float64 = 0.74
    intrinsic::Float64 = 2.5
end

"""
    Taylor()

Taylor (1932) dilute-emulsion limit: μ = μ_c [1 + 2.5 φ (μ_d + 0.4 μ_c)/(μ_d + μ_c)].
Valid only for φ ≲ 0.05.
"""
struct Taylor <: EmulsionViscosity end

emulsion_viscosity(::Brinkman, μc, μd, φ) = μc * (1 - φ)^(-2.5)
emulsion_viscosity(m::KriegerDougherty, μc, μd, φ) =
    μc * (1 - min(φ, 0.999m.φmax) / m.φmax)^(-m.intrinsic * m.φmax)
emulsion_viscosity(::Taylor, μc, μd, φ) = μc * (1 + 2.5φ * (μd + 0.4μc) / (μd + μc))

abstract type InversionModel end

"""
    Arirachakaran()

Arirachakaran et al. (1989), SPE 18836: water cut at phase inversion
ε_w,I = 0.5 − 0.1108 log10(μ_o / 1 mPa·s), clamped to [0.05, 0.95].
Fitted for oil viscosities 4.7–115 mPa·s; extrapolated for lighter oils.
"""
struct Arirachakaran <: InversionModel end

"FixedInversion(wc): user-specified inversion water cut."
struct FixedInversion <: InversionModel
    wc::Float64
end

inversion_watercut(::Arirachakaran, sys::OilWater) =
    clamp(0.5 - 0.1108log10(sys.oil.μ / 1e-3), 0.05, 0.95)
inversion_watercut(m::FixedInversion, sys::OilWater) = m.wc

"""
    continuous_phase(inv, sys, wc) -> :water or :oil

Continuous phase of a dispersion at water cut `wc`: water above the inversion point.
"""
continuous_phase(inv::InversionModel, sys::OilWater, wc) =
    wc >= inversion_watercut(inv, sys) ? :water : :oil

"No-slip mixture density."
mixture_density(sys::OilWater, wc) = wc * sys.water.ρ + (1 - wc) * sys.oil.ρ

"""
    mixture_viscosity(visc, inv, sys, wc) -> μ

Effective viscosity of a homogeneous dispersion at water cut `wc`.
"""
function mixture_viscosity(visc::EmulsionViscosity, inv::InversionModel, sys::OilWater, wc)
    wc == 1 && return sys.water.μ
    wc == 0 && return sys.oil.μ
    if continuous_phase(inv, sys, wc) === :water
        return emulsion_viscosity(visc, sys.water.μ, sys.oil.μ, 1 - wc)
    else
        return emulsion_viscosity(visc, sys.oil.μ, sys.water.μ, wc)
    end
end
