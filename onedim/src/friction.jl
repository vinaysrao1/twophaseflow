# Darcy friction factor correlations, f_D = 4 f_Fanning.
# Δp = f_D (L/D) ρ V²/2 ;  τ_wall = f_F ρ V²/2.

abstract type FrictionModel end

"""
    Churchill()

Churchill (1977), Chem. Eng. 84(24):91. A single expression covering laminar,
transitional and turbulent (smooth and rough) flow. Continuous in Re, which keeps
the root-finding in the two-fluid model well behaved. Default.
"""
struct Churchill <: FrictionModel end

"""
    Colebrook(; Re_laminar=2300)

Colebrook–White (1939) for Re ≥ `Re_laminar`, Hagen–Poiseuille 64/Re below.
Discontinuous at `Re_laminar`.
"""
Base.@kwdef struct Colebrook <: FrictionModel
    Re_laminar::Float64 = 2300.0
end

"""
    Haaland(; Re_laminar=2300)

Haaland (1983) explicit approximation to Colebrook (±2 %), laminar below `Re_laminar`.
"""
Base.@kwdef struct Haaland <: FrictionModel
    Re_laminar::Float64 = 2300.0
end

"Hagen–Poiseuille 64/Re at any Re (verification only)."
struct Laminar <: FrictionModel end

"Zero wall friction (verification only: recovers Bernoulli)."
struct NoFriction <: FrictionModel end

"""
    darcy(model, Re, rr) -> f_D

Darcy friction factor at Reynolds number `Re` and relative roughness `rr = ε/D`.
"""
function darcy end

function darcy(::Churchill, Re, rr)
    Re > 0 || return Inf
    A = (2.457 * log(1 / ((7 / Re)^0.9 + 0.27rr)))^16
    B = (37530 / Re)^16
    return 8 * ((8 / Re)^12 + (A + B)^(-1.5))^(1 / 12)
end

function darcy(m::Colebrook, Re, rr)
    Re > 0 || return Inf
    Re < m.Re_laminar && return 64 / Re
    # Iterate on x = 1/√f starting from Haaland.
    x = 1 / sqrt(darcy(Haaland(0.0), Re, rr))
    for _ in 1:50
        xn = -2log10(rr / 3.7 + 2.51x / Re)
        abs(xn - x) < 1e-13 * abs(x) && return 1 / xn^2
        x = xn
    end
    return 1 / x^2
end

function darcy(m::Haaland, Re, rr)
    Re > 0 || return Inf
    Re < m.Re_laminar && return 64 / Re
    return (-1.8log10((rr / 3.7)^1.11 + 6.9 / Re))^-2
end

darcy(::Laminar, Re, rr) = Re > 0 ? 64 / Re : Inf
darcy(::NoFriction, Re, rr) = 0.0

"Fanning friction factor f_D/4."
fanning(m::FrictionModel, Re, rr) = darcy(m, Re, rr) / 4

"""
    wall_shear(model, ρ, μ, u, Dh, ε) -> τ

Wall shear stress [Pa] for mean velocity `u` in a duct of hydraulic diameter `Dh`,
signed with `u`. Written so that u → 0 gives τ → 0 even though f → ∞.
"""
function wall_shear(m::FrictionModel, ρ, μ, u, Dh, ε)
    u == 0 && return 0.0
    Re = ρ * abs(u) * Dh / μ
    return fanning(m, Re, ε / Dh) * ρ * u * abs(u) / 2
end
