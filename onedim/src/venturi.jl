# Venturi metering quantities.

"""
    iso5167_venturi_C(finish=:machined) -> (C, uncertainty, Re_D_range)

Discharge coefficient of the classical venturi tube per ISO 5167-4, by convergent finish:
`:as_cast` (0.984, ±0.7 %, 2×10⁵–2×10⁶), `:machined` (0.995, ±1 %, 2×10⁵–1×10⁶),
`:welded` (rough-welded sheet iron, 0.985, ±1.5 %, 2×10⁵–2×10⁶).
Outside the Re_D range ISO gives no value (C is lower at low Re).
"""
function iso5167_venturi_C(finish::Symbol=:machined)
    finish === :as_cast && return (C=0.984, uncertainty=0.007, Re_range=(2e5, 2e6))
    finish === :machined && return (C=0.995, uncertainty=0.010, Re_range=(2e5, 1e6))
    finish === :welded && return (C=0.985, uncertainty=0.015, Re_range=(2e5, 2e6))
    throw(ArgumentError("finish must be :as_cast, :machined or :welded"))
end

"""
    venturi_flowrate(Δp, ρ, D, β; C=0.995) -> Q [m³/s]

Volumetric flow rate from the upstream–throat differential pressure (ISO 5167 equation,
incompressible). For oil–water use the homogeneous mixture density.
"""
venturi_flowrate(Δp, ρ, D, β; C=0.995) = C * (π * (β * D)^2 / 4) * sqrt(2Δp / (ρ * (1 - β^4)))

"Ideal (C = 1) upstream–throat differential pressure for flow rate `Q`."
venturi_ideal_dp(Q, ρ, D, β) = ρ / 2 * (Q / (π * (β * D)^2 / 4))^2 * (1 - β^4)

"""
    VenturiReport

`dp` upstream–throat differential [Pa]; `dp_ideal` Bernoulli value; `C` effective
discharge coefficient; `permanent_loss` unrecovered pressure loss relative to a straight
pipe [Pa] and `loss_fraction = permanent_loss/dp`; `ρm` no-slip density; `Re_D` pipe
Reynolds number of the mixture; `β`.
"""
struct VenturiReport
    dp::Float64
    dp_ideal::Float64
    C::Float64
    permanent_loss::Float64
    loss_fraction::Float64
    ρm::Float64
    Re_D::Float64
    β::Float64
end

"""
    venturi_report(sol; homogeneous=Homogeneous()) -> VenturiReport

Needs taps `:upstream`, `:throat`, `:outlet` (as created by `iso_venturi`).
"""
function venturi_report(sol::Solution; homogeneous::Homogeneous=Homogeneous())
    pipe = sol.pipe
    xu, xt, xo = tap(pipe, :upstream), tap(pipe, :throat), tap(pipe, :outlet)
    D, d = diameter(pipe, xu), diameter(pipe, xt)
    β = d / D
    Q = total_rate(sol.rates)
    wc = watercut(sol.rates)
    ρm = mixture_density(sol.sys, wc)
    μm = sol.rates.Qo == 0 ? sol.sys.water.μ : sol.rates.Qw == 0 ? sol.sys.oil.μ :
         mixture_viscosity(homogeneous.viscosity, homogeneous.inversion, sol.sys, wc)
    dp = pressure_at(sol, xu) - pressure_at(sol, xt)
    dpi = venturi_ideal_dp(Q, ρm, D, β)
    C = Q / venturi_flowrate(dp, ρm, D, β; C=1.0)
    # straight-pipe gradient at the outlet (fully developed, diameter D)
    G = sol.dpdx[end]
    perm = pressure_at(sol, xu) - pressure_at(sol, xo) + G * (xo - xu)
    ReD = ρm * (Q / (π * D^2 / 4)) * D / μm
    return VenturiReport(dp, dpi, C, perm, perm / dp, ρm, ReD, β)
end

function Base.show(io::IO, ::MIME"text/plain", r::VenturiReport)
    @printf(io, "Venturi (β = %.3f, Re_D = %.3g)\n", r.β, r.Re_D)
    @printf(io, "  Δp upstream–throat = %.1f Pa (ideal %.1f Pa)\n", r.dp, r.dp_ideal)
    @printf(io, "  effective C        = %.4f\n", r.C)
    @printf(io, "  permanent loss     = %.1f Pa (%.1f %% of Δp)", r.permanent_loss, 100r.loss_fraction)
end
