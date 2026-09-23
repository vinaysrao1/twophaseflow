# Local (form) losses in conical sections.
#
# A 1D model with wall friction alone under-predicts the loss in a diffuser, where the
# boundary layer thickens under the adverse pressure gradient and may separate. The
# classical treatment is a loss coefficient K referenced to the small-end velocity:
#     Δp_loss = K ρ V_small² / 2.
# Crane TP-410 (from Gibson's data) for a gradual enlargement, included angle θ, β = d/D:
#     θ ≤ 45°:  K = 2.6 sin(θ/2) (1 − β²)²
#     θ > 45°:  K = (1 − β²)²                    (Borda–Carnot, sudden expansion)
# K is the total loss of the section, so wall friction is not added on top of it.
#
# Contractions of venturi quality lose little beyond wall friction (ISO 5167-4 gives
# C ≈ 0.985–0.995); they are treated with wall friction only.

abstract type LocalLosses end

"""
    CraneDiffuser(; factor=1.0)

Crane TP-410 loss coefficient in diverging cones (default), multiplied by `factor`.
For a 15° venturi diffuser with β = 0.5 it gives a permanent loss of ≈ 22 % of the
venturi Δp, at the upper end of the 5–20 % range quoted for classical venturis: Crane's
data are for fully developed inlet flow, whereas the short venturi throat delivers a
thinner boundary layer. Use `factor < 1` to calibrate against a specific meter.
"""
Base.@kwdef struct CraneDiffuser <: LocalLosses
    factor::Float64 = 1.0
end

"No local losses; wall friction only everywhere (verification)."
struct NoLocalLosses <: LocalLosses end

"""
    loss_coefficient(model, segment) -> K or nothing

Loss coefficient (small-end velocity) of `segment`, or `nothing` when the segment is
treated with wall friction only.
"""
loss_coefficient(::NoLocalLosses, s::Segment) = nothing
loss_coefficient(::CraneDiffuser, s::Cylinder) = nothing
function loss_coefficient(m::CraneDiffuser, s::Cone)
    s.D2 > s.D1 || return nothing
    θ = included_angle(s)
    b2 = (s.D1 / s.D2)^2
    return m.factor * (θ <= deg2rad(45) ? 2.6sin(θ / 2) * (1 - b2)^2 : (1 - b2)^2)
end
