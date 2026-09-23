# Steady, fully developed stratified oil–water flow: two-fluid model.
#
# Water (denser) occupies a circular segment at the bottom. With φ the half-angle
# subtended at the pipe axis by the wetted water perimeter (0 < φ < π), R = D/2:
#     h   = R (1 − cos φ)                interface height
#     A_w = R² (φ − sin φ cos φ)         water area
#     S_w = 2 R φ,  S_o = 2 R (π − φ)    wall perimeters
#     S_i = 2 R sin φ                    interface width
#
# Momentum, per unit length (Taitel & Dukler 1976; Brauner & Moalem Maron 1989;
# Trallero 1995), with τ_i > 0 when oil moves faster than water:
#     −A_w dp/dx − τ_w S_w + τ_i S_i = 0
#     −A_o dp/dx − τ_o S_o − τ_i S_i = 0
# Eliminating dp/dx gives the holdup equation
#     F(φ) = τ_o S_o/A_o − τ_w S_w/A_w + τ_i S_i (1/A_o + 1/A_w) = 0,
# and summing gives dp/dx = −(τ_w S_w + τ_o S_o)/A.
#
# Hydraulic diameters: the faster phase "sees" the interface as a wall
# (Agrawal et al. 1973; Taitel & Dukler 1976):
#     u_o > u_w:  D_o = 4A_o/(S_o + S_i),  D_w = 4A_w/S_w
#     otherwise:  D_w = 4A_w/(S_w + S_i),  D_o = 4A_o/S_o
# This makes F discontinuous where u_o = u_w. A sign change across that point is a
# "no-slip lock" and is returned as a solution with `locked = true`.

abstract type InterfacialFriction end

"""
    FasterPhase()

Interfacial shear computed with the friction factor and density of the faster phase,
τ_i = f_F,fast ρ_fast (u_o − u_w)|u_o − u_w|/2. Default.
"""
struct FasterPhase <: InterfacialFriction end

"""
    ConstantInterfacial(fi)

Constant Fanning interfacial friction factor `fi`, density of the faster phase.
"""
struct ConstantInterfacial <: InterfacialFriction
    fi::Float64
end

"""
    segment_geometry(φ, D) -> (A_w, A_o, S_w, S_o, S_i, h)

Geometry of the water segment for wetted half-angle `φ`.
"""
function segment_geometry(φ, D)
    R = D / 2
    Aw = R^2 * (φ - sin(φ) * cos(φ))
    Ao = π * R^2 - Aw
    return (Aw=Aw, Ao=Ao, Sw=2R * φ, So=2R * (π - φ), Si=2R * sin(φ), h=R * (1 - cos(φ)))
end

"Evaluate stresses and the holdup residual at half-angle φ."
function stratified_eval(φ, sys::OilWater, Qo, Qw, D, ε, fm::FrictionModel, im::InterfacialFriction)
    g = segment_geometry(φ, D)
    uw, uo = Qw / g.Aw, Qo / g.Ao
    oil_faster = uo > uw
    Dw = oil_faster ? 4g.Aw / g.Sw : 4g.Aw / (g.Sw + g.Si)
    Do = oil_faster ? 4g.Ao / (g.So + g.Si) : 4g.Ao / g.So
    τw = wall_shear(fm, sys.water.ρ, sys.water.μ, uw, Dw, ε)
    τo = wall_shear(fm, sys.oil.ρ, sys.oil.μ, uo, Do, ε)
    τi = interfacial_shear(im, fm, sys, uo, uw, Do, Dw, ε, oil_faster)
    F = τo * g.So / g.Ao - τw * g.Sw / g.Aw + τi * g.Si * (1 / g.Ao + 1 / g.Aw)
    dpdx = -(τw * g.Sw + τo * g.So) / (g.Aw + g.Ao)
    return (F=F, dpdx=dpdx, uw=uw, uo=uo, τw=τw, τo=τo, τi=τi, geom=g, oil_faster=oil_faster)
end

function interfacial_shear(::FasterPhase, fm, sys, uo, uw, Do, Dw, ε, oil_faster)
    Δu = uo - uw
    Δu == 0 && return 0.0
    if oil_faster
        f = fanning(fm, sys.oil.ρ * abs(uo) * Do / sys.oil.μ, ε / Do)
        ρ = sys.oil.ρ
    else
        f = fanning(fm, sys.water.ρ * abs(uw) * Dw / sys.water.μ, ε / Dw)
        ρ = sys.water.ρ
    end
    return f * ρ * Δu * abs(Δu) / 2
end

function interfacial_shear(m::ConstantInterfacial, fm, sys, uo, uw, Do, Dw, ε, oil_faster)
    Δu = uo - uw
    ρ = oil_faster ? sys.oil.ρ : sys.water.ρ
    return m.fi * ρ * Δu * abs(Δu) / 2
end

"""
    stratified_solutions(sys, Qo, Qw, D; ε=0, friction=Churchill(), interfacial=FasterPhase(),
                         nscan=2000) -> Vector of NamedTuples

All solutions of the stratified holdup equation, ordered by increasing water level.
Each has fields `φ, Hw, h, dpdx, uw, uo, locked`. Multiple solutions can exist
(Ullmann et al. 2003); the caller chooses. Requires Qo > 0 and Qw > 0.
"""
function stratified_solutions(sys::OilWater, Qo, Qw, D; ε=0.0, friction::FrictionModel=Churchill(),
                              interfacial::InterfacialFriction=FasterPhase(), nscan::Int=2000)
    (Qo > 0 && Qw > 0) || throw(ArgumentError("stratified model needs both phases flowing"))
    F(φ) = stratified_eval(φ, sys, Qo, Qw, D, ε, friction, interfacial).F
    # Scan uniformly in φ; F → −∞ as φ → 0 and +∞ as φ → π, so a root always exists.
    δ = 1e-6
    φs = range(δ, π - δ; length=nscan)
    Fs = F.(φs)
    sols = NamedTuple[]
    for k in 1:nscan-1
        (isfinite(Fs[k]) && isfinite(Fs[k+1])) || continue
        sign(Fs[k]) == sign(Fs[k+1]) && Fs[k] != 0 && continue
        φ = Fs[k] == 0 ? φs[k] : brent(F, φs[k], φs[k+1]; xtol=1e-13)
        e = stratified_eval(φ, sys, Qo, Qw, D, ε, friction, interfacial)
        # A root of a continuous branch has |F| ≪ the bracketing values; a sign change
        # across the u_o = u_w switch of hydraulic diameters does not.
        scale = max(abs(Fs[k]), abs(Fs[k+1]))
        locked = abs(e.F) > 1e-3 * scale
        if locked
            # Place the solution exactly at u_o = u_w (no slip) within the bracket.
            G(φ) = (g = segment_geometry(φ, D); Qo / g.Ao - Qw / g.Aw)
            φ = brent(G, φs[k], φs[k+1]; xtol=1e-13)
            e = stratified_eval(φ, sys, Qo, Qw, D, ε, friction, interfacial)
            # At the lock both hydraulic-diameter conventions apply; report the mean
            # of the two one-sided pressure gradients.
            e2 = stratified_eval(nextfloat(φ, 1000), sys, Qo, Qw, D, ε, friction, interfacial)
            e = merge(e, (dpdx=(e.dpdx + e2.dpdx) / 2,))
        end
        push!(sols, (φ=φ, Hw=e.geom.Aw / (e.geom.Aw + e.geom.Ao), h=e.geom.h, dpdx=e.dpdx,
                     uw=e.uw, uo=e.uo, locked=locked))
    end
    isempty(sols) && error("no stratified solution found (Qo=$Qo, Qw=$Qw, D=$D)")
    return sols
end
