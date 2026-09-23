# CSV output and parameter sweeps (fully developed flow in a straight pipe).

"""
    write_csv(path, sol::Solution)
    write_csv(path, table::NamedTuple)

Write a solution profile, or a NamedTuple of equal-length vectors, as CSV.
"""
function write_csv(path::AbstractString, sol::Solution)
    tbl = (x_m=sol.x, D_m=sol.D, p_Pa=sol.p, Hw=sol.Hw, uw_m_s=sol.uw, uo_m_s=sol.uo,
           dpdx_fric_Pa_m=sol.dpdx, pattern=string.(sol.pattern), nroots=sol.nroots)
    return write_csv(path, tbl)
end

function write_csv(path::AbstractString, tbl::NamedTuple)
    cols = keys(tbl)
    n = length(first(tbl))
    all(length(c) == n for c in tbl) || throw(ArgumentError("columns differ in length"))
    open(path, "w") do io
        println(io, join(string.(cols), ","))
        for i in 1:n
            println(io, join((csvcell(tbl[c][i]) for c in cols), ","))
        end
    end
    return path
end
csvcell(v::AbstractFloat) = @sprintf("%.8g", v)
csvcell(v) = string(v)

"""
    sweep(sys, D; Um, wc, model=Mechanistic(), friction=Churchill(), ε=0.0) -> NamedTuple

Fully developed straight-pipe results over all combinations of mixture velocity `Um`
[m/s] and water cut `wc`. Columns: Um, wc, Q_m3h, dpdx_Pa_m, Hw, slip (Hw/wc), pattern,
nroots.
"""
function sweep(sys::OilWater, D; Um, wc, model::TwoPhaseModel=Mechanistic(),
               friction::FrictionModel=Churchill(), ε=0.0)
    A = π * D^2 / 4
    rows = [(u, w) for u in Um for w in wc]
    Us, Ws, Qs, G, H, S, P, N = (Float64[] for _ in 1:6)..., String[], Int[]
    for (u, w) in rows
        r = Rates(Qtotal=u * A, wc=w)
        s = local_state(model, sys, r.Qo, r.Qw, D; ε, friction)
        push!(Us, u); push!(Ws, w); push!(Qs, to_m3h(u * A))
        push!(G, s.dpdx); push!(H, s.Hw); push!(S, w > 0 ? s.Hw / w : NaN)
        push!(P, string(s.pattern)); push!(N, s.nroots)
    end
    return (Um=Us, wc=Ws, Q_m3h=Qs, dpdx_Pa_m=G, Hw=H, slip=S, pattern=P, nroots=N)
end

"""
    pattern_map(sys, D; Um, wc, criterion=HinzeBarnea(), inversion=Arirachakaran())
        -> Matrix{FlowPattern}  (rows: Um, columns: wc)
"""
function pattern_map(sys::OilWater, D; Um, wc, criterion::DispersionCriterion=HinzeBarnea(),
                     inversion::InversionModel=Arirachakaran(), friction::FrictionModel=Churchill())
    A = π * D^2 / 4
    return [predict_pattern(criterion, inversion, sys, (1 - w) * u * A, w * u * A, D; friction)[1]
            for u in Um, w in wc]
end
