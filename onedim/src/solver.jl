# Steady 1D marching solver.
#
# Mixture momentum over a control volume of the pipe (horizontal, steady, incompressible):
#     d(A J)/dx = −A dp/dx − τ S − (local form losses)
# with J = ρ_w α_w u_w² + ρ_o α_o u_o² the momentum flux. Hence
#     dp/dx = −(1/A) d(A J)/dx + (dp/dx)_friction − (dp/dx)_local.
# Discretisation between nodes i and i+1:
#     p_{i+1} = p_i − (A_{i+1}J_{i+1} − A_i J_i)/Ā + Δx·[G_i + G_{i+1}]/2 − Δx·L
# where G = (dp/dx)_friction and Ā = 2A_iA_{i+1}/(A_i + A_{i+1}) (harmonic mean). With this
# Ā the acceleration term reproduces Bernoulli exactly for single-phase and homogeneous flow,
# so a frictionless venturi recovers the ideal pressures to round-off.
# Local losses in diffusers are spread uniformly over the section length, and wall friction
# is not added where a loss coefficient applies (the coefficient is the section's total loss).

"""
    Solution

Result of `solve`. Node arrays: `x`, `D`, `p` [Pa, relative to inlet], `Hw`, `uw`, `uo`,
`dpdx` (frictional gradient), `pattern`. Also the inputs `pipe`, `sys`, `rates`.
"""
struct Solution
    x::Vector{Float64}
    D::Vector{Float64}
    p::Vector{Float64}
    Hw::Vector{Float64}
    uw::Vector{Float64}
    uo::Vector{Float64}
    dpdx::Vector{Float64}
    pattern::Vector{FlowPattern}
    nroots::Vector{Int}
    pipe::PipeGeometry
    sys::OilWater
    rates::Rates
end

"Axial grid: segment boundaries, taps, and uniform subdivision with spacing ≤ dx."
function build_grid(pipe::PipeGeometry, dx)
    xs = Float64[]
    for (i, s) in enumerate(pipe.segments)
        x0, L = pipe.xstart[i], seglength(s)
        n = max(4, ceil(Int, L / dx))
        append!(xs, x0 .+ L .* (0:n-1) ./ n)
    end
    push!(xs, pipe_length(pipe))
    append!(xs, values(pipe.taps))
    sort!(xs)
    # merge points closer than a tiny tolerance
    out = [xs[1]]
    for x in xs[2:end]
        x - out[end] > 1e-12 * max(1.0, abs(x)) && push!(out, x)
    end
    return out
end

"""
    solve(pipe, sys::OilWater, rates::Rates; model=Mechanistic(), friction=Churchill(),
          losses=CraneDiffuser(), dx=pipe_D/20) -> Solution
    solve(pipe, fluid::Fluid, Q; kwargs...) -> Solution          (single phase)

March the steady 1D momentum balance from inlet (p = 0) to outlet.
"""
function solve(pipe::PipeGeometry, sys::OilWater, rates::Rates; model::TwoPhaseModel=Mechanistic(),
               friction::FrictionModel=Churchill(), losses::LocalLosses=CraneDiffuser(),
               dx::Real=minimum(segmax_diameter, pipe.segments) / 20)
    x = build_grid(pipe, dx)
    n = length(x)
    ε = pipe.roughness
    Qo, Qw = rates.Qo, rates.Qw
    Q = Qo + Qw

    D = [diameter(pipe, xi) for xi in x]
    A = @. π * D^2 / 4
    p = zeros(n)
    st = Vector{LocalState}(undef, n)

    # interval k = [x_k, x_{k+1}] lies in segment seg(k) (evaluated at the midpoint)
    segof(k) = segment_index(pipe, (x[k] + x[k+1]) / 2)
    state(k, sidx) = local_state(model, sys, Qo, Qw, D[k]; ε, friction,
                                 segment_tag=pipe.segments[sidx].tag)
    prev_seg = 0
    left = state(1, segof(1))
    st[1] = left
    for k in 1:n-1
        si = segof(k)
        seg = pipe.segments[si]
        if si != prev_seg && k > 1
            # model may change across a segment boundary: re-evaluate the node state and
            # apply the resulting momentum-flux jump as a pressure jump
            new_left = state(k, si)
            p[k] -= new_left.J - left.J
            left = new_left
        end
        right = state(k + 1, si)
        Δx = x[k+1] - x[k]
        Ā = 2A[k] * A[k+1] / (A[k] + A[k+1])
        accel = (A[k+1] * right.J - A[k] * left.J) / Ā
        K = loss_coefficient(losses, seg)
        if K === nothing
            fric = Δx * (left.dpdx + right.dpdx) / 2
        else
            Asmall = π * min(segdiameter(seg, 0.0), segdiameter(seg, seg.L))^2 / 4
            ρm = mixture_density(sys, Qw / Q)
            fric = -Δx * K * ρm * (Q / Asmall)^2 / 2 / seg.L
        end
        p[k+1] = p[k] - accel + fric
        st[k+1] = right
        left = right
        prev_seg = si
    end
    # nodes at a segment boundary report the downstream state already (right of interval k-1
    # equals left of interval k up to the model switch); keep what the march produced
    return Solution(x, D, p, [s.Hw for s in st], [s.uw for s in st], [s.uo for s in st],
                    [s.dpdx for s in st], [s.pattern for s in st], [s.nroots for s in st],
                    pipe, sys, rates)
end

segmax_diameter(s::Cylinder) = s.D
segmax_diameter(s::Cone) = max(s.D1, s.D2)

function solve(pipe::PipeGeometry, fluid::Fluid, Q::Real; kwargs...)
    sys = OilWater(oil=fluid, water=fluid, σ=0.0)
    return solve(pipe, sys, Rates(0.0, Q); kwargs...)
end

"Pressure [Pa] at axial position `x` (linear interpolation) or at a named tap."
function pressure_at(sol::Solution, x::Real)
    i = searchsortedlast(sol.x, x)
    i >= length(sol.x) && return sol.p[end]
    i < 1 && return sol.p[1]
    t = (x - sol.x[i]) / (sol.x[i+1] - sol.x[i])
    return (1 - t) * sol.p[i] + t * sol.p[i+1]
end
pressure_at(sol::Solution, name::Symbol) = pressure_at(sol, tap(sol.pipe, name))

"Pressure drop p(a) − p(b) between two positions or taps (default: inlet → outlet)."
pressure_drop(sol::Solution, a=sol.x[1], b=sol.x[end]) = pressure_at(sol, a) - pressure_at(sol, b)

function Base.show(io::IO, ::MIME"text/plain", s::Solution)
    Q = total_rate(s.rates)
    @printf(io, "PipeFlow1D.Solution: L = %.3f m, %d nodes\n", s.x[end], length(s.x))
    @printf(io, "  Q = %.4g m³/h (%.0f bbl/d), water cut = %.3f\n", to_m3h(Q), to_bpd(Q), watercut(s.rates))
    @printf(io, "  Δp inlet→outlet = %.2f Pa (%.4f bar)\n", pressure_drop(s), pressure_drop(s) / 1e5)
    pats = unique(s.pattern)
    print(io, "  patterns: ", join(string.(pats), ", "))
    @printf(io, "\n  water holdup at inlet = %.4f", s.Hw[1])
end
