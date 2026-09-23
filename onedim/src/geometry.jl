# Axisymmetric pipe made of cylindrical and conical segments.

abstract type Segment end

"""
    Cylinder(L, D; tag=:pipe)

Straight section of length `L` and diameter `D`.
"""
struct Cylinder <: Segment
    L::Float64
    D::Float64
    tag::Symbol
end
Cylinder(L, D; tag=:pipe) = Cylinder(L, D, tag)

"""
    Cone(L, D1, D2; tag)

Conical section of length `L` going from diameter `D1` to `D2`
(contraction if `D2 < D1`, diffuser if `D2 > D1`).
"""
struct Cone <: Segment
    L::Float64
    D1::Float64
    D2::Float64
    tag::Symbol
end
Cone(L, D1, D2; tag=(D2 < D1 ? :convergent : :divergent)) = Cone(L, D1, D2, tag)

seglength(s::Segment) = s.L
segdiameter(s::Cylinder, ξ) = s.D
segdiameter(s::Cone, ξ) = s.D1 + (s.D2 - s.D1) * ξ / s.L

"Included angle [rad] of a cone."
included_angle(s::Cone) = 2atan(abs(s.D2 - s.D1) / (2s.L))

"""
    PipeGeometry(segments; roughness=0.0, taps=Dict{Symbol,Float64}())

Pipe assembled from `segments` in flow order (x = 0 at inlet). `roughness` is the
absolute wall roughness ε [m]. `taps` maps names to axial positions [m] at which the
solver reports pressure (e.g. venturi pressure tappings).
"""
struct PipeGeometry
    segments::Vector{Segment}
    xstart::Vector{Float64}
    roughness::Float64
    taps::Dict{Symbol,Float64}
end

function PipeGeometry(segments::AbstractVector{<:Segment}; roughness=0.0, taps=Dict{Symbol,Float64}())
    isempty(segments) && throw(ArgumentError("pipe needs at least one segment"))
    segs = Segment[segments...]
    xs = cumsum([0.0; [seglength(s) for s in segs[1:end-1]]])
    p = PipeGeometry(segs, xs, Float64(roughness), Dict{Symbol,Float64}(taps))
    for (k, x) in p.taps
        0 <= x <= pipe_length(p) || throw(ArgumentError("tap $k at x=$x lies outside the pipe"))
    end
    return p
end

pipe_length(p::PipeGeometry) = p.xstart[end] + seglength(p.segments[end])

"Index of the segment containing `x` (interior points; ties go to the downstream segment)."
function segment_index(p::PipeGeometry, x)
    i = searchsortedlast(p.xstart, x)
    return clamp(i, 1, length(p.segments))
end

"Pipe diameter [m] at axial position `x`."
function diameter(p::PipeGeometry, x)
    i = segment_index(p, x)
    return segdiameter(p.segments[i], x - p.xstart[i])
end

"Flow area [m²] at `x`."
area(p::PipeGeometry, x) = π * diameter(p, x)^2 / 4

"Axial position of the named tap."
tap(p::PipeGeometry, name::Symbol) = p.taps[name]

"""
    straight_pipe(D, L; roughness=0.0)

Straight horizontal pipe. Taps `:inlet` and `:outlet`.
"""
straight_pipe(D, L; roughness=0.0) =
    PipeGeometry([Cylinder(L, D)]; roughness, taps=Dict(:inlet => 0.0, :outlet => Float64(L)))

"""
    iso_venturi(D; β=0.5, L_up=10D, L_down=20D, convergent_angle=21, divergent_angle=15,
                roughness=0.0)

Classical venturi tube to ISO 5167-4 installed in a straight pipe of diameter `D`:

  upstream pipe (`L_up`) → entrance cylinder (length D) → convergent (included angle
  21°) → throat (length d = βD) → divergent (included angle 7–15°) → downstream pipe
  (`L_down`).

Taps: `:upstream` (0.5 D upstream of the convergent, in the entrance cylinder),
`:throat` (mid-throat), `:inlet`, `:outlet`, and `:divergent_end`.
Angles are in degrees.
"""
function iso_venturi(D; β=0.5, L_up=10D, L_down=20D, convergent_angle=21.0,
                     divergent_angle=15.0, roughness=0.0)
    0.3 <= β <= 0.75 || @warn "ISO 5167-4 covers 0.3 ≤ β ≤ 0.75; got β = $β"
    7 <= divergent_angle <= 15 || @warn "ISO 5167-4 divergent angle is 7°–15°; got $divergent_angle°"
    d = β * D
    Lc = (D - d) / (2tand(convergent_angle / 2))
    Ld = (D - d) / (2tand(divergent_angle / 2))
    segs = Segment[
        Cylinder(L_up, D; tag=:pipe),
        Cylinder(D, D; tag=:entrance),
        Cone(Lc, D, d; tag=:convergent),
        Cylinder(d, d; tag=:throat),
        Cone(Ld, d, D; tag=:divergent),
        Cylinder(L_down, D; tag=:pipe),
    ]
    x_conv = L_up + D
    x_throat = x_conv + Lc
    x_div_end = x_throat + d + Ld
    taps = Dict(:inlet => 0.0,
                :upstream => x_conv - 0.5D,
                :throat => x_throat + d / 2,
                :divergent_end => x_div_end,
                :outlet => x_div_end + L_down)
    return PipeGeometry(segs; roughness, taps)
end
