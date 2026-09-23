"""
    Fluid(name, ρ, μ)

Incompressible Newtonian liquid: density `ρ` [kg/m³], dynamic viscosity `μ` [Pa·s].
"""
struct Fluid
    name::String
    ρ::Float64
    μ::Float64
end

"Fresh water at 20 °C."
const WATER = Fluid("water", 998.2, 1.0e-3)

"Very light crude oil (~40 °API) at pipeline temperature. Adjust to your crude assay."
const LIGHT_CRUDE = Fluid("light crude", 820.0, 3.0e-3)

"""
    OilWater(; oil=LIGHT_CRUDE, water=WATER, σ=0.025)

Oil–water system: the two liquids plus oil–water interfacial tension `σ` [N/m].
"""
Base.@kwdef struct OilWater
    oil::Fluid = LIGHT_CRUDE
    water::Fluid = WATER
    σ::Float64 = 0.025
end

"""
    Rates(Qo, Qw)
    Rates(; Qtotal, wc)

Volumetric production rates [m³/s] of oil `Qo` and water `Qw`.
The keyword form takes the total rate and water cut `wc = Qw/(Qo+Qw)`.
"""
struct Rates
    Qo::Float64
    Qw::Float64
    function Rates(Qo, Qw)
        (Qo >= 0 && Qw >= 0) || throw(ArgumentError("rates must be non-negative"))
        Qo + Qw > 0 || throw(ArgumentError("total rate must be positive"))
        new(Qo, Qw)
    end
end

function Rates(; Qtotal, wc)
    0 <= wc <= 1 || throw(ArgumentError("water cut must be in [0, 1]"))
    Rates((1 - wc) * Qtotal, wc * Qtotal)
end

total_rate(r::Rates) = r.Qo + r.Qw
watercut(r::Rates) = r.Qw / (r.Qo + r.Qw)
