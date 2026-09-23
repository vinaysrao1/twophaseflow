"""
    PipeFlow1D

Steady, one-dimensional, incompressible model of single-phase and oil–water flow
in a horizontal pipe with an arbitrary piecewise-conical diameter profile
(straight pipe, ISO 5167-4 classical venturi, ...).

SI units throughout (m, s, kg, Pa, Pa·s, N/m).
"""
module PipeFlow1D

using Printf

include("units.jl")
include("fluids.jl")
include("rootfind.jl")
include("friction.jl")
include("geometry.jl")
include("losses.jl")
include("mixture.jl")
include("stratified.jl")
include("patterns.jl")
include("models.jl")
include("solver.jl")
include("venturi.jl")
include("io.jl")

export
    # units
    m3h, bpd, to_m3h, to_bpd,
    # fluids
    Fluid, OilWater, Rates, WATER, LIGHT_CRUDE, watercut, total_rate,
    # friction
    FrictionModel, Churchill, Colebrook, Haaland, Laminar, NoFriction, darcy, fanning,
    # geometry
    Segment, Cylinder, Cone, PipeGeometry, straight_pipe, iso_venturi,
    pipe_length, diameter, area, tap,
    # losses
    LocalLosses, CraneDiffuser, NoLocalLosses, loss_coefficient,
    # mixture
    EmulsionViscosity, Brinkman, KriegerDougherty, Taylor,
    InversionModel, Arirachakaran, FixedInversion, inversion_watercut, continuous_phase,
    mixture_density, mixture_viscosity,
    # stratified
    InterfacialFriction, FasterPhase, ConstantInterfacial, stratified_solutions, segment_geometry,
    # patterns
    FlowPattern, SINGLE_WATER, SINGLE_OIL, STRATIFIED, DISPERSED_OW, DISPERSED_WO,
    DispersionCriterion, HinzeBarnea, predict_pattern,
    # models
    TwoPhaseModel, Homogeneous, StratifiedTwoFluid, Mechanistic, LocalState, local_state,
    # solver
    Solution, solve, pressure_at, pressure_drop,
    # venturi
    VenturiReport, venturi_report, iso5167_venturi_C, venturi_flowrate, venturi_ideal_dp,
    # io
    write_csv, sweep, pattern_map

end
