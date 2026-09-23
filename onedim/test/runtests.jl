using PipeFlow1D
using PipeFlow1D: brent, stratified_eval, wall_shear, build_grid
using Test

const D = 0.1
const A = π * D^2 / 4

@testset "PipeFlow1D" begin

@testset "units and inputs" begin
    @test to_m3h(m3h(12.5)) ≈ 12.5
    @test to_bpd(bpd(1000)) ≈ 1000
    @test bpd(1) ≈ 0.158987294928 / 86400
    r = Rates(Qtotal=0.01, wc=0.3)
    @test r.Qw ≈ 0.003 && r.Qo ≈ 0.007
    @test watercut(r) ≈ 0.3 && total_rate(r) ≈ 0.01
    @test_throws ArgumentError Rates(-1.0, 1.0)
    @test_throws ArgumentError Rates(0.0, 0.0)
    @test_throws ArgumentError Rates(Qtotal=1.0, wc=1.2)
end

@testset "Brent root finder" begin
    @test brent(x -> x^3 - 2, 0.0, 2.0) ≈ cbrt(2) atol = 1e-12
    @test brent(cos, 1.0, 2.0) ≈ π / 2 atol = 1e-12
    @test_throws ArgumentError brent(x -> x^2 + 1, -1.0, 1.0)
end

@testset "friction factors" begin
    for Re in (10.0, 500.0, 2000.0)
        @test darcy(Laminar(), Re, 0.0) ≈ 64 / Re
        @test darcy(Colebrook(), Re, 0.0) ≈ 64 / Re
        @test darcy(Churchill(), Re, 0.0) ≈ 64 / Re rtol = 0.01
    end
    # Colebrook satisfies its own implicit equation
    for Re in (1e4, 1e5, 1e6, 1e7), rr in (0.0, 1e-4, 1e-3)
        f = darcy(Colebrook(), Re, rr)
        @test 1 / sqrt(f) ≈ -2log10(rr / 3.7 + 2.51 / (Re * sqrt(f))) rtol = 1e-10
        @test darcy(Churchill(), Re, rr) ≈ f rtol = 0.03
        @test darcy(Haaland(), Re, rr) ≈ f rtol = 0.03
    end
    # Moody-chart reference values (smooth pipe)
    @test darcy(Colebrook(), 1e5, 0.0) ≈ 0.01799 rtol = 2e-3
    @test darcy(Colebrook(), 1e6, 0.0) ≈ 0.01165 rtol = 2e-3
    @test darcy(NoFriction(), 1e5, 0.0) == 0
    # Churchill is continuous through transition
    fs = [darcy(Churchill(), Re, 0.0) for Re in 1500:10:5000]
    @test maximum(abs.(diff(fs)) ./ fs[1:end-1]) < 0.02
    # shear → 0 as velocity → 0 (laminar limit 8μu/D)
    @test wall_shear(Churchill(), 1000.0, 1e-3, 1e-6, 0.1, 0.0) ≈ 8e-3 * 1e-6 / 0.1 rtol = 0.01
    @test wall_shear(Churchill(), 1000.0, 1e-3, 0.0, 0.1, 0.0) == 0
end

@testset "geometry" begin
    p = straight_pipe(D, 50.0; roughness=4.5e-5)
    @test pipe_length(p) == 50.0
    @test diameter(p, 17.3) == D
    @test p.roughness == 4.5e-5

    β = 0.5
    v = iso_venturi(D; β, L_up=1.0, L_down=2.0)
    d = β * D
    Lc = (D - d) / (2tand(10.5))
    Ld = (D - d) / (2tand(7.5))
    @test pipe_length(v) ≈ 1.0 + D + Lc + d + Ld + 2.0
    @test diameter(v, tap(v, :upstream)) == D
    @test diameter(v, tap(v, :throat)) ≈ d
    @test tap(v, :upstream) ≈ 1.0 + D - D / 2
    @test tap(v, :throat) ≈ 1.0 + D + Lc + d / 2
    # continuous diameter profile
    xs = range(0, pipe_length(v); length=20001)
    Ds = [diameter(v, x) for x in xs]
    @test maximum(abs.(diff(Ds))) < 1e-4
    @test minimum(Ds) ≈ d
    @test_throws ArgumentError PipeGeometry([Cylinder(1.0, D)]; taps=Dict(:x => 2.0))
    # grid contains segment boundaries and taps
    g = build_grid(v, 0.01)
    @test issorted(g) && tap(v, :throat) in g && (1.0 + D) in g
end

@testset "local losses" begin
    @test loss_coefficient(CraneDiffuser(), Cylinder(1.0, D)) === nothing
    @test loss_coefficient(CraneDiffuser(), Cone(0.1, D, 0.05)) === nothing
    c = Cone((D - 0.05) / (2tand(7.5)), 0.05, D)
    @test loss_coefficient(CraneDiffuser(), c) ≈ 2.6sind(7.5) * (1 - 0.25)^2 rtol = 1e-10
    @test loss_coefficient(CraneDiffuser(), Cone(0.001, 0.05, D)) ≈ (1 - 0.25)^2   # sudden
    @test loss_coefficient(NoLocalLosses(), c) === nothing
    @test loss_coefficient(CraneDiffuser(factor=0.5), c) ≈ 0.5loss_coefficient(CraneDiffuser(), c)
end

@testset "single-phase straight pipe" begin
    L = 100.0
    for (Q, fm) in ((m3h(28.0), Colebrook()), (m3h(5.0), Churchill()), (1e-5, Laminar()))
        s = solve(straight_pipe(D, L), WATER, Q; friction=fm)
        U = Q / A
        Re = WATER.ρ * U * D / WATER.μ
        @test pressure_drop(s) ≈ darcy(fm, Re, 0.0) * L / D * WATER.ρ * U^2 / 2 rtol = 1e-10
        @test all(==(SINGLE_WATER), s.pattern)
    end
    # Hagen–Poiseuille: Δp = 128 μ L Q / (π D⁴)
    Q = 1e-5
    s = solve(straight_pipe(D, L), WATER, Q; friction=Laminar())
    @test pressure_drop(s) ≈ 128 * WATER.μ * L * Q / (π * D^4) rtol = 1e-10
    # linear pressure profile
    @test pressure_at(s, L / 2) ≈ pressure_drop(s) * -0.5 rtol = 1e-10
end

@testset "frictionless venturi recovers Bernoulli exactly" begin
    v = iso_venturi(D; L_up=0.5, L_down=0.5)
    Q = m3h(30.0)
    s = solve(v, WATER, Q; friction=NoFriction(), losses=NoLocalLosses())
    for k in eachindex(s.x)
        U1, U = Q / A, Q / (π * s.D[k]^2 / 4)
        @test s.p[k] ≈ WATER.ρ / 2 * (U1^2 - U^2) atol = 1e-8
    end
    @test pressure_drop(s) ≈ 0 atol = 1e-8
    r = venturi_report(s)
    @test r.C ≈ 1 rtol = 1e-10
    @test r.dp ≈ venturi_ideal_dp(Q, WATER.ρ, D, 0.5) rtol = 1e-10
    @test r.permanent_loss ≈ 0 atol = 1e-8
end

@testset "venturi with friction and diffuser loss" begin
    v = iso_venturi(D)
    for Um in (1.0, 3.0, 5.0)      # Re_D = 1e5 … 5e5
        s = solve(v, WATER, Um * A)
        r = venturi_report(s)
        @test 0.975 < r.C < 1.0                 # ISO 5167-4: 0.984–0.995
        @test 0.05 < r.loss_fraction < 0.30      # typical 5–20 %, 15° diffuser upper end
    end
    # 7° diffuser loses less than 15°
    s7 = solve(iso_venturi(D; divergent_angle=7.0), WATER, 2A)
    s15 = solve(iso_venturi(D; divergent_angle=15.0), WATER, 2A)
    @test venturi_report(s7).permanent_loss < venturi_report(s15).permanent_loss
    # metering round trip with the model's own C
    s = solve(v, WATER, 2A)
    r = venturi_report(s)
    @test venturi_flowrate(r.dp, WATER.ρ, D, 0.5; C=r.C) ≈ 2A rtol = 1e-10
    # grid convergence
    s1 = solve(v, WATER, 2A; dx=0.005)
    s2 = solve(v, WATER, 2A; dx=0.0025)
    @test venturi_report(s1).dp ≈ venturi_report(s2).dp rtol = 1e-4
    @test pressure_drop(s1) ≈ pressure_drop(s2) rtol = 1e-4
    @test iso5167_venturi_C(:machined).C == 0.995
end

sys = OilWater()

@testset "mixture properties" begin
    @test mixture_density(sys, 0.3) ≈ 0.3 * sys.water.ρ + 0.7 * sys.oil.ρ
    wci = inversion_watercut(Arirachakaran(), sys)
    @test wci ≈ 0.5 - 0.1108log10(3.0)
    @test continuous_phase(Arirachakaran(), sys, 0.9) === :water
    @test continuous_phase(Arirachakaran(), sys, 0.1) === :oil
    @test mixture_viscosity(Brinkman(), FixedInversion(0.5), sys, 0.8) ≈ sys.water.μ * 0.8^-2.5
    @test mixture_viscosity(Brinkman(), FixedInversion(0.5), sys, 0.2) ≈ sys.oil.μ * 0.8^-2.5
    @test mixture_viscosity(Brinkman(), Arirachakaran(), sys, 1.0) == sys.water.μ
    @test mixture_viscosity(Brinkman(), Arirachakaran(), sys, 0.0) == sys.oil.μ
    # viscosity peaks at inversion
    wcs = 0.0:0.01:1.0
    μs = [mixture_viscosity(Brinkman(), Arirachakaran(), sys, w) for w in wcs]
    @test abs(wcs[argmax(μs)] - wci) <= 0.011
    # Taylor → Einstein 2.5 for rigid drops
    @test PipeFlow1D.emulsion_viscosity(Taylor(), 1.0, 1e9, 0.01) ≈ 1.025 rtol = 1e-6
    @test PipeFlow1D.emulsion_viscosity(KriegerDougherty(), 1.0, 1.0, 0.01) ≈ 1.025 rtol = 2e-3
end

@testset "stratified two-fluid model" begin
    g = segment_geometry(π / 2, D)
    @test g.Aw ≈ A / 2 && g.Si ≈ D && g.h ≈ D / 2
    g = segment_geometry(π - 1e-9, D)
    @test g.Aw ≈ A rtol = 1e-8

    for Um in (0.2, 0.5, 1.0), wc in (0.05, 0.3, 0.6, 0.95)
        Qo, Qw = (1 - wc) * Um * A, wc * Um * A
        sols = stratified_solutions(sys, Qo, Qw, D)
        @test !isempty(sols)
        for s in sols
            # mass conservation
            gg = segment_geometry(s.φ, D)
            @test s.uw * gg.Aw ≈ Qw rtol = 1e-10
            @test s.uo * gg.Ao ≈ Qo rtol = 1e-10
            if !s.locked
                # both phase momentum equations give the same dp/dx
                e = stratified_eval(s.φ, sys, Qo, Qw, D, 0.0, Churchill(), FasterPhase())
                dw = (-e.τw * gg.Sw + e.τi * gg.Si) / gg.Aw
                do_ = (-e.τo * gg.So - e.τi * gg.Si) / gg.Ao
                @test dw ≈ do_ rtol = 1e-6
                @test dw ≈ s.dpdx rtol = 1e-6
            end
            @test s.dpdx < 0
        end
    end

    # More viscous oil moves slower → accumulates (oil holdup > input oil fraction)
    visc = OilWater(oil=Fluid("viscous oil", 850.0, 50e-3))
    s = stratified_solutions(visc, 0.5 * 0.3A, 0.5 * 0.3A, D)[1]
    @test 1 - s.Hw > 0.5
    # Less viscous water... with equal viscosities the slip is small
    same = OilWater(oil=Fluid("water-like", 998.2, 1e-3), water=WATER)
    s = stratified_solutions(same, 0.5A, 0.5A, D)[1]
    @test s.Hw ≈ 0.5 atol = 0.02

    # Vanishing phase: pressure gradient tends to the single-phase value
    Q = 1.0 * A
    so = stratified_solutions(sys, Q * (1 - 1e-5), Q * 1e-5, D)[1]
    oil_only = local_state(Homogeneous(), sys, Q, 0.0, D)
    @test so.dpdx ≈ oil_only.dpdx rtol = 0.02
    sw = stratified_solutions(sys, Q * 1e-5, Q * (1 - 1e-5), D)[end]
    water_only = local_state(Homogeneous(), sys, 0.0, Q, D)
    @test sw.dpdx ≈ water_only.dpdx rtol = 0.02

    @test_throws ArgumentError stratified_solutions(sys, 0.0, 1e-3, D)
end

@testset "flow patterns" begin
    # low velocity: stratified; high velocity: dispersed; continuous phase by inversion
    @test predict_pattern(HinzeBarnea(), Arirachakaran(), sys, 0.2A * 0.5, 0.2A * 0.5, D)[1] == STRATIFIED
    @test predict_pattern(HinzeBarnea(), Arirachakaran(), sys, 3A * 0.2, 3A * 0.8, D)[1] == DISPERSED_OW
    @test predict_pattern(HinzeBarnea(), Arirachakaran(), sys, 3A * 0.8, 3A * 0.2, D)[1] == DISPERSED_WO
    @test predict_pattern(HinzeBarnea(), Arirachakaran(), sys, 0.0, A, D)[1] == SINGLE_WATER
    @test predict_pattern(HinzeBarnea(), Arirachakaran(), sys, A, 0.0, D)[1] == SINGLE_OIL
    # monotone: once dispersed, stays dispersed at higher velocity
    Um = 0.1:0.05:4.0
    m = pattern_map(sys, D; Um, wc=[0.2, 0.5, 0.8])
    for j in axes(m, 2)
        disp = m[:, j] .!= STRATIFIED
        k = findfirst(disp)
        @test k !== nothing && all(disp[k:end])
        @test 0.5 < Um[k] < 3.0     # transition of the order of 1–2 m/s for light oil, D = 10 cm
    end
end

@testset "two-phase models" begin
    Q = 0.5A
    for wc in (0.2, 0.7)
        r = Rates(Qtotal=Q, wc=wc)
        h = local_state(Homogeneous(), sys, r.Qo, r.Qw, D)
        @test h.Hw == wc
        @test h.J ≈ mixture_density(sys, wc) * (Q / A)^2
        st = local_state(StratifiedTwoFluid(), sys, r.Qo, r.Qw, D)
        @test st.pattern == STRATIFIED
        # momentum flux consistent with holdup
        @test st.J ≈ sys.water.ρ * st.Hw * st.uw^2 + sys.oil.ρ * (1 - st.Hw) * st.uo^2 rtol = 1e-10
    end
    # single-phase limits of the mechanistic model
    a = local_state(Mechanistic(), sys, 0.0, Q, D)
    b = local_state(Homogeneous(), OilWater(), 0.0, Q, D)
    @test a.dpdx == b.dpdx && a.pattern == SINGLE_WATER
    # tiny water cut → close to oil only
    near = local_state(Mechanistic(), sys, Q * (1 - 1e-6), Q * 1e-6, D)
    oil = local_state(Mechanistic(), sys, Q, 0.0, D)
    @test near.dpdx ≈ oil.dpdx rtol = 0.02
end

@testset "two-phase straight pipe and venturi" begin
    L = 20.0
    r = Rates(Qtotal=0.5A, wc=0.4)
    s = solve(straight_pipe(D, L), sys, r)
    st = local_state(Mechanistic(), sys, r.Qo, r.Qw, D)
    @test pressure_drop(s) ≈ -st.dpdx * L rtol = 1e-10
    @test all(==(STRATIFIED), s.pattern)

    v = iso_venturi(D)
    for wc in (0.1, 0.5, 0.9), Um in (0.5, 2.0)
        r = Rates(Qtotal=Um * A, wc=wc)
        s = solve(v, sys, r)
        rep = venturi_report(s)
        @test 0.95 < rep.C < 1.0
        @test rep.permanent_loss > 0
        # inside the venturi the default model is homogeneous
        k = findfirst(>(tap(v, :throat)), s.x)
        @test s.pattern[k] in (DISPERSED_OW, DISPERSED_WO)
        @test s.Hw[k] ≈ wc
        # frictionless homogeneous venturi recovers Bernoulli on mixture density
        s0 = solve(v, sys, r; model=Homogeneous(), friction=NoFriction(), losses=NoLocalLosses())
        @test venturi_report(s0).C ≈ 1 rtol = 1e-10
    end
end

@testset "sweep and CSV" begin
    t = sweep(sys, D; Um=[0.5, 2.0], wc=[0.0, 0.5, 1.0])
    @test length(t.Um) == 6
    @test t.pattern[1] == "SINGLE_OIL" && t.pattern[3] == "SINGLE_WATER"
    @test all(<(0), t.dpdx_Pa_m)
    path = tempname() * ".csv"
    write_csv(path, t)
    lines = readlines(path)
    @test length(lines) == 7 && startswith(lines[1], "Um,wc,")
    s = solve(straight_pipe(D, 1.0), WATER, 0.01)
    write_csv(path, s)
    @test length(readlines(path)) == length(s.x) + 1
end

end
