# Single-phase water and light crude in a straight 10 cm pipe and through an ISO venturi.
# Run from onedim/:  julia --project=. examples/single_phase.jl
using PipeFlow1D

D = 0.10
outdir = mkpath(joinpath(@__DIR__, "output"))

println("Straight pipe, L = 100 m, fully developed friction")
println(rpad("fluid", 14), rpad("Q [m³/h]", 10), rpad("U [m/s]", 9), rpad("Re", 10), rpad("f_D", 9), "Δp [kPa]")
for fl in (WATER, LIGHT_CRUDE), Qh in (5.0, 28.3, 85.0)
    Q = m3h(Qh)
    s = solve(straight_pipe(D, 100.0), fl, Q)
    U = Q / (π * D^2 / 4)
    Re = fl.ρ * U * D / fl.μ
    println(rpad(fl.name, 14), rpad(Qh, 10), rpad(round(U, digits=3), 9), rpad(round(Int, Re), 10),
            rpad(round(darcy(Churchill(), Re, 0.0), digits=5), 9), round(pressure_drop(s) / 1e3, digits=3))
end

println("\nISO 5167-4 venturi, β = 0.5, water at 1 m/s")
v = iso_venturi(D; β=0.5)
s = solve(v, WATER, m3h(28.3))
display(venturi_report(s)); println()
write_csv(joinpath(outdir, "venturi_water_profile.csv"), s)
println("profile written to examples/output/venturi_water_profile.csv")
