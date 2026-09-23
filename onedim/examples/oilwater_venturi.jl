# Oil–water through an ISO 5167-4 venturi (β = 0.5) in a 10 cm pipe.
# Pressure profile, venturi differential pressure, effective discharge coefficient and
# permanent loss versus water cut; and flow-rate inference from Δp (metering).
# Run from onedim/:  julia --project=. examples/oilwater_venturi.jl
using PipeFlow1D

D, β = 0.10, 0.5
sys = OilWater()
v = iso_venturi(D; β)
outdir = mkpath(joinpath(@__DIR__, "output"))

Q = m3h(40.0)    # ≈ 1.41 m/s mixture velocity, ≈ 6000 bbl/d
println("Q = 40 m³/h through ISO venturi, β = $β")
println(rpad("wc", 6), rpad("pattern(up)", 14), rpad("Δp [Pa]", 10), rpad("C_eff", 8),
        rpad("loss [%]", 10), rpad("Q_meter/Q", 10))
wcs = 0.0:0.1:1.0
rows = (wc=Float64[], dp_Pa=Float64[], C=Float64[], loss_frac=Float64[], Re_D=Float64[], Q_ratio=Float64[])
for w in wcs
    s = solve(v, sys, Rates(Qtotal=Q, wc=w))
    r = venturi_report(s)
    # a meter using ISO C = 0.995 and the homogeneous density
    Qm = venturi_flowrate(r.dp, r.ρm, D, β; C=iso5167_venturi_C(:machined).C)
    println(rpad(w, 6), rpad(string(s.pattern[1]), 14), rpad(round(r.dp, digits=1), 10),
            rpad(round(r.C, digits=4), 8), rpad(round(100r.loss_fraction, digits=1), 10), round(Qm / Q, digits=4))
    push!(rows.wc, w); push!(rows.dp_Pa, r.dp); push!(rows.C, r.C)
    push!(rows.loss_frac, r.loss_fraction); push!(rows.Re_D, r.Re_D); push!(rows.Q_ratio, Qm / Q)
    w == 0.3 && write_csv(joinpath(outdir, "venturi_oilwater_wc30_profile.csv"), s)
end
write_csv(joinpath(outdir, "venturi_oilwater_vs_wc.csv"), rows)
println("\nC_eff < ISO value at low Re_D reflects the high emulsion viscosity near inversion.")
println("tables written to examples/output/")
