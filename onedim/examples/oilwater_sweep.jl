# Oil–water in a straight horizontal 10 cm pipe: pressure gradient, in-situ water holdup and
# flow pattern over mixture velocity and water cut (fully developed flow).
# Run from onedim/:  julia --project=. examples/oilwater_sweep.jl
using PipeFlow1D

D = 0.10
sys = OilWater()            # light crude (820 kg/m³, 3 mPa·s) + water, σ = 0.025 N/m
outdir = mkpath(joinpath(@__DIR__, "output"))

Um = [0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
wc = 0.0:0.05:1.0
t = sweep(sys, D; Um, wc)
write_csv(joinpath(outdir, "oilwater_sweep.csv"), t)

println("Inversion water cut (Arirachakaran): ", round(inversion_watercut(Arirachakaran(), sys), digits=3))
println("\n dp/dx [Pa/m]  (rows: Um [m/s], columns: water cut)")
print(rpad("Um\\wc", 8)); foreach(w -> print(lpad(w, 7)), wc[1:2:end]); println()
for u in Um
    print(rpad(u, 8))
    for w in wc[1:2:end]
        i = findfirst(k -> t.Um[k] == u && t.wc[k] == w, eachindex(t.Um))
        print(lpad(round(-t.dpdx_Pa_m[i], digits=1), 7))
    end
    println()
end

println("\nFlow pattern map (S = stratified, o = oil-in-water, w = water-in-oil, - = single phase)")
Umap = 0.1:0.1:3.0
wmap = 0.0:0.05:1.0
m = pattern_map(sys, D; Um=Umap, wc=wmap)
sym = Dict(STRATIFIED => 'S', DISPERSED_OW => 'o', DISPERSED_WO => 'w', SINGLE_OIL => '-', SINGLE_WATER => '-')
for i in reverse(axes(m, 1))
    println(rpad(round(Umap[i], digits=1), 5), join(sym[p] for p in m[i, :]))
end
println(" "^5, "wc: 0 → 1")
write_csv(joinpath(outdir, "pattern_map.csv"),
          (Um=[u for u in Umap for w in wmap], wc=[w for u in Umap for w in wmap],
           pattern=[string(m[i, j]) for i in axes(m, 1) for j in axes(m, 2)]))
println("\ntables written to examples/output/")
