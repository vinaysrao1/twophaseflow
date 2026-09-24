# Laminar Poiseuille flow in a periodic pipe with WaterLily BDIM: flow-rate error vs resolution.
using WaterLily, Printf
function run(R; ν_over_R=0.1, T=Float64)
    c = R + 4; n = 2R + 8
    ν = ν_over_R * R
    G = 4ν / R^2                       # gives u_max = 1 for exact Poiseuille
    body = AutoBody((x, t) -> R - sqrt((x[2] - c)^2 + (x[3] - c)^2))
    sim = Simulation((8, n, n), (0, 0, 0), R; U=1, ν, body, perdir=(1,),
                     g=(i, x, t) -> i == 1 ? G : zero(G), T, mem=Array)
    sim_step!(sim, 25.0)               # 25 R/U ≈ 2.5 R²/ν: decay e^-14
    u = sim.flow.u[3, :, :, 1]
    Q = sum(u[2:end-1, 2:end-1])
    Qex = π * R^2 / 2
    Reff = (8ν * Q / (π * G))^(1 / 4)
    @printf("R = %3d cells: Q/Q_exact = %.4f   u_max = %.4f   R_eff - R = %+.3f cells   steps = %d\n",
            R, Q / Qex, maximum(u), Reff - R, length(sim.flow.Δt))
end
for R in (8, 16, 32); run(R); end
