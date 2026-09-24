using WaterLily, Printf
function tp(n; T=Float32)
    R = n ÷ 2 - 4; c = n / 2
    body = AutoBody((x, t) -> R - sqrt((x[2] - c)^2 + (x[3] - c)^2))
    sim = Simulation((2n, n, n), (1, 0, 0), R; U=1, ν=1e-3, body, perdir=(1,), T, mem=Array)
    sim_step!(sim; remeasure=false)
    k = 20; t0 = time(); for _ in 1:k; sim_step!(sim; remeasure=false); end; dt = time() - t0
    N = 2n * n * n
    @printf("grid %d×%d×%d (%.2e cells), threads=%d: %.3e cell-steps/s\n", 2n, n, n, N, Threads.nthreads(), N * k / dt)
end
tp(48); tp(96)
