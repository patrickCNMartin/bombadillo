# Autocorrelation summary of a time series: first negative local minimum of the
# ACF (half period) and the following maximum (full period).
function acf_period(x::AbstractVector{<:Real}; maxlag = 800)
    y = Float64.(x) .- mean(x)
    v = sum(abs2, y)
    a = [sum(@view(y[1:end-l]) .* @view(y[1+l:end])) / v for l in 0:maxlag]
    i = findfirst(l -> a[l] < 0 && a[l] < a[l-1] && a[l] <= a[l+1], 3:maxlag)
    i === nothing && return (acf_min = minimum(a), acf_max = NaN, period = 0)
    lmin = i + 2
    w = lmin:min(maxlag + 1, 3lmin)
    lmax = w[argmax(a[w])]
    return (acf_min = a[lmin], acf_max = a[lmax], period = lmax - 1)
end

# Advance one tick processing chunks in reverse order (what a different thread
# schedule could do); must give the same result as step!.
function step_reversed!(sim)
    for k in reverse(eachindex(sim.rngs))
        B.step_cells!(sim.nxt, sim.cur, sim.field, sim.net, sim.kin, sim.programs,
                      sim.program, sim.rngs[k], B._chunk_range(sim, k))
    end
    sim.cur, sim.nxt = sim.nxt, sim.cur
    sim.tick += 1
    return sim
end

@testset "repressilator oscillates" begin
    net = repressilator_network()
    kin = Kinetics(net)
    for seed in 1:3
        sim = Simulation(net, kin, 1; seed)
        set_level!(sim, :P, "lacI", 1, 200)
        run!(sim, 500)                                   # burn-in
        tr = run!(sim, 6000; record = :P)
        @test size(tr.data) == (3, 1, 6001)
        @test tr.ticks[1] == 500 && tr.ticks[end] == 6500
        for g in 1:3
            s = acf_period(tr.data[g, 1, :])
            @test s.acf_min < -0.3          # anticorrelated half a period later
            @test s.acf_max > 0.2           # correlated again one period later
            @test 80 <= s.period <= 250     # ≈ 2 h, as in Elowitz & Leibler 2000
        end
        # ring order: peaks of lacI → cI → tetR ... check the phase lag sign
        x = Float64.(tr.data[1, 1, :]); y = Float64.(tr.data[2, 1, :])
        x .-= mean(x); y .-= mean(y)
        xc(l) = sum(@view(x[1:end-l]) .* @view(y[1+l:end]))
        # lacI represses tetR, so tetR follows lacI's troughs: tetR(t+l) is
        # anticorrelated with lacI(t) at small positive lag
        @test xc(15) < 0
    end
end

@testset "modality cascade: one tick per layer" begin
    # X activates Y at the binding layer; Z opens Y's chromatin
    edges = [Edge("X", "Y"; sign = 1, K = 1, n = 1, module_name = "bind"),
             Edge("Z", "W"; sign = 1, K = 1, n = 1, layer = :chromatin, module_name = "chrom")]
    net = RegNetwork(edges)
    # deterministic limits: certain ON, large transcription/translation, no leak
    kin = Kinetics(net; k_on = 1.0, k_off = 0.0, alpha = 1000, beta = 10, leak = 0,
                   delta_r = 0, delta_p = 0, kappa = 1.0)
    progs = Programs(net)
    progs.A0[net.gene_index["W"], 1] = 1f0
    sim = Simulation(net, kin, progs, [1])
    y, w = net.gene_index["Y"], net.gene_index["W"]
    # start with W closed so the chromatin edge has to open it
    sim.cur.A[w, 1] = 0f0
    set_level!(sim, :P, ["X", "Z"], 1, 30000)
    states = [cell(sim, 1)]
    for _ in 1:4
        step!(sim)
        push!(states, cell(sim, 1))
    end
    # binding input at t=0 → promoter ON at t=1 → RNA at t=2 → protein at t=3
    @test !states[1].S[y] && states[2].S[y]
    @test states[2].R[y] == 0 && states[3].R[y] > 0
    @test states[3].P[y] == 0 && states[4].P[y] > 0
    # chromatin input at t=0 → open at t=1 → promoter at t=2 → RNA at t=3
    @test states[1].A[w] == 0 && states[2].A[w] > 0.99
    @test !states[2].S[w] && states[3].S[w]
    @test states[3].R[w] == 0 && states[4].R[w] > 0
end

@testset "programs gate modules" begin
    edges = [Edge("a", "b"; sign = 1, K = 5, module_name = "on"),
             Edge("c", "d"; sign = 1, K = 5, module_name = "other")]
    net = RegNetwork(edges)
    progs = Programs(net, ["p1" => ["on"], "p2" => ["on", "other"]])
    sim = Simulation(net, Kinetics(net), progs, [1, 2]; seed = 4)
    set_level!(sim, :P, ["a", "c"], [1, 2], 500)
    run!(sim, 200)
    d = net.gene_index["d"]
    @test sim.cur.P[d, 1] == 0                 # module absent and chromatin closed
    @test sim.cur.P[d, 2] > 0
    @test sim.cur.P[net.gene_index["b"], 1] > 0
    @test cell(sim, 1).program == "p1" && cell(sim, 2).program == "p2"
end

@testset "reproducibility and thread-schedule invariance" begin
    net = repressilator_network()
    kin = Kinetics(net)
    mk() = (s = Simulation(net, kin, 700; seed = 9, chunk = 100);
            set_level!(s, :P, "lacI", 1:50, 300); s)
    a, b = mk(), mk()
    run!(a, 50)
    for _ in 1:50
        step_reversed!(b)
    end
    @test a.cur.P == b.cur.P && a.cur.R == b.cur.R && a.cur.S == b.cur.S && a.cur.A == b.cur.A
    @test a.tick == b.tick == 50

    c = Simulation(net, kin, 700; seed = 10, chunk = 100)
    set_level!(c, :P, "lacI", 1:50, 300)
    run!(c, 50)
    @test c.cur.P != a.cur.P

    # reset! replays exactly
    p_before = copy(a.cur.P)
    reset!(a); set_level!(a, :P, "lacI", 1:50, 300); run!(a, 50)
    @test a.cur.P == p_before
end

@testset "Int16 saturation" begin
    net = RegNetwork([Edge("a", "b"; sign = 1)])
    kin = Kinetics(net; k_on = 1, k_off = 0, alpha = 5000, beta = 50, delta_r = 0, delta_p = 0)
    sim = Simulation(net, kin, 2)
    run!(sim, 30)
    @test all(==(typemax(Int16)), sim.cur.P)
    @test all(==(typemax(Int16)), sim.cur.R)
end
