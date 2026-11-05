using JET, AllocCheck

@testset "kernel quality" begin
    net = repressilator_network()
    sim = Simulation(net, Kinetics(net), 300; chunk = 128)
    set_level!(sim, :P, "lacI", 1:20, 200)
    args = (sim.nxt, sim.cur, sim.field, sim.net, sim.kin, sim.programs,
            sim.program, sim.rngs[1], 1:128)
    argtypes = Tuple{map(typeof, args)...}

    # type stability: no runtime dispatch in the kernel or samplers
    @test_opt target_modules = (Bombadillo,) B.step_cells!(args...)
    @test_opt target_modules = (Bombadillo,) B.rand_poisson(Xoshiro(1), 35.0)
    @test_opt target_modules = (Bombadillo,) B.rand_binomial(Xoshiro(1), 500, 0.3)
    @test_call target_modules = (Bombadillo,) B.step_cells!(args...)

    # zero heap allocations, proven statically (throw paths ignored)
    @test isempty(check_allocs(B.step_cells!, argtypes))
    @test isempty(check_allocs(B.rand_poisson, (Xoshiro, Float64)))
    @test isempty(check_allocs(B.rand_binomial, (Xoshiro, Int, Float64)))

    # and measured: a full single-threaded step allocates nothing
    step!(sim)
    @test (@allocated step!(sim)) == 0
end
