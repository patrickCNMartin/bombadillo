@testset "samplers" begin
    rng = Xoshiro(11)
    N = 100_000
    # mean within 5 standard errors, variance within 5 %
    for λ in (0.0, 0.4, 6.0, 9.99, 10.0, 55.0, 3000.0)
        x = [B.rand_poisson(rng, λ) for _ in 1:N]
        @test all(>=(0), x)
        @test abs(mean(x) - λ) <= 5 * sqrt(max(λ, 1e-9) / N) + 1e-12
        λ > 0 && @test isapprox(var(x), λ; rtol = 0.05)
    end
    for (n, p) in ((0, 0.3), (7, 0.0), (7, 1.0), (12, 0.3), (60, 0.07),
                   (500, 0.3), (5000, 0.07), (400, 0.95), (32767, 0.5))
        x = [B.rand_binomial(rng, n, p) for _ in 1:N]
        @test all(k -> 0 <= k <= n, x)
        μ, σ² = n * p, n * p * (1 - p)
        @test abs(mean(x) - μ) <= 5 * sqrt(σ² / N) + 1e-12
        σ² > 0 && @test isapprox(var(x), σ²; rtol = 0.05)
    end
    @test B.logfactorial(0) == 0.0
    @test B.logfactorial(10) ≈ log(factorial(10))
    @test B.logfactorial(20) ≈ log(factorial(big(20))) rtol = 1e-12
    @test B.logfactorial(170) ≈ log(factorial(big(170))) rtol = 1e-12
end
