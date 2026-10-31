using Comrade, MCMCDiagnosticTools, Random
using Test

@testset "ComradeMCMCDiagnosticToolsExt" begin
    rng = Random.Xoshiro(3)
    n, k = 200, 3
    function chain(μ)
        θs = [(; σ = μ + randn(rng), a = randn(rng, 2, 3), ρ = (randn(rng), randn(rng))) for _ in 1:n]
        return PosteriorSamples(θs, nothing)
    end
    ps = reduce(hcat, [chain(0.0) for _ in 1:k])
    cube(x) = reshape(x, size(x, 1), size(x, 2), 1)

    e = ess(ps)
    @test e.σ ≈ only(ess(cube(ps.σ)))
    @test size(e.a) == (2, 3)
    @test e.a[2, 3] ≈ only(ess(cube(getindex.(ps.a, 2, 3))))
    @test e.ρ[2] ≈ only(ess(cube(last.(ps.ρ))))
    @test ess(ps; kind = :tail).σ ≈ only(ess(cube(ps.σ); kind = :tail))
    @test rhat(ps).σ ≈ only(rhat(cube(ps.σ)))
    @test mcse(ps).a[1, 2] ≈ only(mcse(cube(getindex.(ps.a, 1, 2))))

    single = chain(0.0)
    @test ess(single).σ ≈ only(ess(reshape(single.σ, n, 1, 1)))
    @test rhat(hcat(chain(0.0), chain(5.0))).σ > 1.5
end
