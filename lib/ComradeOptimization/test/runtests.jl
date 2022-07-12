using Comrade, ComradeOptimization, OptimizationOptimJL, Distributions
using Test

load_ehtim()
include(joinpath(@__DIR__, "../../../test/test_util.jl"))

@testset "ComradeOptimization.jl" begin
    m, vis, amp, lcamp, cphase = load_data()
    prior = test_prior()
    lklhd = RadioLikelihood(lcamp, cphase)
    post = Posterior(lklhd, prior, test_model)

    tpost = ascube(post)
    f = OptimizationFunction(tpost, GalacticOptim.AutoForwardDiff{4}())
    x0 = [ 0.1,
           0.4,
           0.5,
           0.1,
           0.3,
           0.1,
           0.4,
           0.3,
           0.8,
           0.8]
    prob = GalacticOptim.OptimizationProblem(f, x0, nothing, lb=fill(0.001, 10), ub=fill(0.999,10))
    sol = solve(prob, LBFGS(); maxiters=100_000)

    xopt = transform(tpost, sol)
    @test isapprox(xopt.f1/xopt.f2, 2.0, atol=1e-3)
    @test isapprox(xopt.σ1*2*sqrt(2*log(2)), μas2rad(40.0), rtol=1e-3)
    @test isapprox(xopt.σ1*xopt.τ1*2*sqrt(2*log(2)), μas2rad(20.0), rtol=1e-3)
    @test isapprox(-xopt.ξ1, π/3, atol=1e-3)
    @test isapprox(xopt.σ2*2*sqrt(2*log(2)), μas2rad(20.0), atol=1e-3)
    @test isapprox(xopt.σ2*xopt.τ2*2*sqrt(2*log(2)), μas2rad(10.0), rtol=1e-3)
    @test isapprox(-xopt.ξ2, π/6, atol=1e-3)
    @test isapprox(xopt.x, μas2rad(30.0), rtol=1e-3)
    @test isapprox(xopt.y, μas2rad(30.0), rtol=1e-3)

end
