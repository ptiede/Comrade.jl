using Distributions
using ComradeGalactic
using GalacticBBO
using StatsBase
using Plots

load_ehtim()

include(joinpath(@__DIR__, "../test_util.jl"))

@testset "bayes" begin
    _,_, _, amp, cphase = load_data()
    lklhd = RadioLikelihood(amp, cphase)

    prior = test_prior()

    post = Posterior(lklhd, prior, test_model)
    prior_sample(post)
    prior_sample(post, 2)
    tpostf = asflat(post)
    tpostc = ascube(post)
    tpost = Comrade.flatten(post)

    ndim = dimension(post)
    @inferred logdensityof(post, prior_sample(post))
    @inferred logdensityof(tpostf, rand(ndim))
    @inferred logdensityof(tpostc, rand(ndim))


    @test dimension(post) == dimension(tpost)
    @test dimension(post) == dimension(tpostf)
    @test dimension(post) == dimension(tpostc)

    f = OptimizationFunction(tpostc, GalacticOptim.AutoForwardDiff{4}())
    prob = GalacticOptim.OptimizationProblem(f, prior_sample(tpostc), nothing, lb=fill(0.001, 10), ub=fill(0.999,10))
    sol = solve(prob, BBO_adaptive_de_rand_1_bin(); maxiters=100_000)
#
    xopt = transform(tpostc, sol)
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
