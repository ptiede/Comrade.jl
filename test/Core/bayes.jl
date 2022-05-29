using Distributions
using AdvancedHMC
using Pathfinder
using ComradeGalactic
using GalacticBBO
using StatsBase
using Plots
using GalaticOptimJL

include(joinpath(@__DIR__, "test_util.jl"))

@testset "bayes" begin
    _,_, _, amp, cphase = load_data()
    lklhd = RadioLikelihood(amp, cphase)

    prior = test_prior()

    post = Posterior(lklhd, prior, test_model)
    prior_sample(post, 1)
    prior_sample(post, 2)
    tpostf = asflat(post)
    tpostc = ascube(post)
    tpost = Comrade.flatten(post)

    @inferred logdensityof(post, prior_sample(post, 1)[1])
    @inferred logdensityof(tpostf, rand(ndim))
    @inferred logdensityof(tpostc, rand(ndim))


    @test dimension(post) == dimension(tpost)
    @test dimension(post) == dimension(tpostf)
    @test dimension(post) == dimension(tpostc)

    f,tr = OptimizationFunction(tpostc, GalacticOptim.AutoForwardDiff{4}())
    prob = GalacticOptim.OptimizationProblem(f, prior_sample(tpostc), nothing, lb=fill(-5.0, 10), ub=fill(5.0,10))
    sol = solve(prob, BBO_adaptive_de_rand_1_bin(), tr; maxiters=500_000)
#
    xtest = transform(tr, sol)
    @test isapprox(xopt.f1-xopt.f2, 0.5, atol=1e-3)
    @test isapprox(xopt.σ1*2*sqrt(2*log(2)), μas2rad(40.0), rtol=1e-3)
    @test isapprox(xopt.σ1*xopt.τ1*2*sqrt(2*log(2)), μas2rad(20.0), rtol=1e-3)
    @test isapprox(-xopt.ξ1, π/3, atol=1e-3)
    @test isapprox(xopt.σ2*2*sqrt(2*log(2)), μas2rad(20.0), atol=1e-3)
    @test isapprox(xopt.σ2*xopt.τ2*2*sqrt(2*log(2)), μas2rad(10.0), rtol=1e-3)
    @test isapprox(-xopt.ξ2, π/6, atol=1e-3)
    @test isapprox(xopt.x, μas2rad(30.0), rtol=1e-3)
    @test isapprox(xopt.y, μas2rad(30.0), rtol=1e-3)


    #nchain, nstats = sample(post, Nested(dimension(post), 2_000); dlogz=0.01, progress=true)
    #@test isapprox(collect(values(nchain[end])), collect(values(xopt)), rtol=1e-2)
    #echain = Comrade.TupleVector(sample(nchain, Weights(nstats.weights), 500_000))
    #mn = Comrade.rmap(mean, echain)[1]
    #sn = Comrade.rmap(std, echain)[1]


end
