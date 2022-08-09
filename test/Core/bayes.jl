using Distributions
using ComradeOptimization
using OptimizationOptimJL
using StatsBase
using Plots

load_ehtim()


@testset "bayes" begin
    _,vis, amp, lcamp, cphase = load_data()
    lklhd = RadioLikelihood(lcamp, cphase)

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

    f = OptimizationFunction(tpostc, Optimization.AutoForwardDiff{4}())
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
    prob = OptimizationProblem(f, x0, nothing, lb=fill(0.001, 10), ub=fill(0.999,10))
    sol = solve(prob, LBFGS(); maxiters=100_000)

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

    mopt = test_model(xopt)

    @testset "Plot model" begin
        plot(mopt, vis)
        plot(mopt, amp)
        plot(mopt, lcamp)
        plot(mopt, cphase)
    end

    @testset "Plot residuals" begin
        residual(mopt, vis)
        residual(mopt, amp)
        residual(mopt, lcamp)
        residual(mopt, cphase)
    end



end
