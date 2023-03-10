using Distributions
using ComradeOptimization
using OptimizationOptimJL
using StatsBase
using Plots
using LogDensityProblems
using LogDensityProblemsAD

load_ehtim()


@testset "bayes" begin
    _,vis, amp, lcamp, cphase = load_data()
    lklhd = RadioLikelihood(test_model, lcamp, cphase)

    prior = test_prior()

    post = Posterior(lklhd, prior)
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
           0.9,
           0.3,
           0.1,
           0.4,
           0.7,
           0.8,
           0.8]
    prob = OptimizationProblem(f, x0, nothing, lb=fill(0.001, 10), ub=fill(0.999,10))
    sol = solve(prob, LBFGS(); maxiters=100_000)

    xopt = transform(tpostc, sol)
    @test isapprox(xopt.f1/xopt.f2, 2.0, atol=1e-3)
    @test isapprox(xopt.σ1*2*sqrt(2*log(2)), μas2rad(40.0), rtol=1e-3)
    @test isapprox(xopt.σ1*xopt.τ1*2*sqrt(2*log(2)), μas2rad(20.0), rtol=1e-3)
    @test isapprox(xopt.ξ1, π/3, atol=1e-3)
    @test isapprox(xopt.σ2*2*sqrt(2*log(2)), μas2rad(20.0), atol=1e-3)
    @test isapprox(xopt.σ2*xopt.τ2*2*sqrt(2*log(2)), μas2rad(10.0), rtol=1e-3)
    @test isapprox(xopt.ξ2, π/6, atol=1e-3)
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

    @testset "LogDensityProblems" begin
        x0 = prior_sample(post)
        f0 = prior_sample(tpostf)
        c0 = prior_sample(tpostc)
        @test LogDensityProblems.logdensity(post, x0)   ≈ logdensityof(post, x0)
        @test LogDensityProblems.logdensity(tpostf, f0) ≈ logdensityof(tpostf, f0)
        @test LogDensityProblems.logdensity(tpostc, c0) ≈ logdensityof(tpostc, c0)

        @test LogDensityProblems.dimension(tpostf) == length(f0)
        @test LogDensityProblems.dimension(tpostc) == length(c0)

        @test LogDensityProblems.capabilities(typeof(post)) === LogDensityProblems.LogDensityOrder{0}()
        @test LogDensityProblems.capabilities(typeof(tpostf)) === LogDensityProblems.LogDensityOrder{0}()
        @test LogDensityProblems.capabilities(typeof(tpostc)) === LogDensityProblems.LogDensityOrder{0}()
    end



end

@testset "RadioLikelihood" begin
    _,vis, amp, lcamp, cphase, coh = load_data()
    tcache = TransformCache(coh)
    lklhd_amp = RadioLikelihood(test_model, amp)
    lklhd_cp = RadioLikelihood(test_model, cphase)
    lklhd_lc = RadioLikelihood(test_model, lcamp)
    lklhd_vis = RadioLikelihood(test_model, vis)
    lklhd_coh = RadioLikelihood(test_model_polarized, (;tcache), coh)


    prior = test_prior()
    post = Posterior(lklhd_amp, prior)
    x0 = prior_sample(post)

    lamp = logdensityof(lklhd_amp, x0)
    lcp  = logdensityof(lklhd_cp, x0)
    llc  = logdensityof(lklhd_lc, x0)
    lvis = logdensityof(lklhd_vis, x0)
    lcoh = logdensityof(lklhd_coh, x0)

    lklhd = MultiRadioLikelihood(lklhd_cp, lklhd_coh)
    show(lklhd)
    l = logdensityof(lklhd, x0)
    @test l ≈ lcp + lcoh

end

using ForwardDiff
using FiniteDifferences
@testset "Bayes Non-analytic ForwardDiff" begin
    _,vis, amp, lcamp, cphase = load_data()

    mfd = central_fdm(5,1)
    @testset "DFT" begin
        mt = (fovx = μas2rad(150.0), fovy=μas2rad(150.0), nx=256, ny=256, alg=DFTAlg())
        lklhd = RadioLikelihood(test_model2, mt, lcamp, cphase)

        prior = test_prior2()

        post = Posterior(lklhd, prior)
        tpostf = asflat(post)
        x0 = prior_sample(tpostf)

        @inferred logdensityof(tpostf, x0)
        ℓ = logdensityof(tpostf)
        gf = ForwardDiff.gradient(ℓ, x0)
        gn = FiniteDifferences.grad(mfd, ℓ, x0)
        @test isapprox(gf, gn[1], atol=1e-1, rtol=1e-5)
    end

    @testset "NFFT" begin
        mt = (fovx = μas2rad(150.0), fovy=μas2rad(150.0), nx=256, ny=256, alg=NFFTAlg(lcamp))
        lklhd = RadioLikelihood(test_model2, mt, lcamp, cphase)

        prior = test_prior2()

        post = Posterior(lklhd, prior)
        tpostf = asflat(post)
        x0 = prior_sample(tpostf)

        # @inferred logdensityof(tpostf, x0)
        ℓ = logdensityof(tpostf)
        gf = ForwardDiff.gradient(ℓ, x0)
        gn = FiniteDifferences.grad(mfd, ℓ, x0)
        @test isapprox(gf, gn[1], atol=1e-1, rtol=1e-5)
    end

    @testset "FFT" begin
        mt = (fovx = 10.0, fovy=10.0, nx=256, ny=256, alg=FFTAlg())
        lklhd = RadioLikelihood(test_model2, mt, lcamp, cphase)

        prior = test_prior2()

        post = Posterior(lklhd, prior)
        tpostf = asflat(post)
        x0 = prior_sample(tpostf)

        @inferred logdensityof(tpostf, x0)
        ℓ = logdensityof(tpostf)
        gf = ForwardDiff.gradient(ℓ, x0)
        gn = FiniteDifferences.grad(mfd, ℓ, x0)
        @test_broken isapprox(gf, gn[1], atol=1e-1, rtol=1e-5)
    end

end
