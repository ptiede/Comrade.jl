using Distributions
using ComradeOptimization
using OptimizationGCMAES
using StatsBase
using Plots
using LogDensityProblems
using LogDensityProblemsAD
using Pyehtim
using Zygote



@testset "bayes" begin
    _,vis, amp, lcamp, cphase = load_data()
    lklhd_cl = RadioLikelihood(test_model, lcamp, cphase)
    lklhd = RadioLikelihood(test_model, vis)
    show(lklhd)
    show(lklhd_cl)

    prior = test_prior()

    post = Posterior(lklhd, prior)
    prior_sample(post)
    prior_sample(post, 2)
    tpostf = asflat(post)
    tpostc = ascube(post)
    show(post)
    show(tpostf)
    show(tpostc)
    tpost = Comrade.flatten(post)

    ndim = dimension(post)
    @inferred logdensityof(post, prior_sample(post))
    @inferred logdensityof(tpostf, rand(ndim))
    @inferred logdensityof(tpostc, rand(ndim))


    @test dimension(post) == dimension(tpost)
    @test dimension(post) == dimension(tpostf)
    @test dimension(post) == dimension(tpostc)

    f = OptimizationFunction(tpostf, Optimization.AutoZygote())
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
    prob = OptimizationProblem(f, x0, nothing; lb=fill(-5.0, ndim), ub = fill(5.0, ndim))
    sol = solve(prob, GCMAESOpt(); maxiters=10_000)

    xopt = transform(tpostf, sol)
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
    tcache = ResponseCache(coh)
    lklhd_amp = RadioLikelihood(test_model, amp)
    lklhd_cp = RadioLikelihood(test_model, cphase)
    lklhd_lc = RadioLikelihood(test_model, lcamp)
    lklhd_vis = RadioLikelihood(test_model, vis)
    lklhd_coh = RadioLikelihood(test_skymodel_polarized, test_instrumentmodel_polarized, coh;
                                skymeta = (;lp = 0.1), instrumentmeta=(;tcache))


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

@testset "simulate_obs" begin
    _,vis, amp, lcamp, cphase, coh = load_data()
    tcache = ResponseCache(coh)
    lklhd_amp = RadioLikelihood(test_model, amp)
    lklhd_cp = RadioLikelihood(test_model, cphase)
    lklhd_lc = RadioLikelihood(test_model, lcamp)
    lklhd_vis = RadioLikelihood(test_model, vis)
    lklhd_coh = RadioLikelihood(test_skymodel_polarized, test_instrumentmodel_polarized, coh;
                                skymeta = (;lp = 0.1), instrumentmeta=(;tcache))


    prior = test_prior()

    simulate_observation(Posterior(lklhd_amp, prior), rand(prior))
    simulate_observation(Posterior(lklhd_cp,  prior), rand(prior))
    simulate_observation(Posterior(lklhd_lc,  prior), rand(prior))
    simulate_observation(Posterior(lklhd_vis, prior), rand(prior))
    simulate_observation(Posterior(lklhd_coh, prior), rand(prior))

    lklhd_all = RadioLikelihood(test_model, amp, cphase, lcamp, vis)

    simulate_observation(Posterior(lklhd_all, prior), rand(prior))

end


using ForwardDiff
using FiniteDifferences
@testset "Bayes Non-analytic ForwardDiff" begin
    _,vis, amp, lcamp, cphase = load_data()

    mfd = central_fdm(5,1)
    @testset "DFT" begin
        g = imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256)
        mt = (g=g, alg=DFTAlg())
        lklhd = RadioLikelihood(test_model2, vis; skymeta=mt)

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
        g = imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256)
        mt = (g=g, alg=DFTAlg())
        lklhd = RadioLikelihood(test_model2, vis; skymeta=mt)

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
        g = imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256)
        mt = (g=g, alg=DFTAlg())
        lklhd = RadioLikelihood(test_model2, vis; skymeta=mt)

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

end
