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
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256)
    skym = SkyModel(test_model, test_prior(), g)
    post_cl = VLBIPosterior(skym, lcamp, cphase)
    post    = VLBIPosterior(skym, vis)

    prior = test_prior()

    prior_sample(post)
    tpostf = asflat(post)
    tpostc = ascube(post)
    show(post)
    show(tpostf)
    show(tpostc)
    # tpost = Comrade.flatten(post)

    ndim = dimension(post)
    @inferred logdensityof(post, prior_sample(post))
    @inferred logdensityof(tpostf, prior_sample(tpostf))
    @inferred logdensityof(tpostc, prior_sample(tpostc))

    x = prior_sample(tpostf)
    @test logdensityof(tpostf, x) == tpostf(x)
    x = prior_sample(tpostc)
    @test logdensityof(tpostc, x) == tpostc(x)

    @test dataproducts(post) == (post.data)


    x = prior_sample(post)
    @test skymodel(post, x) == test_model(x.sky, nothing)
    @test Comrade.idealvisibilities(skymodel(post), x) == visibilitymap(test_model(x.sky, nothing), post.skymodel.grid)
    @test Comrade.forward_model(post, x) == visibilitymap(test_model(x.sky, nothing), post.skymodel.grid)


    @test dimension(post) == length(post.prior)
    @test dimension(tpostf) == length(prior_sample(tpostf))
    @test dimension(tpostc) == length(prior_sample(tpostc))

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
    @test isapprox(xopt.sky.f1/xopt.sky.f2, 2.0, atol=1e-3)
    @test isapprox(xopt.sky.σ1*2*sqrt(2*log(2)), μas2rad(40.0), rtol=1e-3)
    @test isapprox(xopt.sky.σ1*xopt.sky.τ1*2*sqrt(2*log(2)), μas2rad(20.0), rtol=1e-3)
    @test isapprox(xopt.sky.ξ1, π/3, atol=1e-3)
    @test isapprox(xopt.sky.σ2*2*sqrt(2*log(2)), μas2rad(20.0), atol=1e-3)
    @test isapprox(xopt.sky.σ2*xopt.sky.τ2*2*sqrt(2*log(2)), μas2rad(10.0), rtol=1e-3)
    @test isapprox(xopt.sky.ξ2, π/6, atol=1e-3)
    @test isapprox(xopt.sky.x, μas2rad(30.0), rtol=1e-3)
    @test isapprox(xopt.sky.y, μas2rad(30.0), rtol=1e-3)

    show(IOBuffer(), MIME"text/plain"(), post)
    show(IOBuffer(), MIME"text/plain"(), tpostf)

    mopt = test_model(xopt.sky, nothing)
    @test mopt == Comrade.skymodel(post, xopt)

    @testset "Plot residuals" begin
        post = VLBIPosterior(skym, vis)
        residual(post, xopt)
        post = VLBIPosterior(skym, amp)
        residual(post, xopt)
        post = VLBIPosterior(skym, lcamp)
        residual(post, xopt)
        post = VLBIPosterior(skym, cphase)
        residual(post, xopt)
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

@testset "Polarized" begin
    _,vis, amp, lcamp, cphase, coh = load_data()

    R = JonesR()
    intm = InstrumentModel(R, NamedTuple())
    skym = SkyModel(test_skymodel_polarized, test_prior(), imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256))
    post = VLBIPosterior(skym, intm, coh)
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
