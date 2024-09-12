using Distributions
using Optimization
using OptimizationGCMAES
using StatsBase
using Plots
using LogDensityProblems
using LogDensityProblemsAD
using Pyehtim
using Enzyme



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

    show(IOBuffer(), MIME"text/plain"(), post)
    show(IOBuffer(), MIME"text/plain"(), tpostf)


    f = OptimizationFunction(tpostf, Optimization.AutoEnzyme(;mode=Enzyme.Reverse))
    x0 = transform(tpostf, [ 0.1,
           0.4,
           0.5,
           0.9,
           0.3,
           0.1,
           0.4,
           0.7,
           0.8,
           0.8])
    xopt, sol = comrade_opt(post, GCMAESOpt(); initial_params=x0, maxiters=10_000)

    @test isapprox(xopt.sky.f1/xopt.sky.f2, 2.0, atol=1e-3)
    @test isapprox(xopt.sky.σ1*2*sqrt(2*log(2)), μas2rad(40.0), rtol=1e-3)
    @test isapprox(xopt.sky.σ1*xopt.sky.τ1*2*sqrt(2*log(2)), μas2rad(20.0), rtol=1e-3)
    @test isapprox(xopt.sky.ξ1, π/3, atol=1e-3)
    @test isapprox(xopt.sky.σ2*2*sqrt(2*log(2)), μas2rad(20.0), atol=1e-3)
    @test isapprox(xopt.sky.σ2*xopt.sky.τ2*2*sqrt(2*log(2)), μas2rad(10.0), rtol=1e-3)
    @test isapprox(xopt.sky.ξ2, π/6, atol=1e-3)
    @test isapprox(xopt.sky.x, μas2rad(30.0), rtol=1e-3)
    @test isapprox(xopt.sky.y, μas2rad(30.0), rtol=1e-3)


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

using Enzyme
using FiniteDifferences
@testset "Polarized" begin
    _,vis, amp, lcamp, cphase, coh = load_data()

    R = JonesR()
    intm = InstrumentModel(R)
    skym = SkyModel(test_skymodel_polarized, test_prior(), imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256); metadata=(;lp=0.1))
    post = VLBIPosterior(skym, intm, coh)

    tpost = asflat(post)

    x = prior_sample(tpost)
    gz = Enzyme.gradient(Enzyme.Reverse, Const(tpost), x)
    mfd = central_fdm(5,1)
    gfd, = FiniteDifferences.grad(mfd, tpost, x)
    @test gz ≈ gfd

    R = JonesR()
    Gp = JonesG(x->(exp(x.lg + 1im*x.gp), exp(x.lg + 1im*x.gp)))
    J = JonesSandwich(Gp, R)
    pr = (lg = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 1.0))),
          gp = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 1.0)), refant=SEFDReference(0.0), phase=true))
    intm_coh = InstrumentModel(J, pr)
    skym = SkyModel(test_skymodel_polarized, test_prior(), imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256); metadata=(;lp=0.1))
    post = VLBIPosterior(skym, intm_coh, coh)
    tpost = asflat(post)
    x = prior_sample(tpost)
    residual(post, Comrade.transform(tpost, x))
    gz = Enzyme.gradient(Enzyme.Reverse, Const(tpost), x)
    mfd = central_fdm(5,1)
    gfd, = FiniteDifferences.grad(mfd, tpost, x)
    @test gz ≈ gfd

end

@testset "simulate_obs" begin
    _,vis, amp, lcamp, cphase, coh = load_data()
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256)
    skym = SkyModel(test_model, test_prior(), g)

    post_amp = VLBIPosterior(skym, amp)
    post_cp = VLBIPosterior(skym, cphase)
    post_lc = VLBIPosterior(skym, lcamp)
    post_vis = VLBIPosterior(skym, vis)

    G = SingleStokesGain(x->exp(x.lg + 1im.*x.gp))
    intm = InstrumentModel(G, (lg = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 1.0))),
                               gp = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 1.0)), refant=SEFDReference(0.0), phase=true)))
    post_gvis = VLBIPosterior(skym, vis)

    function test_simobs(post, x)
        obs = simulate_observation(post, x)[begin]
        @test length(obs) == length(post.data[begin])
        obs_nn = simulate_observation(post, x, add_thermal_noise=false)[begin]
        @test Comrade.measurement(obs_nn) == Comrade.likelihood(post.lklhds[1], Comrade.forward_model(post, x)).μ
    end

    test_simobs(post_amp, prior_sample(post_amp))
    test_simobs(post_cp,  prior_sample(post_cp))
    test_simobs(post_lc,  prior_sample(post_lc))
    test_simobs(post_vis, prior_sample(post_vis))
    test_simobs(post_gvis, prior_sample(post_gvis))

    post_all = VLBIPosterior(skym, vis, amp, lcamp, cphase)
    simulate_observation(post_all, prior_sample(post_all))


    R = JonesR()
    Gp = JonesG(x->(exp(x.lg + 1im*x.gp), exp(x.lg + 1im*x.gp)))
    J = JonesSandwich(Gp, R)
    intm_coh = InstrumentModel(J, intm.prior)
    skym = SkyModel(test_skymodel_polarized, test_prior(), imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256); metadata=(;lp=0.1))
    post_coh = VLBIPosterior(skym, intm_coh, coh)
    test_simobs(post_coh, prior_sample(post_coh))
end


@testset "Bayes Non-analytic AD" begin
    _,vis, amp, lcamp, cphase = load_data()

    mfd = central_fdm(5,1)

    G = SingleStokesGain(x->exp(x.lg + 1im.*x.gp))
    intm = InstrumentModel(G, (lg = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 1.0))),
                               gp = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 1.0)), refant=SEFDReference(0.0), phase=true)))

    function test_nonanalytic(intm, algorithm, vis)
        skym = SkyModel(test_model, test_prior(), imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256))
        post = VLBIPosterior(skym, intm, vis)

        tpostf = asflat(post)
        x0 = prior_sample(tpostf)

        @inferred logdensityof(tpostf, x0)
        gz = Enzyme.gradient(Enzyme.Reverse, Const(tpostf), x0)
        gn, = FiniteDifferences.grad(mfd, tpostf, x0)
        @test gz ≈ gn
    end
    @testset "DFT" begin
        test_nonanalytic(intm, DFTAlg(), vis)
    end

    @testset "NFFT" begin
        test_nonanalytic(intm, NFFTAlg(), vis)
    end

end
