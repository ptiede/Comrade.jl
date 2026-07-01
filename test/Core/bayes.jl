using Distributions
using Optimization
using OptimizationLBFGSB
using StatsBase
using Plots
using LogDensityProblems
using LogDensityProblemsAD
using Pyehtim
using Enzyme


@testset "bayes" begin
    _, vis, amp, lcamp, cphase = load_data()
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256)
    skym = SkyModel(test_model, test_prior(), g)
    post_cl = VLBIPosterior(skym, lcamp, cphase; admode = set_runtime_activity(Enzyme.Reverse))
    post = VLBIPosterior(skym, vis)

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
    @test last(Comrade.idealmaps(Comrade.VisData(), skymodel(post), x)) == visibilitymap(test_model(x.sky, nothing), post.skymodel.grid)
    @test last(Comrade.forward_model(post, x)) == visibilitymap(test_model(x.sky, nothing), post.skymodel.grid)


    @test dimension(post) == length(post.prior)
    @test dimension(tpostf) == length(prior_sample(tpostf))
    @test dimension(tpostc) == length(prior_sample(tpostc))

    show(IOBuffer(), MIME"text/plain"(), post)
    show(IOBuffer(), MIME"text/plain"(), tpostf)


    f = OptimizationFunction(tpostf)
    x0 = transform(
        tpostf, [
            0.1,
            0.4,
            0.5,
            0.9,
            0.3,
            0.1,
            0.4,
            0.7,
            0.8,
            0.8,
        ]
    )
    xopt, sol = comrade_opt(post, LBFGSB(); initial_params = x0, maxiters = 10_000)

    @test isapprox(xopt.sky.f1 / xopt.sky.f2, 2.0, atol = 1.0e-3)
    @test isapprox(xopt.sky.σ1 * 2 * sqrt(2 * log(2)), μas2rad(40.0), rtol = 1.0e-3)
    @test isapprox(xopt.sky.σ1 * xopt.sky.τ1 * 2 * sqrt(2 * log(2)), μas2rad(20.0), rtol = 1.0e-3)
    @test isapprox(xopt.sky.ξ1, π / 3, atol = 1.0e-3)
    @test isapprox(xopt.sky.σ2 * 2 * sqrt(2 * log(2)), μas2rad(20.0), atol = 1.0e-3)
    @test isapprox(xopt.sky.σ2 * xopt.sky.τ2 * 2 * sqrt(2 * log(2)), μas2rad(10.0), rtol = 1.0e-3)
    @test isapprox(xopt.sky.ξ2, π / 6, atol = 1.0e-3)
    @test isapprox(xopt.sky.x, μas2rad(30.0), rtol = 1.0e-3)
    @test isapprox(xopt.sky.y, μas2rad(30.0), rtol = 1.0e-3)


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
        @test LogDensityProblems.logdensity(post, x0) ≈ logdensityof(post, x0)
        @test LogDensityProblems.logdensity(tpostf, f0) ≈ logdensityof(tpostf, f0)
        @test LogDensityProblems.logdensity(tpostc, c0) ≈ logdensityof(tpostc, c0)

        @test LogDensityProblems.dimension(tpostf) == length(f0)
        @test LogDensityProblems.dimension(tpostc) == length(c0)

        @test LogDensityProblems.capabilities(typeof(post)) === LogDensityProblems.LogDensityOrder{0}()
        @test LogDensityProblems.capabilities(typeof(tpostf)) === LogDensityProblems.LogDensityOrder{1}()
        @test LogDensityProblems.capabilities(typeof(tpostc)) === LogDensityProblems.LogDensityOrder{1}()
    end

    @testset "corr image prior" begin
        cprior1 = corr_image_prior(g, 10.0; base = EMRF, order = 1)
        cprior2 = corr_image_prior(g, 10.0; base = EMRF, order = 2)

        @test cprior1 isa VLBIImagePriors.HierarchicalPrior
        @test cprior2 isa VLBIImagePriors.HierarchicalPrior

        bs = beamsize(vis)
        @test corr_image_prior(g, bs).hyperprior == corr_image_prior(g, vis).hyperprior
    end


end

using Enzyme
using FiniteDifferences
@testset "Polarized" begin
    _, vis, amp, lcamp, cphase, coh = load_data()

    R = JonesR()
    intm = InstrumentModel(R)
    skym = SkyModel(test_skymodel_polarized, test_prior(), imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256); metadata = (; lp = 0.1))
    post = VLBIPosterior(skym, intm, coh)

    tpost = asflat(post)

    x = prior_sample(tpost)
    gz, = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), Const(tpost), x)
    mfd = central_fdm(5, 1)
    gfd, = FiniteDifferences.grad(mfd, tpost, x)
    @test gz ≈ gfd

    R = JonesR()
    Gp = JonesG(x -> (exp(complex(x.lg, x.gp)), exp(complex(x.lg, x.gp))))
    J = JonesSandwich(Gp, R)
    pr = (
        lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 1.0))),
        gp = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 1.0)), refant = SEFDReference(0.0), phase = true),
    )
    intm_coh = InstrumentModel(J, pr)
    skym = SkyModel(test_skymodel_polarized, test_prior(), imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256); metadata = (; lp = 0.1))
    post = VLBIPosterior(skym, intm_coh, coh)
    tpost = asflat(post)
    x = prior_sample(tpost)
    fj = instrumentmodel(post, prior_sample(post))
    residual(post, Comrade.transform(tpost, x))
    gz, = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), Const(tpost), x)
    mfd = central_fdm(5, 1)
    gfd, = FiniteDifferences.grad(mfd, tpost, x)
    @test gz ≈ gfd

end

@testset "simulate_obs" begin
    _, vis, amp, lcamp, cphase, coh = load_data()
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256)
    skym = SkyModel(test_model, test_prior(), g)

    post_amp = VLBIPosterior(skym, amp)
    post_cp = VLBIPosterior(skym, cphase)
    post_lc = VLBIPosterior(skym, lcamp)
    post_vis = VLBIPosterior(skym, vis)

    G = SingleStokesGain(x -> exp(complex(x.lg, x.gp)))
    intm = InstrumentModel(
        G, (
            lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 1.0))),
            gp = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 1.0)), refant = SEFDReference(0.0), phase = true),
        )
    )
    post_gvis = VLBIPosterior(skym, vis)

    fwhmfac = 2 * sqrt(2 * log(2))
    x0 = (;
        sky = (
            f1 = 1.0, σ1 = μas2rad(40.0) / fwhmfac, τ1 = 0.5, ξ1 = π / 3,
            f2 = 0.5, σ2 = μas2rad(20.0) / fwhmfac, τ2 = 0.5, ξ2 = π / 6,
            x = μas2rad(30.0), y = μas2rad(30.0),
        ),
    )
    function test_simobs(post, x, int = nothing)
        obs = simulate_observation(post, x)[begin]
        @test length(obs) == length(post.data[begin])
        obs_nn = simulate_observation(post, x, add_thermal_noise = false)[begin]
        @test Comrade.measurement(obs_nn) == Comrade.likelihood(post.lklhds[1], last(Comrade.forward_model(post, x))).μ

        if isnothing(int)
            postsim = VLBIPosterior(skym, obs)
            postsim_nn = VLBIPosterior(skym, obs_nn)
        else
            postsim = VLBIPosterior(skym, int, obs)
            postsim_nn = VLBIPosterior(skym, int, obs_nn)
        end


        c2 = chi2(postsim, x; reduce = true)
        c2nn = chi2(postsim_nn, x; reduce = true)
        @test all(x -> reduce(&, x .< 1.25), c2)
        @test all(x -> reduce(&, x .≈ 0), c2nn)


    end

    test_simobs(post_amp, x0)
    test_simobs(post_cp, x0)
    test_simobs(post_lc, x0)
    test_simobs(post_vis, x0)
    test_simobs(post_gvis, x0)

    post_all = VLBIPosterior(skym, vis, amp, lcamp, cphase)
    simulate_observation(post_all, x0)


    R = JonesR()
    Gp = JonesG(x -> (exp(x.lg + 1im * x.gp), exp(x.lg + 1im * x.gp)))
    J = JonesSandwich(Gp, R)
    intm_coh = InstrumentModel(J, intm.prior)
    skym = SkyModel(test_skymodel_polarized, test_prior(), imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256); metadata = (; lp = 0.1))
    post_coh = VLBIPosterior(skym, intm_coh, coh)
    test_simobs(post_coh, prior_sample(post_coh), intm_coh)
end


@testset "Bayes Non-analytic AD" begin
    _, vis, amp, lcamp, cphase = load_data()

    mfd = central_fdm(5, 1)

    G = SingleStokesGain(x -> exp(x.lg + 1im .* x.gp))
    intm = InstrumentModel(
        G, (
            lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 1.0))),
            gp = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 1.0)), refant = SEFDReference(0.0), phase = true),
        )
    )

    function test_nonanalytic(intm, algorithm, vis)
        skym = SkyModel(test_model, test_prior(), imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256))
        post = VLBIPosterior(skym, intm, vis)

        tpostf = asflat(post)
        x0 = prior_sample(tpostf)

        @inferred logdensityof(tpostf, x0)
        gz, = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), Const(tpostf), x0)
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

@testset "FixedSkyModel" begin
    _, vis, amp, lcamp, cphase = load_data()

    G = SingleStokesGain(x -> exp(x.lg + 1im .* x.gp))
    intm = InstrumentModel(
        G, (
            lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 1.0))),
            gp = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 1.0)), refant = SEFDReference(0.0), phase = true),
        )
    )

    f = test_model
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256)
    skym = SkyModel(f, test_prior(), g)

    x = rand(Comrade.NamedDist(test_prior()))
    m = Comrade.skymodel(skym, x)
    skyf = FixedSkyModel(m, g)

    @testset "With Instrument Model" begin
        post = VLBIPosterior(skyf, intm, vis)
        prior_sample(post)
        tpostf = asflat(post)
        tpostc = ascube(post)

        xf = prior_sample(tpostf)
        xc = prior_sample(tpostc)

        logdensityof(tpostf, xf)
        logdensityof(tpostc, xc)
    end

    @testset "No Instrument Model" begin
        post = VLBIPosterior(skyf, lcamp, cphase)
        prior_sample(post) == (;)
        tpostf = asflat(post)
        tpostc = ascube(post)

        xf = prior_sample(tpostf)
        xc = prior_sample(tpostc)

        l1 = logdensityof(tpostf, xf)
        l2 = logdensityof(tpostc, xc)

        @test l1 ≈ l2
    end


end

@testset "Visibility regrouping" begin
    _, vis, amp, lcamp, cphase = load_data()
    gxy = imagepixels(μas2rad(150.0), μas2rad(150.0), 32, 32)
    g3 = RectiGrid((; X = gxy.X, Y = gxy.Y, Fr = [220.0e9, 230.0e9, 240.0e9]))
    g4 = RectiGrid((; X = gxy.X, Y = gxy.Y, Ti = [0.0, 1.0], Fr = [220.0e9, 230.0e9]))
    skym2 = SkyModel(test_model, test_prior(), gxy)
    skym3 = SkyModel(test_model, test_prior(), g3)
    skym4 = SkyModel(test_model, test_prior(), g4)

    # grouping axes come from the image grid's non-spatial dims, in regular Julia (fastest-first) order
    @test Comrade._grouping_axes(skym2) == ()
    @test Comrade._grouping_axes(skym3) == (:Fr,)
    @test Comrade._grouping_axes(skym4) == (:Ti, :Fr)

    # fake a multifrequency dataset by cycling shuffled frequencies onto the real visibility points
    vismf = deepcopy(vis)
    n = length(vismf)
    frs = [220.0e9, 230.0e9, 240.0e9]
    arrayconfig(vismf)[:Fr] .= [frs[(i % 3) + 1] for i in 0:(n - 1)]
    @test !issorted(arrayconfig(vismf)[:Fr])

    dp = Comrade._regroup_dataproducts(skym3, (vismf,))
    @test issorted(arrayconfig(dp[1])[:Fr])
    # measurement / config are permuted by the same Fr-sorting permutation (kept aligned)
    p = sortperm(collect(arrayconfig(vismf)[:Fr]); alg = Base.Sort.DEFAULT_STABLE)
    @test measurement(dp[1]) == measurement(vismf)[p]
    @test arrayconfig(dp[1])[:U] == arrayconfig(vismf)[:U][p]

    # no-op for a plain X/Y image (nothing to group) -> same object returned
    @test Comrade._regroup_dataproducts(skym2, (vismf,)) === (vismf,)

    # closure products short-circuit: a 3D grid would otherwise regroup the vis, but the presence of
    # a closure product leaves everything untouched
    mixed = (vismf, cphase)
    @test Comrade._regroup_dataproducts(skym3, mixed) === mixed

    # the constructor wires it in: the posterior's stored data is Fr-contiguous
    post3 = VLBIPosterior(skym3, vismf; admode = nothing)
    @test issorted(arrayconfig(dataproducts(post3)[begin])[:Fr])

    # non-ascending grid: the data must follow the grid's STORED plane order (matching the NUFFT's
    # `DimPoints` enumeration), NOT ascending value — otherwise the data blocks would not co-locate
    # with their image planes under sharding. Here the grid is descending, so the regrouped data is too.
    gdesc = RectiGrid((; X = gxy.X, Y = gxy.Y, Fr = [240.0e9, 230.0e9, 220.0e9]))
    skymdesc = SkyModel(test_model, test_prior(), gdesc)
    dpd = Comrade._regroup_dataproducts(skymdesc, (vismf,))
    frd = collect(arrayconfig(dpd[1])[:Fr])
    @test issorted(frd; rev = true)                 # blocks in grid order [240, 230, 220] => descending
    @test frd != sort(frd)                          # an ascending sort would NOT match the grid here
    # explicit: measurement follows the position-within-grid ranking (the fix), not the value sort
    posd = Dict(v => i for (i, v) in enumerate([240.0e9, 230.0e9, 220.0e9]))
    pd = sortperm([posd[v] for v in collect(arrayconfig(vismf)[:Fr])]; alg = Base.Sort.DEFAULT_STABLE)
    @test measurement(dpd[1]) == measurement(vismf)[pd]

    # floating-point round-off in the data frequencies still matches the grid plane (isapprox lookup),
    # rather than routing every point to the absent-coordinate fallback and silently skipping regrouping
    visfp = deepcopy(vismf)
    arrayconfig(visfp)[:Fr] .= collect(arrayconfig(vismf)[:Fr]) .* (1 + 1.0e-12)
    dpfp = Comrade._regroup_dataproducts(skym3, (visfp,))
    @test issorted(collect(arrayconfig(dpfp[1])[:Fr]))
    @test measurement(dpfp[1]) == measurement(visfp)[p]

    # differing-length non-closure products cannot share one permutation -> clear error, not a silent
    # row-drop or BoundsError
    @test_throws ArgumentError Comrade._regroup_dataproducts(skym3, (vismf, vismf[1:(n - 2)]))
end
