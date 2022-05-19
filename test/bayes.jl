using Distributions
using AdvancedHMC
using Pathfinder
using GalacticOptim
using NestedSamplers
using BlackBoxOptim
using StatsBase
using Plots
using Optim

@testset "bayes" begin
    load_ehtim()
    obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "test_data.uvfits"))
    obs.add_scans()
    obsavg = obs.avg_coherent(0.0, scan_avg=true)

    m = ehtim.model.Model()
    m = m.add_gauss(1.0, μas2rad(40.0), μas2rad(20.0), π/3, 0.0, 0.0)
    m = m.add_gauss(0.5, μas2rad(20.0), μas2rad(10.0), π/6, μas2rad(30.0), μas2rad(30.0))
    obsm = m.observe_same_nonoise(obsavg)

    amp = extract_lcamp(obsm)
    cphase = extract_cphase(obsm, cut_trivial=true)

    lklhd = RadioLikelihood(amp, cphase)

    function model(θ)
        m1 = θ.f1*rotated(stretched(Gaussian(), θ.σ1*θ.τ1, θ.σ1), θ.ξ1)
        m2 = θ.f2*rotated(stretched(Gaussian(), θ.σ2*θ.τ2, θ.σ2), θ.ξ2)
        return m1 + shifted(m2, θ.x, θ.y)
    end

    prior = (f1=Uniform(0.8, 1.2),
             σ1 = Uniform(μas2rad(1.0), μas2rad(40.0)),
             τ1 = Uniform(0.35, 0.65),
             ξ1 = Uniform(-π/2, π/2),
             f2 = Uniform(0.3, 0.7),
             σ2 = Uniform(μas2rad(1.0), μas2rad(40.0)),
             τ2 = Uniform(0.35, 0.65),
             ξ2 = Uniform(-π/2, π/2),
             x = Uniform(-μas2rad(40.0), μas2rad(40.0)),
             y = Uniform(-μas2rad(40.0), μas2rad(40.0))
            )

    post = Posterior(lklhd, prior, model)
    tpostf = asflat(post)
    tpostc = ascube(post)
    tpost = Comrade.flatten(post)

    @test dimension(post) == dimension(tpost)
    @test dimension(post) == dimension(tpostf)
    @test dimension(post) == dimension(tpostc)

    f,tr = OptimizationFunction(post, GalacticOptim.AutoForwardDiff{4}())
    prob = GalacticOptim.OptimizationProblem(f, rand(post.prior), nothing, lb=fill(-10.0, 10), ub=fill(10.0,10))
    sol, xopt = solve(prob, BBO_adaptive_de_rand_1_bin(), tr; maxiters=500_000)

    xtest = transform(tr, sol.u)
    @test collect(values(xtest)) ≈ collect(values(xopt))
    @test isapprox(xopt.f1-xopt.f2, 0.5, atol=1e-3)
    @test isapprox(xopt.σ1*2*sqrt(2*log(2)), μas2rad(40.0), rtol=1e-3)
    @test isapprox(xopt.σ1*xopt.τ1*2*sqrt(2*log(2)), μas2rad(20.0), rtol=1e-3)
    @test isapprox(-xopt.ξ1, π/3, atol=1e-3)
    @test isapprox(xopt.σ2*2*sqrt(2*log(2)), μas2rad(20.0), atol=1e-3)
    @test isapprox(xopt.σ2*xopt.τ2*2*sqrt(2*log(2)), μas2rad(10.0), rtol=1e-3)
    @test isapprox(-xopt.ξ2, π/6, atol=1e-3)
    @test isapprox(xopt.x, μas2rad(30.0), rtol=1e-3)
    @test isapprox(xopt.y, μas2rad(30.0), rtol=1e-3)


    nchain, nstats = sample(post, Nested(dimension(post), 2_000); dlogz=0.01, progress=true)
    @test isapprox(collect(values(nchain[end])), collect(values(xopt)), rtol=1e-2)
    echain = Comrade.TupleVector(sample(nchain, Weights(nstats.weights), 500_000))
    mn = Comrade.rmap(mean, echain)[1]
    sn = Comrade.rmap(std, echain)[1]

    res = pathfinder(post, 10; init_params=xopt)
    q2, ϕ2, logqϕ2 = multipathfinder(post, 10, init_params=nchain[end:-1:end-5])


    hchain, hstats, res = sample(post, HMC(metric=DiagEuclideanMetric(dimension(post))),
                            10_000; nadapts=2000, init_params=nchain[end], progress=true)
    mh = Comrade.rmap(mean, hchain)[1]
    sh = Comrade.rmap(std, hchain)[1]

    @test isapprox(collect(values(mh)), collect(values(mn)), rtol=1e-3)
    @test isapprox(collect(values(sh)), collect(values(sn)), rtol=1e-1)

end
