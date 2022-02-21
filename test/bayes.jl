using Distributions
using AdvancedHMC
using Pathfinder
using GalacticOptim
using NestedSamplers
using BlackBoxOptim
using StatsBase

@testset "bayes" begin
    load_ehtim()
    obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "test_data.uvfits"))
    obs.add_scans()
    obsavg = obs.avg_coherent(0.0, scan_avg=true)

    m = ehtim.model.Model()
    m = m.add_gauss(1.0, μas2rad(40.0), μas2rad(20.0), π/3, 0.0, 0.0)
    obsm = m.observe_same_nonoise(obsavg)

    amp = extract_amp(obsm; debias=true)
    cphase = extract_cphase(obsm, count="min")

    lklhd = RadioLikelihood(amp, cphase)

    function model(θ)
        θ.f*rotated(stretched(Gaussian(), θ.σ*θ.τ, θ.σ), θ.ξ)
    end

    prior = (f=Uniform(0.8, 1.2),
             σ = Uniform(μas2rad(10.0), μas2rad(40.0)),
             τ = Uniform(0.01, 0.99),
             ξ = Uniform(-π/2, π/2)
            )

    post = Posterior(lklhd, prior, model)
    tpostf = asflat(post)
    tpostc = ascube(post)
    tpost = Comrade.flatten(post)

    @test dimension(post) == dimension(tpost)
    @test dimension(post) == dimension(tpostf)
    @test dimension(post) == dimension(tpostc)

    f,tr = OptimizationFunction(post, GalacticOptim.AutoForwardDiff{4}())
    prob = GalacticOptim.OptimizationProblem(f, rand(post.prior), nothing, lb=fill(-10.0, 4), ub=fill(10.0,4))
    sol, xopt = solve(prob, BBO_adaptive_de_rand_1_bin(), tr; maxiters=50_000)

    xtest = transform(tr, sol.u)
    @test collect(values(xtest)) ≈ collect(values(xopt))
    @test isapprox(xopt.f, 1.0, atol=1e-3)
    @test isapprox(xopt.σ*2*sqrt(2*log(2)), μas2rad(40.0), atol=1e-3)
    @test isapprox(xopt.σ*xopt.τ*2*sqrt(2*log(2)), μas2rad(20.0), atol=1e-3)
    @test isapprox(-xopt.ξ, π/3, atol=1e-3)

    q1, ϕ1, logqϕ1 = multipathfinder(post, 10; init_params=2)
    q2, ϕ2, logqϕ2 = multipathfinder(post, 10, init_params=rand(post.prior, 2))

    nchain, nstats = sample(post, Nested(dimension(post), 2_000); dlogz=1.0)
    @test isapprox(collect(values(nchain[end])), collect(values(xopt)), rtol=1e-2)
    echain = Comrade.TupleVector(sample(nchain, Weights(nstats.weights), 100_000))
    mn = Comrade.rmap(mean, echain)[1]
    sn = Comrade.rmap(std, echain)[1]

    hchain, hstats = sample(post, HMC(metric=DiagEuclideanMetric(dimension(post))),
                            10_000; nadapts=2000, init_params=nchain[end], progress=false)
    mh = Comrade.rmap(mean, hchain)[1]
    sh = Comrade.rmap(std, hchain)[1]

    @test isapprox(collect(values(mh)), collect(values(mn)), rtol=1e-3)
    @test isapprox(collect(values(sh)), collect(values(sn)), rtol=1e-1)

end
