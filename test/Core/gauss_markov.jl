using Distributions
using VLBIFiles
using LinearAlgebra
using Random
using Statistics
using FiniteDifferences
using Enzyme
using LogDensityProblems
import TransformVariables as TV

_gm_test_data() = load_test_table()

function _gm_test_sky()
    tm(θ, _, _) = θ.f1 * stretched(Gaussian(), θ.σ1, θ.σ1)
    prior = (f1 = VLBIUniform(0.5, 1.5), σ1 = VLBIUniform(μas2rad(1.0), μas2rad(100.0)))
    g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 64, 64)
    return SkyModel(tm, prior, g)
end

@testset "GaussMarkov instrument priors" begin
    rng = Random.Xoshiro(42)


    dvis = _gm_test_data()
    skym = _gm_test_sky()

    @instrument function gmint_fitted()
        return @jones begin
            lg ~ ArrayPrior(GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = Exponential(0.1), τ = Exponential(2.0))))
            gp ~ ArrayPrior(GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = 2.0, τ = InverseGamma(3.0, 6.0))); refant = SensitiveReference(0.0))
            return Jones(exp(complex(lg, gp)))
        end
    end

    @testset "posterior with fitted hyperparameters" begin
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_fitted(), admode = set_runtime_activity(Enzyme.Reverse))
        s = prior_sample(rng, post)
        @test s.instrument.lg isa NamedTuple{(:params, :hyperparams)}
        @test s.instrument.lg.params isa Comrade.StationArray
        # every station gets its own hyperparameters, nested under the station name
        @test keys(s.instrument.lg.hyperparams) == Tuple(Comrade.stations(dvis))
        @test all(h -> keys(h) == (:σ, :τ), values(s.instrument.lg.hyperparams))
        # gp fits only τ; σ is fixed
        @test all(h -> keys(h) == (:τ,), values(s.instrument.gp.hyperparams))
        @test isfinite(logdensityof(post, s))
        @inferred logdensityof(post, s)

        tpostf = asflat(post)
        xf = prior_sample(rng, tpostf)
        @test isfinite(logdensityof(tpostf, xf))
        @inferred logdensityof(tpostf, xf)

        # transform roundtrip (all component transforms here are bijective)
        y = Comrade.transform(tpostf, xf)
        @test Comrade.inverse(tpostf, y) ≈ xf
        # the reference-fixed phases are present in the transformed params
        @test any(iszero, parent(y.instrument.gp.params))

        # whitened (default) priors transport exactly to the Std spaces
        tpostc = ascube(post)
        xc = prior_sample(rng, tpostc)
        @test isfinite(logdensityof(tpostc, xc))
        yc = Comrade.transform(tpostc, xc)
        @test Comrade.inverse(tpostc, yc) ≈ xc

        # caltable strips the hyperparameters
        ct = caltable(s.instrument.lg)
        @test ct isa VLBIInstrumentModels.CalTable

        # residuals/instrument forward model work on hierarchical samples
        @test Comrade.instrumentmodel(post, s) isa Comrade.DD.AbstractDimArray

        @testset "Enzyme gradient" begin
            f = let tp = tpostf
                x -> logdensityof(tp, x)
            end
            gz, = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), Const(f), xf)
            gfd, = grad(central_fdm(5, 1), f, xf)
            @test gz ≈ gfd rtol = 1.0e-5
            ld, gld = LogDensityProblems.logdensity_and_gradient(tpostf, xf)
            @test ld ≈ f(xf)
            @test gld ≈ gfd rtol = 1.0e-5
        end
    end

    @testset "hierarchical PosteriorSamples postprocessing" begin
        # `(params, hyperparams)` samples pass through PosteriorSamples, `rmap`, `caltable`,
        # and `stationparams`
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_fitted())
        samples = [prior_sample(rng, post) for _ in 1:20]
        chain = Comrade.PosteriorSamples(samples, nothing; metadata = Dict(:sampler => :prior))
        @test length(Comrade.postsamples(chain)) == 20

        mchain = Comrade.rmap(mean, chain)
        schain = Comrade.rmap(std, chain)
        # the hierarchical instrument terms survive reduction as (params, hyperparams) tuples
        @test mchain.instrument.lg.params isa Comrade.StationArray
        @test keys(mchain.instrument.lg.hyperparams) == Tuple(Comrade.stations(dvis))
        @test Comrade.stationparams(mchain.instrument.lg) === mchain.instrument.lg.params

        # caltable strips the hyperparameters from the reduced chain
        @test caltable(mchain.instrument.lg) isa VLBIInstrumentModels.CalTable
        @test caltable(abs.(mchain.instrument.lg.params)) isa VLBIInstrumentModels.CalTable

        # instrumentmodel accepts the reduced hierarchical chain (stationparams strips hyperparams)
        @test instrumentmodel(post, mchain) isa Comrade.DD.AbstractDimArray
        # combining mean and std: stationparams unwraps both before the elementwise op
        combined = map(
            (x, y) -> Comrade.stationparams(x) .+ Comrade.stationparams(y),
            mchain.instrument, schain.instrument,
        )
        @test instrumentmodel(post, (; instrument = combined)) isa Comrade.DD.AbstractDimArray
    end

    @testset "empty-group PosteriorSamples (zero-field StructArray guard)" begin
        # an empty NamedTuple leaf (a group with no fitted parameters)
        raw = [(a = i * 1.0, empty = NamedTuple()) for i in 1:5]
        ps = Comrade.PosteriorSamples(raw, nothing)
        @test length(Comrade.postsamples(ps)) == 5
        @test ps[3].a == 3.0
        @test ps[1].empty == NamedTuple()
    end

    @testset "derived element type" begin
        # eltype is derived from the process/override parameters, not pinned to Float64
        gmp(p) = GaussMarkovStationPrior(ScanSeg(); Ti = p)
        @test VLBIInstrumentModels._working_type((a = gmp(Comrade.OrnsteinUhlenbeck(σ = 0.1, τ = 2.0)),)) === Float64
        @test VLBIInstrumentModels._working_type((a = gmp(Comrade.OrnsteinUhlenbeck(σ = 0.1f0, τ = 2.0f0, μ = 0.0f0)),)) === Float32
        # fitted (Distribution) fields contribute their own eltype
        @test VLBIInstrumentModels._working_type((a = gmp(Comrade.OrnsteinUhlenbeck(σ = Exponential(0.1), τ = 2.0)),)) === Float64
        # the init contributes to the working type as well
        @test VLBIInstrumentModels._working_type(
            (a = GaussMarkovStationPrior(ScanSeg(); Ti = Comrade.OrnsteinUhlenbeck(σ = 0.1f0, τ = 2.0f0, μ = 0.0f0), init = FixedInit(0.0)),)
        ) === Float64
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_fitted())
        @test eltype(post.prior.instrument.lg) === Float64
        @test eltype(parent(prior_sample(rng, post).instrument.lg.params)) === Float64
    end

    @testset "fixed hyperparameters degrade to plain prior" begin
        @instrument function gmint_fixed()
            return @jones begin
                lg ~ ArrayPrior(GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = 0.1, τ = 2.0)))
                gp ~ ArrayPrior(IIDStationPrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))); refant = SensitiveReference(0.0))
                return Jones(exp(complex(lg, gp)))
            end
        end
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_fixed())
        s = prior_sample(rng, post)
        @test s.instrument.lg isa Comrade.StationArray
        @test isfinite(logdensityof(post, s))
        tp = asflat(post)
        @test isfinite(logdensityof(tp, prior_sample(rng, tp)))
        @test caltable(s.instrument.lg) isa VLBIInstrumentModels.CalTable
    end

    @testset "station overrides" begin
        # IID override under a GaussMarkov default
        @instrument function gmint_iidover()
            return @jones begin
                lg ~ ArrayPrior(
                    GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = Exponential(0.1), τ = 2.0));
                    LM = IIDStationPrior(ScanSeg(), VLBIGaussian(0.0, 1.0))
                )
                return Jones(exp(lg))
            end
        end
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_iidover())
        s = prior_sample(rng, post)
        # the IID-override station carries no hyperparameters; the GM stations each get their own
        @test :LM ∉ keys(s.instrument.lg.hyperparams)
        @test all(h -> keys(h) == (:σ,), values(s.instrument.lg.hyperparams))
        @test isfinite(logdensityof(post, s))
        tp = asflat(post)
        @test isfinite(logdensityof(tp, prior_sample(rng, tp)))

        # GaussMarkov override with its own fitted hyperparameters, nested under the station
        @instrument function gmint_gmover()
            return @jones begin
                lg ~ ArrayPrior(
                    GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = Exponential(0.1), τ = Exponential(2.0)));
                    LM = GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = Exponential(1.0), τ = 2.0))
                )
                return Jones(exp(lg))
            end
        end
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_gmover())
        s = prior_sample(rng, post)
        @test keys(s.instrument.lg.hyperparams) == Tuple(Comrade.stations(dvis))
        @test keys(s.instrument.lg.hyperparams.AA) == (:σ, :τ)
        # the override station fixes τ, so only σ is fitted there
        @test keys(s.instrument.lg.hyperparams.LM) == (:σ,)
        @test isfinite(logdensityof(post, s))
        tp = asflat(post)
        @test isfinite(logdensityof(tp, prior_sample(rng, tp)))

        # an IID override on the reference station has every point fixed, so its chain term
        # is an empty sum
        @instrument function gmint_iidrefover()
            return @jones begin
                gp ~ ArrayPrior(
                    GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = 2.0, τ = 1.0));
                    AA = IIDStationPrior(ScanSeg(), VLBIGaussian(0.0, 1.0)),
                    refant = SingleReference(:AA, 0.0)
                )
                return Jones(exp(1im * gp))
            end
        end
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_iidrefover())
        s = prior_sample(rng, post)
        @test isfinite(logdensityof(post, s))
        tp = asflat(post)
        @test isfinite(logdensityof(tp, prior_sample(rng, tp)))
    end

    @testset "whitened flat density identity" begin
        # Under whitening the flat prior is standard normal in z: changing only z shifts
        # the flat log density by −Δ‖z‖²/2. This checks the coloring map and its Jacobian.
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_fitted())
        nstation = length(Comrade.stations(dvis))
        for (name, nhp) in ((:lg, 2nstation), (:gp, nstation))  # number of hyperparameter coordinates
            obsprior = getproperty(post.prior.instrument, name)
            node = Comrade.transport_node(obsprior, Comrade.TVFlat())
            dim = TV.dimension(node)
            f = x -> begin
                y, ℓj, _ = TV.transform_with(TV.LogJac(), node, x, 1)
                return logpdf(obsprior, y) + ℓj
            end
            x1 = randn(rng, dim)
            x2 = copy(x1)
            x2[(nhp + 1):end] .= randn(rng, dim - nhp)
            Δz² = sum(abs2, @view(x2[(nhp + 1):end])) - sum(abs2, @view(x1[(nhp + 1):end]))
            @test f(x1) - f(x2) ≈ Δz² / 2
        end
    end

    @testset "cube pushforward matches flat pushforward (fixed hp)" begin
        # For a pure-Markov whitened prior with fixed hyperparameters the cube transport
        # at latents u must coincide with the flat transport at z = Φ⁻¹(u) coordinatewise.
        @instrument function gmint_fixref()
            return @jones begin
                gp ~ ArrayPrior(GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = 2.0, τ = 1.0)); refant = SensitiveReference(0.0))
                return Jones(exp(1im * gp))
            end
        end
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_fixref())
        op = post.prior.instrument.gp
        nf = Comrade.transport_node(op, Comrade.TVFlat())
        nc = Comrade.transport_node(op, Comrade.StdUniform())
        @test TV.dimension(nf) == Comrade.PT.dimension(nc)
        u = rand(rng, TV.dimension(nf))
        z = quantile.(Normal(), u)
        yf = TV.transform(nf, z)
        ycube, _ = Comrade.PT.pfwd_step(nc, u, 1)
        @test parent(yf) ≈ parent(ycube)
        # and both roundtrip
        @test TV.inverse(nf, yf) ≈ z
        xb = Vector{Float64}(undef, Comrade.PT.dimension(nc))
        Comrade.PT.pback_step!(xb, 1, nc, ycube)
        @test xb ≈ u
    end

    @testset "hierarchical cube pushforward moments" begin
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_fitted())
        obsprior = post.prior.instrument.lg
        nc = Comrade.transport_node(obsprior, Comrade.StdUniform())
        N = 5000
        ys = map(1:N) do _
            first(Comrade.PT.pfwd_step(nc, rand(rng, Comrade.PT.dimension(nc)), 1))
        end
        yd = map(_ -> rand(rng, obsprior), 1:N)
        # hyperparameter and gain moments agree between the cube pushforward and rand
        @test mean(y -> y.hyperparams.AA.σ, ys) ≈ mean(y -> y.hyperparams.AA.σ, yd) rtol = 0.1
        @test mean(y -> y.hyperparams.AA.τ, ys) ≈ mean(y -> y.hyperparams.AA.τ, yd) rtol = 0.1
        @test mean(y -> std(parent(y.params)), ys) ≈ mean(y -> std(parent(y.params)), yd) rtol = 0.1
    end

    @testset "bounded IID override under GM default" begin
        # IID overrides use their own transform, which keeps bounded supports bounded
        @instrument function gmint_boundover()
            return @jones begin
                lg ~ ArrayPrior(
                    GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = Exponential(0.1), τ = 2.0));
                    LM = IIDStationPrior(ScanSeg(), VLBIUniform(-1.0, 1.0))
                )
                return Jones(exp(lg))
            end
        end
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_boundover())
        tp = asflat(post)
        for _ in 1:5
            # Far out (3σ in every latent) a hyperparameter can overflow the gains, so the
            # log density is checked at prior-scale draws and the transform far out.
            @test isfinite(logdensityof(tp, randn(rng, dimension(tp))))
            x = 3 .* randn(rng, dimension(tp))
            y = Comrade.transform(tp, x)
            @test all(v -> -1 ≤ v ≤ 1, parent(y.instrument.lg.params)[VLBIInstrumentModels.lookup(post.prior.instrument.lg.stationmap).LM])
            @test Comrade.inverse(tp, y) ≈ x
        end
        tpc = ascube(post)
        @test isfinite(logdensityof(tpc, prior_sample(rng, tpc)))
    end

    @testset "centered option" begin
        @instrument function gmint_centered()
            return @jones begin
                lg ~ ArrayPrior(GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = Exponential(0.1), τ = Exponential(2.0)), centered = true))
                return Jones(exp(lg))
            end
        end
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_centered())
        s = prior_sample(rng, post)
        @test isfinite(logdensityof(post, s))
        tp = asflat(post)
        xf = prior_sample(rng, tp)
        @test isfinite(logdensityof(tp, xf))
        y = Comrade.transform(tp, xf)
        @test Comrade.inverse(tp, y) ≈ xf
        # centered chains need the target density, which Std transports never evaluate
        @test_throws ArgumentError ascube(post)
    end

    @testset "hyperparameter recovery (conditional)" begin
        # Gains simulated from a known OU process: the hyperparameter posterior conditional
        # on the true gains peaks near the truth. A joint MAP over gains and hyperparameters
        # is degenerate for a hierarchical model, so it is not a valid check.
        rrng = Random.Xoshiro(7)
        @instrument function gmint_rec()
            return @jones begin
                lg ~ ArrayPrior(GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = Exponential(0.3), τ = Exponential(2.0))))
                return Jones(exp(lg))
            end
        end
        post0 = VLBIPosterior(dvis; sky = skym, instrument = gmint_rec())
        σtrue, τtrue = 0.15, 1.5
        # every station has its own hyperparameters; use the same truth for all of them
        hpof(σ, τ) = station_tuple(dvis, (σ = σ, τ = τ))
        obsprior = post0.prior.instrument.lg
        xg = Vector{Float64}(undef, length(obsprior.dists))
        VLBIInstrumentModels._rand_chains!(rrng, xg, obsprior.dists, hpof(σtrue, τtrue))
        gtruth = Comrade.StationArray(xg, obsprior.stationmap)
        θtrue = (
            sky = (f1 = 1.0, σ1 = μas2rad(40.0)),
            instrument = (lg = (params = gtruth, hyperparams = hpof(σtrue, τtrue)),),
        )
        @test isfinite(logdensityof(post0, θtrue))
        obs = simulate_observation(rrng, post0, θtrue)
        post = VLBIPosterior(obs[1]; sky = skym, instrument = gmint_rec())
        condld(σ, τ) = logdensityof(
            post,
            (sky = θtrue.sky, instrument = (lg = (params = gtruth, hyperparams = hpof(σ, τ)),))
        )
        σs = range(0.05, 0.4; length = 40)
        τs = range(0.2, 6.0; length = 40)
        L = [condld(σ, τ) for σ in σs, τ in τs]
        imax = argmax(L)
        @test 0.5σtrue < σs[imax[1]] < 2σtrue
        @test 0.3τtrue < τs[imax[2]] < 3τtrue
    end

    @testset "FixedInit chains" begin
        # FixedInit(0.0) fixes every station's chain to zero at its first time, removing the
        # level degeneracy with a separate circular offset term
        @instrument function gmint_anch()
            return @jones begin
                gp0 ~ ArrayPrior(IIDStationPrior(TrackSeg(), DiagonalVonMises(0.0, inv(π^2))))
                dgp ~ ArrayPrior(
                    GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = 1.0, τ = InverseGamma(3.0, 3.0)), init = FixedInit(0.0));
                    refant = SensitiveReference(0.0)
                )
                return Jones(exp(1im * (gp0 + dgp)))
            end
        end
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_anch())
        s = prior_sample(rng, post)
        chinds = VLBIInstrumentModels.lookup(post.prior.instrument.dgp.stationmap)
        # every station's first point is exactly zero, in rand and through the transform
        @test all(inds -> iszero(parent(s.instrument.dgp.params)[first(inds)]), chinds)
        @test isfinite(logdensityof(post, s))
        tp = asflat(post)
        xf = prior_sample(rng, tp)
        @test isfinite(logdensityof(tp, xf))
        y = Comrade.transform(tp, xf)
        @test Comrade.inverse(tp, y) ≈ xf
        @test all(inds -> iszero(parent(y.instrument.dgp.params)[first(inds)]), chinds)

        # ascube works for FixedInit chains (checked without the circular offset term,
        # which refuses Std transports)
        @instrument function gmint_anchonly()
            return @jones begin
                dgp ~ ArrayPrior(
                    GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = 1.0, τ = InverseGamma(3.0, 3.0)), init = FixedInit(0.0));
                    refant = SensitiveReference(0.0)
                )
                return Jones(exp(1im * dgp))
            end
        end
        postc = VLBIPosterior(dvis; sky = skym, instrument = gmint_anchonly())
        tpc = ascube(postc)
        xc = prior_sample(rng, tpc)
        @test isfinite(logdensityof(tpc, xc))
        yc = Comrade.transform(tpc, xc)
        @test Comrade.inverse(tpc, yc) ≈ xc

        # FixedInit + centered + fitted hyperparameters
        @instrument function gmint_anchcent()
            return @jones begin
                dgp ~ ArrayPrior(
                    GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = 1.0, τ = InverseGamma(3.0, 3.0)), init = FixedInit(0.0), centered = true);
                    refant = SensitiveReference(0.0)
                )
                return Jones(exp(1im * dgp))
            end
        end
        postcc = VLBIPosterior(dvis; sky = skym, instrument = gmint_anchcent())
        scc = prior_sample(rng, postcc)
        chindscc = VLBIInstrumentModels.lookup(postcc.prior.instrument.dgp.stationmap)
        @test all(inds -> iszero(parent(scc.instrument.dgp.params)[first(inds)]), chindscc)
        tpcc = asflat(postcc)
        xcc = prior_sample(rng, tpcc)
        @test isfinite(logdensityof(tpcc, xcc))
        ycc = Comrade.transform(tpcc, xcc)
        @test Comrade.inverse(tpcc, ycc) ≈ xcc
        @test all(inds -> iszero(parent(ycc.instrument.dgp.params)[first(inds)]), chindscc)
    end

    @testset "FixedInit value and refant conflict" begin
        # a nonzero init value propagates exactly through rand and the transform
        @instrument function gmint_fixval()
            return @jones begin
                lg ~ ArrayPrior(GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = 0.5, τ = 2.0), init = FixedInit(0.3)))
                return Jones(exp(lg))
            end
        end
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_fixval())
        s = prior_sample(rng, post)
        chinds = VLBIInstrumentModels.lookup(post.prior.instrument.lg.stationmap)
        @test all(inds -> parent(s.instrument.lg)[first(inds)] == 0.3, chinds)
        @test isfinite(logdensityof(post, s))
        tp = asflat(post)
        y = Comrade.transform(tp, prior_sample(rng, tp))
        @test all(inds -> parent(y.instrument.lg)[first(inds)] == 0.3, chinds)

        # a reference value that disagrees with FixedInit at the first point is an error
        @instrument function gmint_conflict()
            return @jones begin
                gp ~ ArrayPrior(GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = 2.0, τ = 1.0), init = FixedInit(0.3)); refant = SensitiveReference(0.0))
                return Jones(exp(1im * gp))
            end
        end
        @test_throws ArgumentError VLBIPosterior(dvis; sky = skym, instrument = gmint_conflict())

        # agreeing values compose: the reference chain's first point is fixed by both
        @instrument function gmint_agree()
            return @jones begin
                gp ~ ArrayPrior(GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = 2.0, τ = 1.0), init = FixedInit(0.0)); refant = SensitiveReference(0.0))
                return Jones(exp(1im * gp))
            end
        end
        post2 = VLBIPosterior(dvis; sky = skym, instrument = gmint_agree())
        @test isfinite(logdensityof(post2, prior_sample(rng, post2)))
    end

    @testset "GaussianInit posterior" begin
        # wide explicit start together with a fitted hyperparameter
        @instrument function gmint_ginit()
            return @jones begin
                lg ~ ArrayPrior(GaussMarkovStationPrior(ScanSeg(); Ti = OrnsteinUhlenbeck(σ = Exponential(0.1), τ = 2.0), init = GaussianInit(0.0, 1.0)))
                return Jones(exp(lg))
            end
        end
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_ginit())
        s = prior_sample(rng, post)
        @test isfinite(logdensityof(post, s))
        tp = asflat(post)
        xf = prior_sample(rng, tp)
        @test isfinite(logdensityof(tp, xf))
        y = Comrade.transform(tp, xf)
        @test Comrade.inverse(tp, y) ≈ xf
        tpc = ascube(post)
        xc = prior_sample(rng, tpc)
        @test isfinite(logdensityof(tpc, xc))
        yc = Comrade.transform(tpc, xc)
        @test Comrade.inverse(tpc, yc) ≈ xc
    end

    @testset "BrownianMotion posterior" begin
        # random-walk drift from a fixed start with a fitted diffusion coefficient
        @instrument function gmint_bm()
            return @jones begin
                lg ~ ArrayPrior(GaussMarkovStationPrior(ScanSeg(); Ti = BrownianMotion(D = Exponential(0.1)), init = FixedInit(0.0)))
                return Jones(exp(lg))
            end
        end
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_bm())
        s = prior_sample(rng, post)
        chinds = VLBIInstrumentModels.lookup(post.prior.instrument.lg.stationmap)
        @test all(inds -> iszero(parent(s.instrument.lg.params)[first(inds)]), chinds)
        @test all(h -> keys(h) == (:D,), values(s.instrument.lg.hyperparams))
        @test isfinite(logdensityof(post, s))
        tp = asflat(post)
        xf = prior_sample(rng, tp)
        @test isfinite(logdensityof(tp, xf))
        y = Comrade.transform(tp, xf)
        @test Comrade.inverse(tp, y) ≈ xf
        tpc = ascube(post)
        xc = prior_sample(rng, tpc)
        @test isfinite(logdensityof(tpc, xc))

        @testset "Enzyme gradient" begin
            f = let tp = tp
                x -> logdensityof(tp, x)
            end
            gz, = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), Const(f), xf)
            gfd, = grad(central_fdm(5, 1), f, xf)
            @test gz ≈ gfd rtol = 1.0e-5
        end
    end

    @testset "WrappedBrownian phase posterior" begin
        # a circular phase prior; UniformInit covers the per-track offset, so no separate
        # offset term is needed
        @instrument function gmint_wb()
            return @jones begin
                gp ~ ArrayPrior(
                    GaussMarkovStationPrior(ScanSeg(); Ti = WrappedBrownian(τ = InverseGamma(3.0, 6.0)), init = UniformInit());
                    refant = SensitiveReference(0.0)
                )
                return Jones(exp(1im * gp))
            end
        end
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_wb())
        s = prior_sample(rng, post)
        @test all(h -> keys(h) == (:τ,), values(s.instrument.gp.hyperparams))
        @test isfinite(logdensityof(post, s))

        # the model is invariant under a 2π shift of any free phase
        freeinds = setdiff(1:length(post.prior.instrument.gp.dists), post.prior.instrument.gp.dists.fixedinds)
        s2 = deepcopy(s)
        parent(s2.instrument.gp.params)[first(freeinds)] += 2π
        @test logdensityof(post, s2) ≈ logdensityof(post, s)

        tp = asflat(post)
        # one coordinate per free phase, plus one τ per station and the sky's dimensions
        dchain = post.prior.instrument.gp.dists
        nfreew = length(dchain) - length(dchain.fixedinds)
        nstations = length(keys(s.instrument.gp.hyperparams))
        skydim = LogDensityProblems.dimension(asflat(VLBIPosterior(dvis; sky = skym)))
        @test LogDensityProblems.dimension(tp) == skydim + nstations + nfreew

        xf = prior_sample(rng, tp)
        @test isfinite(logdensityof(tp, xf))
        # the coloring wraps, so the roundtrip is checked as `transform ∘ inverse`
        y = Comrade.transform(tp, xf)
        @test all(θ -> -π - 1.0e-12 ≤ θ ≤ π + 1.0e-12, parent(y.instrument.gp.params))
        y2 = Comrade.transform(tp, Comrade.inverse(tp, y))
        @test parent(y2.instrument.gp.params) ≈ parent(y.instrument.gp.params)
        @test all(
            map(
                (a, b) -> a.τ ≈ b.τ,
                values(y2.instrument.gp.hyperparams), values(y.instrument.gp.hyperparams)
            )
        )
        # the flat lift is proper: the sheet weights break the 2π periodicity
        jc = LogDensityProblems.dimension(tp)
        ℓ0 = logdensityof(tp, xf)
        shift(δ) = (z = copy(xf); z[jc] += δ; logdensityof(tp, z))
        # Only the symmetry breaking is checked, not its sign: `jc` need not be a phase
        # coordinate, and a prior sample need not sit on the favored sheet.
        @test !isapprox(shift(2π), ℓ0)
        @test !isapprox(shift(-2π), ℓ0)
        # the density decays far from the origin
        @test shift(40π) < ℓ0 - 10
        @test shift(-40π) < ℓ0 - 10
        # wrapped chains refuse the Std transports
        @test_throws ArgumentError ascube(post)

        @testset "Enzyme gradient" begin
            f = let tp = tp
                x -> logdensityof(tp, x)
            end
            gz, = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), Const(f), xf)
            gfd, = grad(central_fdm(5, 1), f, xf)
            @test gz ≈ gfd rtol = 1.0e-5
        end

        @testset "angle-embedding option" begin
            @instrument function gmint_wbe()
                return @jones begin
                    gp ~ ArrayPrior(
                        GaussMarkovStationPrior(ScanSeg(); Ti = WrappedBrownian(τ = InverseGamma(3.0, 6.0)), init = UniformInit(), param = AngleEmbedded());
                        refant = SensitiveReference(0.0)
                    )
                    return Jones(exp(1im * gp))
                end
            end
            poste = VLBIPosterior(dvis; sky = skym, instrument = gmint_wbe())
            se = prior_sample(rng, poste)
            @test isfinite(logdensityof(poste, se))
            tpe = asflat(poste)
            # each free phase contributes two latent reals, same D and sky dims
            @test LogDensityProblems.dimension(tpe) == skydim + nstations + 2 * nfreew
            xe = prior_sample(rng, tpe)
            @test isfinite(logdensityof(tpe, xe))
            # the embedding drops the latent radius, so the roundtrip is checked on angles
            ye = Comrade.transform(tpe, xe)
            ye2 = Comrade.transform(tpe, Comrade.inverse(tpe, ye))
            @test parent(ye2.instrument.gp.params) ≈ parent(ye.instrument.gp.params)
            # the flat target depends only on the angle, so scaling a latent pair is free
            @test_throws ArgumentError ascube(poste)

            fe = let tpe = tpe
                x -> logdensityof(tpe, x)
            end
            gze, = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), Const(fe), xe)
            gfde, = grad(central_fdm(5, 1), fe, xe)
            @test gze ≈ gfde rtol = 1.0e-5
        end
    end

    @testset "WrappedOrnsteinUhlenbeck ratio-phase posterior" begin
        # a mean-reverting circular process: σ bounds how far the phase strays from its mean
        @instrument function gmint_wou()
            return @jones begin
                gp ~ ArrayPrior(
                    GaussMarkovStationPrior(
                        ScanSeg(); Ti = WrappedOrnsteinUhlenbeck(σ = Exponential(0.3), τ = InverseGamma(3.0, 6.0))
                    );
                    refant = SensitiveReference(0.0)
                )
                return Jones(exp(1im * gp))
            end
        end
        post = VLBIPosterior(dvis; sky = skym, instrument = gmint_wou())
        s = prior_sample(rng, post)
        @test all(h -> keys(h) == (:σ, :τ), values(s.instrument.gp.hyperparams))
        @test isfinite(logdensityof(post, s))
        # prior draws land on the circle
        @test all(θ -> -π - 1.0e-8 ≤ θ ≤ π + 1.0e-8, parent(s.instrument.gp.params))

        # the model is invariant under a 2π shift of any free phase
        freeinds = setdiff(1:length(post.prior.instrument.gp.dists), post.prior.instrument.gp.dists.fixedinds)
        s2 = deepcopy(s)
        parent(s2.instrument.gp.params)[first(freeinds)] += 2π
        @test logdensityof(post, s2) ≈ logdensityof(post, s)

        tp = asflat(post)
        dchain = post.prior.instrument.gp.dists
        nfreew = length(dchain) - length(dchain.fixedinds)
        nstations = length(keys(s.instrument.gp.hyperparams))
        skydim = LogDensityProblems.dimension(asflat(VLBIPosterior(dvis; sky = skym)))
        # one coordinate per free phase, plus σ and τ per station
        @test LogDensityProblems.dimension(tp) == skydim + 2 * nstations + nfreew

        xf = prior_sample(rng, tp)
        @test isfinite(logdensityof(tp, xf))
        y = Comrade.transform(tp, xf)
        @test all(θ -> -π - 1.0e-12 ≤ θ ≤ π + 1.0e-12, parent(y.instrument.gp.params))
        y2 = Comrade.transform(tp, Comrade.inverse(tp, y))
        @test parent(y2.instrument.gp.params) ≈ parent(y.instrument.gp.params)

        # the flat lift is proper: it decays far out instead of repeating forever
        jc = LogDensityProblems.dimension(tp)
        ℓ0 = logdensityof(tp, xf)
        shift(δ) = (z = copy(xf); z[jc] += δ; logdensityof(tp, z))
        @test !isapprox(shift(2π), ℓ0)
        @test shift(40π) < ℓ0 - 10
        @test shift(-40π) < ℓ0 - 10
        @test_throws ArgumentError ascube(post)

        @testset "Enzyme gradient" begin
            f = let tp = tp
                x -> logdensityof(tp, x)
            end
            gz, = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), Const(f), xf)
            gfd, = grad(central_fdm(5, 1), f, xf)
            @test gz ≈ gfd rtol = 1.0e-5
        end

        @testset "angle-embedding option" begin
            @instrument function gmint_woue()
                return @jones begin
                    gp ~ ArrayPrior(
                        GaussMarkovStationPrior(
                            ScanSeg(); Ti = WrappedOrnsteinUhlenbeck(σ = 0.3, τ = InverseGamma(3.0, 6.0)),
                            param = AngleEmbedded()
                        );
                        refant = SensitiveReference(0.0)
                    )
                    return Jones(exp(1im * gp))
                end
            end
            poste = VLBIPosterior(dvis; sky = skym, instrument = gmint_woue())
            se = prior_sample(rng, poste)
            @test isfinite(logdensityof(poste, se))
            tpe = asflat(poste)
            @test LogDensityProblems.dimension(tpe) == skydim + nstations + 2 * nfreew
            xe = prior_sample(rng, tpe)
            @test isfinite(logdensityof(tpe, xe))
            ye = Comrade.transform(tpe, xe)
            ye2 = Comrade.transform(tpe, Comrade.inverse(tpe, ye))
            @test parent(ye2.instrument.gp.params) ≈ parent(ye.instrument.gp.params)

            fe = let tpe = tpe
                x -> logdensityof(tpe, x)
            end
            gze, = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), Const(fe), xe)
            gfde, = grad(central_fdm(5, 1), fe, xe)
            @test gze ≈ gfde rtol = 1.0e-5
        end
    end

end
