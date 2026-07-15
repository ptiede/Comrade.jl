using Distributions
using VLBIFiles
using LinearAlgebra
using Random
using Statistics
using FiniteDifferences
using Enzyme
using LogDensityProblems
import TransformVariables as TV

function _gm_test_data()
    path = joinpath(@__DIR__, "..", "test_data.uvfits")
    uvd = VLBIFiles.load(VLBIFiles.UVData, path)
    return extract_table(uvd, Visibilities(; time_average = VLBI.GapBasedScans()))
end

function _gm_test_sky()
    tm(θ, meta) = θ.f1 * stretched(Gaussian(), θ.σ1, θ.σ1)
    prior = (f1 = VLBIUniform(0.5, 1.5), σ1 = VLBIUniform(μas2rad(1.0), μas2rad(100.0)))
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 64, 64)
    return SkyModel(tm, prior, g)
end

@testset "GaussMarkov instrument priors" begin
    rng = Random.Xoshiro(42)

    @testset "OrnsteinUhlenbeck process" begin
        p = Comrade.OrnsteinUhlenbeck(σ = 0.3, τ = 1.7, μ = 0.2)
        μ, P = Comrade.stationary_moments(p)
        @test μ == 0.2
        @test P ≈ 0.3^2
        Φ, Q = Comrade.transition_moments(p, 0.5)
        @test Φ ≈ exp(-0.5 / 1.7)
        @test Q ≈ 0.3^2 * (1 - Φ^2)
        @test isempty(Comrade.hyperprior(p))
        @test Comrade.materialize(p, (;)) === p

        ph = Comrade.OrnsteinUhlenbeck(σ = Exponential(0.1), τ = 2.0)
        @test keys(Comrade.hyperprior(ph)) == (:σ,)
        pm = Comrade.materialize(ph, (σ = 0.25,))
        @test pm.σ == 0.25 && pm.τ == 2.0
        @inferred Comrade.materialize(ph, (σ = 0.25,))

        @test_throws ArgumentError Comrade.OrnsteinUhlenbeck(σ = -1.0, τ = 1.0)
        @test_throws ArgumentError Comrade.OrnsteinUhlenbeck(σ = 1.0, τ = 0.0)
        @test_throws ArgumentError Comrade.OrnsteinUhlenbeck(σ = 1.0, τ = 1.0, μ = Normal())
    end

    # dense reference: stationary OU covariance over irregular times
    n = 10
    ts = cumsum(rand(rng, n) .+ 0.05)
    inds = collect(1:n)
    p = Comrade.OrnsteinUhlenbeck(σ = 0.3, τ = 1.7, μ = 0.2)
    Σ = [0.3^2 * exp(-abs(ts[i] - ts[j]) / 1.7) for i in 1:n, j in 1:n]
    dmv = MvNormal(fill(0.2, n), Symmetric(Σ))
    x = randn(rng, n)

    @testset "chain logpdf vs dense" begin
        @test Comrade._gm_chain_logpdf(p, x, inds, ts) ≈ logpdf(dmv, x)
        @inferred Comrade._gm_chain_logpdf(p, x, inds, ts)
    end

    @testset "extreme correlation time stays finite" begin
        # regression: Δt/τ underflow rounded Φ to exactly 1, so Q = σ²(1-Φ²) was 0 and
        # the chain logpdf was NaN (0/0 + log 0) while the whitening divided by s = 0
        pext = Comrade.OrnsteinUhlenbeck(σ = 0.3, τ = 1.0e20, μ = 0.0)
        Φ, Q = Comrade.transition_moments(pext, 1.0e-3)
        @test Φ < 1
        @test Q > 0
        @test !isnan(Comrade._gm_chain_logpdf(pext, x, inds, ts))
        m, s = Comrade._bridge_moments(0.0, 0.3^2, Φ, 0.1, Φ, -0.2)
        @test isfinite(m)
        @test s > 0
    end

    @testset "exact conditioning identity" begin
        # scattered, consecutive, endpoints, and fully fixed subsets
        for fp in ([2, 5, 6, 10], [1, 10], [1, 2], [9, 10], collect(1:n))
            frp = setdiff(1:n, fp)
            ℓc = Comrade._gm_chain_logpdf(p, x, inds, ts) -
                Comrade._gm_chain_logpdf(p, x, inds[fp], ts[fp])
            if isempty(frp)
                @test abs(ℓc) < 1.0e-12
            else
                Sff = Σ[frp, frp]
                Sfx = Σ[frp, fp]
                Sxx = Σ[fp, fp]
                mc = fill(0.2, length(frp)) .+ Sfx * (Sxx \ (x[fp] .- 0.2))
                Sc = Symmetric(Sff .- Sfx * (Sxx \ Sfx'))
                @test ℓc ≈ logpdf(MvNormal(mc, Sc), x[frp])
            end
        end
    end

    @testset "conditional rand (bridge)" begin
        fixedpos = [2, 5, 6, 10]
        freepos = setdiff(1:n, fixedpos)
        spec = Comrade.MarkovChainSpec(p, Val(nothing), inds, ts, fixedpos, false)
        d = Comrade.GaussMarkovChainDist((AA = spec,), fixedpos, x[fixedpos], n)
        N = 100_000
        draws = Matrix{Float64}(undef, n, N)
        for k in 1:N
            Comrade._rand_chains!(rng, view(draws, :, k), d, (;))
        end
        @test all(k -> draws[fixedpos, k] == x[fixedpos], 1:100)
        Σff = Σ[freepos, freepos]
        Σfx = Σ[freepos, fixedpos]
        Σxx = Σ[fixedpos, fixedpos]
        μc = fill(0.2, length(freepos)) .+ Σfx * (Σxx \ (x[fixedpos] .- 0.2))
        Σc = Symmetric(Σff .- Σfx * (Σxx \ Σfx'))
        @test isapprox(vec(mean(draws[freepos, :], dims = 2)), μc; atol = 1.0e-2)
        @test isapprox(cov(draws[freepos, :], dims = 2), Matrix(Σc); atol = 1.0e-2)

        # unconditional rand/logpdf consistency
        spec0 = Comrade.MarkovChainSpec(p, Val(nothing), inds, ts, Int[], false)
        d0 = Comrade.GaussMarkovChainDist((AA = spec0,), Int[], Float64[], n)
        z = rand(rng, d0)
        @test logpdf(d0, z) ≈ logpdf(dmv, z)
    end

    dvis = _gm_test_data()
    skym = _gm_test_sky()

    @instrument function gmint_fitted()
        return @jones begin
            lg ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = Exponential(0.1), τ = Exponential(2.0))))
            gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 2.0, τ = InverseGamma(3.0, 6.0))); refant = SEFDReference(0.0))
            return SingleStokesGain(exp(complex(lg, gp)))
        end
    end

    @testset "posterior with fitted hyperparameters" begin
        post = VLBIPosterior(skym, gmint_fitted(), dvis; admode = set_runtime_activity(Enzyme.Reverse))
        s = prior_sample(rng, post)
        @test s.instrument.lg isa NamedTuple{(:params, :hyperparams)}
        @test s.instrument.lg.params isa Comrade.SiteArray
        # every site gets its own hyperparameters, nested under the site name
        @test keys(s.instrument.lg.hyperparams) == Tuple(Comrade.sites(dvis))
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
        @test ct isa Comrade.CalTable

        # residuals/instrument forward model work on hierarchical samples
        @test Comrade.instrumentmodel(post, s) isa Comrade.SiteArray

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

    @testset "fixed hyperparameters degrade to plain prior" begin
        @instrument function gmint_fixed()
            return @jones begin
                lg ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 0.1, τ = 2.0)))
                gp ~ ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))); refant = SEFDReference(0.0))
                return SingleStokesGain(exp(complex(lg, gp)))
            end
        end
        post = VLBIPosterior(skym, gmint_fixed(), dvis)
        s = prior_sample(rng, post)
        @test s.instrument.lg isa Comrade.SiteArray
        @test isfinite(logdensityof(post, s))
        tp = asflat(post)
        @test isfinite(logdensityof(tp, prior_sample(rng, tp)))
        @test caltable(s.instrument.lg) isa Comrade.CalTable
    end

    @testset "site overrides" begin
        # IID override under a GaussMarkov default
        @instrument function gmint_iidover()
            return @jones begin
                lg ~ ArrayPrior(
                    GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = Exponential(0.1), τ = 2.0));
                    LM = IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 1.0))
                )
                return SingleStokesGain(exp(lg))
            end
        end
        post = VLBIPosterior(skym, gmint_iidover(), dvis)
        s = prior_sample(rng, post)
        # the IID-override site carries no hyperparameters; the GM sites each get their own
        @test :LM ∉ keys(s.instrument.lg.hyperparams)
        @test all(h -> keys(h) == (:σ,), values(s.instrument.lg.hyperparams))
        @test isfinite(logdensityof(post, s))
        tp = asflat(post)
        @test isfinite(logdensityof(tp, prior_sample(rng, tp)))

        # GaussMarkov override with its own fitted hyperparameters, nested under the site
        @instrument function gmint_gmover()
            return @jones begin
                lg ~ ArrayPrior(
                    GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = Exponential(0.1), τ = Exponential(2.0)));
                    LM = GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = Exponential(1.0), τ = 2.0))
                )
                return SingleStokesGain(exp(lg))
            end
        end
        post = VLBIPosterior(skym, gmint_gmover(), dvis)
        s = prior_sample(rng, post)
        @test keys(s.instrument.lg.hyperparams) == Tuple(Comrade.sites(dvis))
        @test keys(s.instrument.lg.hyperparams.AA) == (:σ, :τ)
        # the override site fixes τ, so only σ is fitted there
        @test keys(s.instrument.lg.hyperparams.LM) == (:σ,)
        @test isfinite(logdensityof(post, s))
        tp = asflat(post)
        @test isfinite(logdensityof(tp, prior_sample(rng, tp)))

        # regression: an IID override on the reference site has every point fixed, so its
        # chain term is a sum over an empty index set (previously threw an ArgumentError)
        @instrument function gmint_iidrefover()
            return @jones begin
                gp ~ ArrayPrior(
                    GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 2.0, τ = 1.0));
                    AA = IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 1.0)),
                    refant = SingleReference(:AA, 0.0)
                )
                return SingleStokesGain(exp(1im * gp))
            end
        end
        post = VLBIPosterior(skym, gmint_iidrefover(), dvis)
        s = prior_sample(rng, post)
        @test isfinite(logdensityof(post, s))
        tp = asflat(post)
        @test isfinite(logdensityof(tp, prior_sample(rng, tp)))
    end

    @testset "whitened flat density identity" begin
        # Under exact whitening the flat density of the transported prior restricted to
        # the z block is standard normal: perturbing only the z coordinates changes the
        # flat log-density by exactly −Δ‖z‖²/2. This pins down both the coloring map and
        # its Jacobian, including the hyperparameter coupling and refant conditioning.
        post = VLBIPosterior(skym, gmint_fitted(), dvis)
        nsite = length(Comrade.sites(dvis))
        for (name, nhp) in ((:lg, 2nsite), (:gp, nsite))  # number of hyperparameter coordinates
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
                gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 2.0, τ = 1.0)); refant = SEFDReference(0.0))
                return SingleStokesGain(exp(1im * gp))
            end
        end
        post = VLBIPosterior(skym, gmint_fixref(), dvis)
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
        post = VLBIPosterior(skym, gmint_fitted(), dvis)
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
        # regression: IID overrides must use their own transform (an identity transform
        # would expose latents outside a bounded support)
        @instrument function gmint_boundover()
            return @jones begin
                lg ~ ArrayPrior(
                    GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = Exponential(0.1), τ = 2.0));
                    LM = IIDSitePrior(ScanSeg(), VLBIUniform(-1.0, 1.0))
                )
                return SingleStokesGain(exp(lg))
            end
        end
        post = VLBIPosterior(skym, gmint_boundover(), dvis)
        tp = asflat(post)
        for _ in 1:5
            x = 3 .* randn(rng, dimension(tp))
            @test isfinite(logdensityof(tp, x))
            y = Comrade.transform(tp, x)
            @test all(v -> -1 ≤ v ≤ 1, parent(y.instrument.lg.params)[Comrade.lookup(post.prior.instrument.lg.sitemap).LM])
            @test Comrade.inverse(tp, y) ≈ x
        end
        tpc = ascube(post)
        @test isfinite(logdensityof(tpc, prior_sample(rng, tpc)))
    end

    @testset "centered option" begin
        @instrument function gmint_centered()
            return @jones begin
                lg ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = Exponential(0.1), τ = Exponential(2.0)); centered = true))
                return SingleStokesGain(exp(lg))
            end
        end
        post = VLBIPosterior(skym, gmint_centered(), dvis)
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
        # Simulate gains from a known OU process and check the conditional hyperparameter
        # posterior at the true gains peaks near the truth. (The *joint* MAP is not a
        # valid check here: maximizing a centered hierarchical model over params and
        # hyperparams jointly is a degenerate variance-component estimate. Full NUTS
        # recovery was verified out-of-band; it is too slow for CI.)
        rrng = Random.Xoshiro(7)
        @instrument function gmint_rec()
            return @jones begin
                lg ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = Exponential(0.3), τ = Exponential(2.0))))
                return SingleStokesGain(exp(lg))
            end
        end
        post0 = VLBIPosterior(skym, gmint_rec(), dvis)
        σtrue, τtrue = 0.15, 1.5
        # every site has its own hyperparameters; use the same truth for all of them
        hpof(σ, τ) = site_tuple(dvis, (σ = σ, τ = τ))
        obsprior = post0.prior.instrument.lg
        xg = Vector{Float64}(undef, length(obsprior.dists))
        Comrade._rand_chains!(rrng, xg, obsprior.dists, hpof(σtrue, τtrue))
        gtruth = Comrade.SiteArray(xg, obsprior.sitemap)
        θtrue = (
            sky = (f1 = 1.0, σ1 = μas2rad(40.0)),
            instrument = (lg = (params = gtruth, hyperparams = hpof(σtrue, τtrue)),),
        )
        @test isfinite(logdensityof(post0, θtrue))
        obs = simulate_observation(rrng, post0, θtrue)
        post = VLBIPosterior(skym, gmint_rec(), obs[1])
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

    @testset "anchored chains" begin
        # anchored = true conditions every site's chain to zero at its first time, killing
        # the level redundancy with a separate (circular) offset term
        @instrument function gmint_anch()
            return @jones begin
                gp0 ~ ArrayPrior(IIDSitePrior(TrackSeg(), DiagonalVonMises(0.0, inv(π^2))))
                dgp ~ ArrayPrior(
                    GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 1.0, τ = InverseGamma(3.0, 3.0)); anchored = true);
                    refant = SEFDReference(0.0)
                )
                return SingleStokesGain(exp(1im * (gp0 + dgp)))
            end
        end
        post = VLBIPosterior(skym, gmint_anch(), dvis)
        s = prior_sample(rng, post)
        chinds = Comrade.lookup(post.prior.instrument.dgp.sitemap)
        # every site's first point is exactly zero, in rand and through the transform
        @test all(inds -> iszero(parent(s.instrument.dgp.params)[first(inds)]), chinds)
        @test isfinite(logdensityof(post, s))
        tp = asflat(post)
        xf = prior_sample(rng, tp)
        @test isfinite(logdensityof(tp, xf))
        y = Comrade.transform(tp, xf)
        @test Comrade.inverse(tp, y) ≈ xf
        @test all(inds -> iszero(parent(y.instrument.dgp.params)[first(inds)]), chinds)

        # anchored chains are whitened like any other; ascube works (checked without the
        # circular offset term, which correctly refuses Std transport)
        @instrument function gmint_anchonly()
            return @jones begin
                dgp ~ ArrayPrior(
                    GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 1.0, τ = InverseGamma(3.0, 3.0)); anchored = true);
                    refant = SEFDReference(0.0)
                )
                return SingleStokesGain(exp(1im * dgp))
            end
        end
        postc = VLBIPosterior(skym, gmint_anchonly(), dvis)
        tpc = ascube(postc)
        xc = prior_sample(rng, tpc)
        @test isfinite(logdensityof(tpc, xc))
        yc = Comrade.transform(tpc, xc)
        @test Comrade.inverse(tpc, yc) ≈ xc

        # anchored + centered + fitted hyperparameters (the recommended phase setup)
        @instrument function gmint_anchcent()
            return @jones begin
                dgp ~ ArrayPrior(
                    GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 1.0, τ = InverseGamma(3.0, 3.0)); anchored = true, centered = true);
                    refant = SEFDReference(0.0)
                )
                return SingleStokesGain(exp(1im * dgp))
            end
        end
        postcc = VLBIPosterior(skym, gmint_anchcent(), dvis)
        scc = prior_sample(rng, postcc)
        chindscc = Comrade.lookup(postcc.prior.instrument.dgp.sitemap)
        @test all(inds -> iszero(parent(scc.instrument.dgp.params)[first(inds)]), chindscc)
        tpcc = asflat(postcc)
        xcc = prior_sample(rng, tpcc)
        @test isfinite(logdensityof(tpcc, xcc))
        ycc = Comrade.transform(tpcc, xcc)
        @test Comrade.inverse(tpcc, ycc) ≈ xcc
        @test all(inds -> iszero(parent(ycc.instrument.dgp.params)[first(inds)]), chindscc)
    end

    @testset "phase=true errors" begin
        @instrument function gmint_phase()
            return @jones begin
                gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 2.0, τ = 1.0)); phase = true)
                return SingleStokesGain(exp(1im * gp))
            end
        end
        @test_throws ArgumentError VLBIPosterior(skym, gmint_phase(), dvis)
    end
end
