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
        @test Comrade._gm_chain_logpdf(p, StationaryInit(), x, inds, ts) ≈ logpdf(dmv, x)
        @inferred Comrade._gm_chain_logpdf(p, StationaryInit(), x, inds, ts)
    end

    @testset "initial priors" begin
        @test_throws ArgumentError GaussianInit(0.0, -1.0)
        @test_throws ArgumentError GaussianInit(0.0, 0.0)
        gi = GaussianInit(0.5, 2.0)
        @test Comrade.initial_moments(gi, p) == (0.5, 4.0)
        @test Comrade.initial_moments(StationaryInit(), p) == Comrade.stationary_moments(p)
        @test FixedInit(0.3).value == 0.3
        # marginal propagation: init moments flow through the transition law
        Φ, Q = Comrade.transition_moments(p, 0.7)
        m, P = Comrade.marginal_moments(gi, p, 0.7)
        @test m ≈ 0.2 + Φ * (0.5 - 0.2)
        @test P ≈ Φ^2 * 4.0 + Q
        @test Comrade.marginal_moments(gi, p, 0.0) == (0.5, 4.0)
        @test Comrade.marginal_moments(StationaryInit(), p, 0.7) == Comrade.stationary_moments(p)
    end

    @testset "extreme correlation time stays finite" begin
        # regression: Δt/τ underflow rounded Φ to exactly 1, so Q = σ²(1-Φ²) was 0 and
        # the chain logpdf was NaN (0/0 + log 0) while the whitening divided by s = 0
        pext = Comrade.OrnsteinUhlenbeck(σ = 0.3, τ = 1.0e20, μ = 0.0)
        Φ, Q = Comrade.transition_moments(pext, 1.0e-3)
        @test Φ < 1
        @test Q > 0
        @test !isnan(Comrade._gm_chain_logpdf(pext, StationaryInit(), x, inds, ts))
        m, s = Comrade._bridge_moments(0.0, Φ * 0.1, Q, Φ, -0.2, Q)
        @test isfinite(m)
        @test s > 0
    end

    @testset "exact conditioning identity" begin
        # scattered, consecutive, endpoints, and fully fixed subsets
        for fp in ([2, 5, 6, 10], [1, 10], [1, 2], [9, 10], collect(1:n))
            frp = setdiff(1:n, fp)
            ℓc = Comrade._gm_chain_logpdf(p, StationaryInit(), x, inds, ts) -
                Comrade._gm_chain_logpdf(p, StationaryInit(), x, inds[fp], ts[fp], ts[fp[1]] - ts[1])
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
        spec = Comrade.MarkovChainSpec(p, StationaryInit(), Val(nothing), inds, ts, fixedpos, false)
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
        spec0 = Comrade.MarkovChainSpec(p, StationaryInit(), Val(nothing), inds, ts, Int[], false)
        d0 = Comrade.GaussMarkovChainDist((AA = spec0,), Int[], Float64[], n)
        z = rand(rng, d0)
        @test logpdf(d0, z) ≈ logpdf(dmv, z)
    end

    @testset "GaussianInit chain vs dense" begin
        # dense reference for a nonstationary start: init N(m0, P0) propagated through
        # the OU transition law, Cov(x_i, x_j) = Φ(|t_j − t_i|) Var(x_min(i,j)) with
        # Var(x_i) = P + Φ(t_i − t_1)² (P0 − P)
        init = GaussianInit(0.5, 2.0)
        P, P0 = 0.3^2, 4.0
        Φt(Δ) = exp(-Δ / 1.7)
        Φ1 = Φt.(ts .- ts[1])
        mg = 0.2 .+ Φ1 .* (0.5 - 0.2)
        Σg = [
            Φt(abs(ts[i] - ts[j])) * (P + Φ1[min(i, j)]^2 * (P0 - P))
                for i in 1:n, j in 1:n
        ]
        dg = MvNormal(mg, Symmetric(Σg))
        @test Comrade._gm_chain_logpdf(p, init, x, inds, ts) ≈ logpdf(dg, x)

        # conditioning identity with a propagated first marginal (first point free)
        fp = [3, 7]
        frp = setdiff(1:n, fp)
        ℓc = Comrade._gm_chain_logpdf(p, init, x, inds, ts) -
            Comrade._gm_chain_logpdf(p, init, x, inds[fp], ts[fp], ts[fp[1]] - ts[1])
        Sfx = Σg[frp, fp]
        Sxx = Σg[fp, fp]
        mc = mg[frp] .+ Sfx * (Sxx \ (x[fp] .- mg[fp]))
        Sc = Symmetric(Σg[frp, frp] .- Sfx * (Sxx \ Sfx'))
        @test ℓc ≈ logpdf(MvNormal(mc, Sc), x[frp])

        # conditional rand exercises the init-marginal prior-from-the-left of the
        # chain-opening free point
        rng2 = Random.Xoshiro(11)
        spec = Comrade.MarkovChainSpec(p, init, Val(nothing), inds, ts, fp, false)
        d = Comrade.GaussMarkovChainDist((AA = spec,), fp, x[fp], n)
        N = 100_000
        draws = Matrix{Float64}(undef, n, N)
        for k in 1:N
            Comrade._rand_chains!(rng2, view(draws, :, k), d, (;))
        end
        @test isapprox(vec(mean(draws[frp, :], dims = 2)), mc; atol = 5.0e-2)
        @test isapprox(cov(draws[frp, :], dims = 2), Matrix(Sc); atol = 1.0e-1)
    end

    @testset "BrownianMotion process" begin
        @test_throws ArgumentError Comrade.BrownianMotion(D = -1.0)
        @test_throws ArgumentError Comrade.BrownianMotion(D = 0.0)
        pb = Comrade.BrownianMotion(D = 0.4)
        @test !Comrade.isstationary(pb)
        Φ, Q = Comrade.transition_moments(pb, 0.5)
        @test Φ == 1
        @test Q ≈ 0.4 * 0.5
        @test_throws ArgumentError Comrade.stationary_moments(pb)
        # a nonstationary process demands an explicit init
        @test_throws ArgumentError GaussMarkovSitePrior(ScanSeg(), pb)
        @test keys(Comrade.hyperprior(Comrade.BrownianMotion(D = Exponential(1.0)))) == (:D,)

        # dense reference: Cov(x_i, x_j) = P0 + D min(t_i − t_1, t_j − t_1)
        init = GaussianInit(0.3, 1.5)
        Σb = [1.5^2 + 0.4 * min(ts[i] - ts[1], ts[j] - ts[1]) for i in 1:n, j in 1:n]
        mb = fill(0.3, n)
        @test Comrade._gm_chain_logpdf(pb, init, x, inds, ts) ≈ logpdf(MvNormal(mb, Symmetric(Σb)), x)

        # conditioning identity at Φ = 1 (Brownian bridge)
        fp = [2, 6]
        frp = setdiff(1:n, fp)
        ℓc = Comrade._gm_chain_logpdf(pb, init, x, inds, ts) -
            Comrade._gm_chain_logpdf(pb, init, x, inds[fp], ts[fp], ts[fp[1]] - ts[1])
        Sfx = Σb[frp, fp]
        Sxx = Σb[fp, fp]
        mc = mb[frp] .+ Sfx * (Sxx \ (x[fp] .- mb[fp]))
        Sc = Symmetric(Σb[frp, frp] .- Sfx * (Sxx \ Sfx'))
        @test ℓc ≈ logpdf(MvNormal(mc, Sc), x[frp])

        # conditional rand through the Φ = 1 bridge
        rng3 = Random.Xoshiro(23)
        spec = Comrade.MarkovChainSpec(pb, init, Val(nothing), inds, ts, fp, false)
        d = Comrade.GaussMarkovChainDist((AA = spec,), fp, x[fp], n)
        N = 100_000
        draws = Matrix{Float64}(undef, n, N)
        for k in 1:N
            Comrade._rand_chains!(rng3, view(draws, :, k), d, (;))
        end
        @test isapprox(vec(mean(draws[frp, :], dims = 2)), mc; atol = 5.0e-2)
        @test isapprox(cov(draws[frp, :], dims = 2), Matrix(Sc); atol = 1.0e-1)
    end

    @testset "WrappedBrownian process" begin
        @test_throws ArgumentError Comrade.WrappedBrownian(D = -1.0)
        pw = Comrade.WrappedBrownian(D = 0.8)
        @test !Comrade.isstationary(pw)
        @test Comrade.is_wrapped(pw)
        @test_throws ArgumentError Comrade.stationary_moments(pw)

        # construction rules
        @test_throws ArgumentError GaussMarkovSitePrior(ScanSeg(), pw)                          # StationaryInit invalid
        @test_throws ArgumentError GaussMarkovSitePrior(ScanSeg(), pw; init = GaussianInit(0.0, 1.0))
        @test_throws ArgumentError GaussMarkovSitePrior(ScanSeg(), p; init = UniformInit())     # circular init on a real-line process
        sp = GaussMarkovSitePrior(ScanSeg(), pw; init = UniformInit())
        @test !sp.centered
        # centered = true opts back into raw-angle coordinates
        spc = GaussMarkovSitePrior(ScanSeg(), pw; init = UniformInit(), centered = true)
        @test spc.centered

        # wrapped-normal density: matches a wide brute-force image sum in both regimes
        for Q in (0.05, 0.5, 2.0, 3.9, 4.1, 8.0, 50.0), Δ in (-3.0, -0.4, 0.0, 1.1, 2.9, 7.0)
            brute = log(sum(k -> exp(-abs2(Δ + 2π * k) / (2Q)), -200:200) / sqrt(2π * Q))
            @test Comrade._wn_logpdf(Δ, Q) ≈ brute rtol = 1.0e-10
        end

        # transitions compose exactly: marginalizing the middle point over the circle
        let Q1 = 0.7, Q2 = 1.3, Δ = 1.9
            g(θm) = exp(Comrade._wn_logpdf(θm, Q1) + Comrade._wn_logpdf(Δ - θm, Q2))
            θs = range(-π, π; length = 20_001)[1:(end - 1)]
            quad = log(sum(g, θs) * step(θs))
            @test quad ≈ Comrade._wn_logpdf(Δ, Q1 + Q2) atol = 1.0e-8
        end

        # exact 2π periodicity of the chain density — the property the wrap buys
        ℓ0 = Comrade._gm_chain_logpdf(pw, UniformInit(), x, inds, ts)
        for i in (1, 4, n)
            xs = copy(x)
            xs[i] += 2π
            @test Comrade._gm_chain_logpdf(pw, UniformInit(), xs, inds, ts) ≈ ℓ0
        end
        @test Comrade._gm_chain_logpdf(pw, UniformInit(), x .+ 2π, inds, ts) ≈ ℓ0

        # unconditional rand: range, uniform start, and the e^{-DΔ/2} coherence decay
        rng4 = Random.Xoshiro(31)
        spec0 = Comrade.MarkovChainSpec(pw, UniformInit(), Val(nothing), inds, ts, Int[], false)
        d0 = Comrade.GaussMarkovChainDist((AA = spec0,), Int[], Float64[], n)
        N = 50_000
        draws = Matrix{Float64}(undef, n, N)
        for k in 1:N
            Comrade._rand_chains!(rng4, view(draws, :, k), d0, (;))
        end
        @test all(abs.(draws) .<= π + 1.0e-8)
        @test abs(mean(cis.(draws[1, :]))) < 0.02
        for (i, j) in ((1, 2), (3, 7), (1, n))
            coh = mean(cis.(draws[j, :] .- draws[i, :]))
            @test abs(coh) ≈ exp(-0.8 * (ts[j] - ts[i]) / 2) atol = 2.0e-2
            @test abs(imag(coh)) < 2.0e-2
        end

        # conditional rand: fixed values hit exactly, in range
        fpw = [2, 6]
        θf = [0.4, -2.9]
        specf = Comrade.MarkovChainSpec(pw, UniformInit(), Val(nothing), inds, ts, fpw, false)
        df = Comrade.GaussMarkovChainDist((AA = specf,), fpw, θf, n)
        drawsf = Matrix{Float64}(undef, n, N)
        for k in 1:N
            Comrade._rand_chains!(rng4, view(drawsf, :, k), df, (;))
        end
        @test all(k -> drawsf[fpw, k] == θf, 1:200)
        @test all(abs.(drawsf[setdiff(1:n, fpw), :]) .<= π + 1.0e-8)

        # single-point bridge matches quadrature, including the winding mixture
        inds3 = [1, 2, 3]
        ts3 = [0.0, 0.4, 1.1]
        θa, θb = 0.3, -2.8
        spec3 = Comrade.MarkovChainSpec(pw, UniformInit(), Val(nothing), inds3, ts3, [1, 3], false)
        d3 = Comrade.GaussMarkovChainDist((AA = spec3,), [1, 3], [θa, θb], 3)
        v = Vector{Float64}(undef, 3)
        dm = map(1:N) do _
            Comrade._rand_chains!(rng4, v, d3, (;))
            v[2]
        end
        Q1, Q2 = 0.8 * 0.4, 0.8 * 0.7
        h(θm) = exp(Comrade._wn_logpdf(θm - θa, Q1) + Comrade._wn_logpdf(θb - θm, Q2))
        θgrid = range(-π, π; length = 4001)[1:(end - 1)]
        cquad = sum(θ -> cis(θ) * h(θ), θgrid) / sum(h, θgrid)
        @test abs(mean(cis.(dm)) - cquad) < 2.0e-2

        # angle-embedding flat transport: two latents per free phase, angles from the
        # atan pairs, and the per-point logjacs are exactly the AngleTransform's
        specw = Comrade.MarkovChainSpec(pw, UniformInit(), Val(nothing), inds, ts, Int[], false)
        dw = Comrade.GaussMarkovChainDist((AA = specw,), Int[], Float64[], n)
        node = Comrade.PT.transport_node(dw, Comrade.PT.TVFlat())
        @test TV.dimension(node) == 2n
        z = randn(rng4, 2n)
        θe, ℓe, _ = TV.transform_with(TV.LogJac(), node, z, 1)
        @test θe ≈ [atan(z[2k - 1], z[2k]) for k in 1:n]
        at = Comrade.PT.angle_transform()
        ℓpair(v) = TV.transform_with(TV.LogJac(), at, v, 1)[2]
        @test ℓe ≈ sum(k -> ℓpair(z[(2k - 1):(2k)]), 1:n)
        # scaling the latents is pure radial freedom: the angles are unchanged
        θs, ℓs, _ = TV.transform_with(TV.LogJac(), node, 2 .* z, 1)
        @test θs ≈ θe
        @test ℓs ≈ sum(k -> ℓpair(2 .* z[(2k - 1):(2k)]), 1:n)
        # inverse embeds at unit radius; transform ∘ inverse is the identity on angles
        zi = zeros(2n)
        TV.inverse_at!(zi, 1, node, θe)
        @test all(hypot(zi[2k - 1], zi[2k]) ≈ 1 for k in 1:n)
        θr, _, _ = TV.transform_with(TV.LogJac(), node, zi, 1)
        @test θr ≈ θe

        # centered = true keeps the raw angles as flat coordinates (identity, zero
        # Jacobian, one coordinate per free phase)
        specc = Comrade.MarkovChainSpec(pw, UniformInit(), Val(nothing), inds, ts, Int[], true)
        dc = Comrade.GaussMarkovChainDist((AA = specc,), Int[], Float64[], n)
        nodec = Comrade.PT.transport_node(dc, Comrade.PT.TVFlat())
        @test TV.dimension(nodec) == n
        θin = 2π .* rand(rng4, n) .- π
        θc, ℓcj, _ = TV.transform_with(TV.LogJac(), nodec, θin, 1)
        @test θc ≈ θin
        @test iszero(ℓcj)
        zc = zeros(n)
        TV.inverse_at!(zc, 1, nodec, θc)
        @test zc ≈ θin
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

    @testset "hierarchical PosteriorSamples postprocessing" begin
        # Exercises the caltable/rmap/siteparams path the StokesI example advertises:
        # samples are `(params, hyperparams)` NamedTuples that must round-trip through
        # PosteriorSamples' StructArray conversion and reduce with `rmap`.
        post = VLBIPosterior(skym, gmint_fitted(), dvis)
        samples = [prior_sample(rng, post) for _ in 1:20]
        chain = Comrade.PosteriorSamples(samples, nothing; metadata = Dict(:sampler => :prior))
        @test length(Comrade.postsamples(chain)) == 20

        mchain = Comrade.rmap(mean, chain)
        schain = Comrade.rmap(std, chain)
        # the hierarchical instrument terms survive reduction as (params, hyperparams) tuples
        @test mchain.instrument.lg.params isa Comrade.SiteArray
        @test keys(mchain.instrument.lg.hyperparams) == Tuple(Comrade.sites(dvis))
        @test Comrade.siteparams(mchain.instrument.lg) === mchain.instrument.lg.params

        # caltable strips the hyperparameters from the reduced chain
        @test caltable(mchain.instrument.lg) isa Comrade.CalTable
        @test caltable(abs.(mchain.instrument.lg.params)) isa Comrade.CalTable

        # instrumentmodel accepts the reduced hierarchical chain (siteparams strips hyperparams)
        @test instrumentmodel(post, mchain) isa Comrade.SiteArray
        # the example's mean/std combine: siteparams unwraps both before the elementwise op
        combined = map(
            (x, y) -> Comrade.siteparams(x) .+ Comrade.siteparams(y),
            mchain.instrument, schain.instrument,
        )
        @test instrumentmodel(post, (; instrument = combined)) isa Comrade.SiteArray
    end

    @testset "empty-group PosteriorSamples (zero-field StructArray guard)" begin
        # regression: an empty NamedTuple leaf (e.g. a sky group with no fitted params)
        # must not error in PosteriorSamples' StructArray conversion (fieldcount guard in
        # `_convert2structarr`).
        raw = [(a = i * 1.0, empty = NamedTuple()) for i in 1:5]
        ps = Comrade.PosteriorSamples(raw, nothing)
        @test length(Comrade.postsamples(ps)) == 5
        @test ps[3].a == 3.0
        @test ps[1].empty == NamedTuple()
    end

    @testset "derived element type" begin
        # eltype is derived from the process/override parameters, not pinned to Float64
        gmp(p) = GaussMarkovSitePrior(ScanSeg(), p)
        @test Comrade._working_type((a = gmp(Comrade.OrnsteinUhlenbeck(σ = 0.1, τ = 2.0)),)) === Float64
        @test Comrade._working_type((a = gmp(Comrade.OrnsteinUhlenbeck(σ = 0.1f0, τ = 2.0f0, μ = 0.0f0)),)) === Float32
        # fitted (Distribution) fields contribute their own eltype
        @test Comrade._working_type((a = gmp(Comrade.OrnsteinUhlenbeck(σ = Exponential(0.1), τ = 2.0)),)) === Float64
        # the init contributes to the working type as well
        @test Comrade._working_type(
            (a = GaussMarkovSitePrior(ScanSeg(), Comrade.OrnsteinUhlenbeck(σ = 0.1f0, τ = 2.0f0, μ = 0.0f0); init = FixedInit(0.0)),)
        ) === Float64
        post = VLBIPosterior(skym, gmint_fitted(), dvis)
        @test eltype(post.prior.instrument.lg) === Float64
        @test eltype(parent(prior_sample(rng, post).instrument.lg.params)) === Float64
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

    @testset "FixedInit chains" begin
        # FixedInit(0.0) conditions every site's chain to zero at its first time, killing
        # the level redundancy with a separate (circular) offset term
        @instrument function gmint_anch()
            return @jones begin
                gp0 ~ ArrayPrior(IIDSitePrior(TrackSeg(), DiagonalVonMises(0.0, inv(π^2))))
                dgp ~ ArrayPrior(
                    GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 1.0, τ = InverseGamma(3.0, 3.0)); init = FixedInit(0.0));
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

        # FixedInit chains are whitened like any other; ascube works (checked without the
        # circular offset term, which correctly refuses Std transport)
        @instrument function gmint_anchonly()
            return @jones begin
                dgp ~ ArrayPrior(
                    GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 1.0, τ = InverseGamma(3.0, 3.0)); init = FixedInit(0.0));
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

        # FixedInit + centered + fitted hyperparameters (the recommended phase setup)
        @instrument function gmint_anchcent()
            return @jones begin
                dgp ~ ArrayPrior(
                    GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 1.0, τ = InverseGamma(3.0, 3.0)); init = FixedInit(0.0), centered = true);
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

    @testset "FixedInit value and refant conflict" begin
        # a nonzero init value propagates exactly through rand and the transform
        @instrument function gmint_fixval()
            return @jones begin
                lg ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 0.5, τ = 2.0); init = FixedInit(0.3)))
                return SingleStokesGain(exp(lg))
            end
        end
        post = VLBIPosterior(skym, gmint_fixval(), dvis)
        s = prior_sample(rng, post)
        chinds = Comrade.lookup(post.prior.instrument.lg.sitemap)
        @test all(inds -> parent(s.instrument.lg)[first(inds)] == 0.3, chinds)
        @test isfinite(logdensityof(post, s))
        tp = asflat(post)
        y = Comrade.transform(tp, prior_sample(rng, tp))
        @test all(inds -> parent(y.instrument.lg)[first(inds)] == 0.3, chinds)

        # a reference scheme that fixes a chain's first point to a different value than
        # FixedInit is a conflict, not a silent precedence
        @instrument function gmint_conflict()
            return @jones begin
                gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 2.0, τ = 1.0); init = FixedInit(0.3)); refant = SEFDReference(0.0))
                return SingleStokesGain(exp(1im * gp))
            end
        end
        @test_throws ArgumentError VLBIPosterior(skym, gmint_conflict(), dvis)

        # agreeing values compose: the reference chain's first point is fixed by both
        @instrument function gmint_agree()
            return @jones begin
                gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 2.0, τ = 1.0); init = FixedInit(0.0)); refant = SEFDReference(0.0))
                return SingleStokesGain(exp(1im * gp))
            end
        end
        post2 = VLBIPosterior(skym, gmint_agree(), dvis)
        @test isfinite(logdensityof(post2, prior_sample(rng, post2)))
    end

    @testset "GaussianInit posterior" begin
        # wide explicit start together with a fitted hyperparameter
        @instrument function gmint_ginit()
            return @jones begin
                lg ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = Exponential(0.1), τ = 2.0); init = GaussianInit(0.0, 1.0)))
                return SingleStokesGain(exp(lg))
            end
        end
        post = VLBIPosterior(skym, gmint_ginit(), dvis)
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
                lg ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), BrownianMotion(D = Exponential(0.1)); init = FixedInit(0.0)))
                return SingleStokesGain(exp(lg))
            end
        end
        post = VLBIPosterior(skym, gmint_bm(), dvis)
        s = prior_sample(rng, post)
        chinds = Comrade.lookup(post.prior.instrument.lg.sitemap)
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
        # a fully circular phase prior: UniformInit absorbs the per-track offset, so no
        # separate circular offset term is needed
        @instrument function gmint_wb()
            return @jones begin
                gp ~ ArrayPrior(
                    GaussMarkovSitePrior(ScanSeg(), WrappedBrownian(D = Exponential(1.0)); init = UniformInit());
                    refant = SEFDReference(0.0)
                )
                return SingleStokesGain(exp(1im * gp))
            end
        end
        post = VLBIPosterior(skym, gmint_wb(), dvis)
        s = prior_sample(rng, post)
        @test all(h -> keys(h) == (:D,), values(s.instrument.gp.hyperparams))
        @test isfinite(logdensityof(post, s))

        # exactly invariant under a 2π shift of any single free phase
        freeinds = setdiff(1:length(post.prior.instrument.gp.dists), post.prior.instrument.gp.dists.fixedinds)
        s2 = deepcopy(s)
        parent(s2.instrument.gp.params)[first(freeinds)] += 2π
        @test logdensityof(post, s2) ≈ logdensityof(post, s)

        tp = asflat(post)
        # each free phase contributes two latent reals (angle embedding), plus one D per
        # site and the sky's dimensions
        dchain = post.prior.instrument.gp.dists
        nfreew = length(dchain) - length(dchain.fixedinds)
        nsites = length(keys(s.instrument.gp.hyperparams))
        skydim = LogDensityProblems.dimension(asflat(VLBIPosterior(skym, dvis)))
        @test LogDensityProblems.dimension(tp) == skydim + nsites + 2 * nfreew

        xf = prior_sample(rng, tp)
        @test isfinite(logdensityof(tp, xf))
        y = Comrade.transform(tp, xf)
        # the angle embedding is not injective (the latent radius is dropped on inverse),
        # so the roundtrip identity holds through the constrained space, not the latent
        y2 = Comrade.transform(tp, Comrade.inverse(tp, y))
        @test parent(y2.instrument.gp.params) ≈ parent(y.instrument.gp.params)
        @test all(
            map(
                (a, b) -> a.D ≈ b.D,
                values(y2.instrument.gp.hyperparams), values(y.instrument.gp.hyperparams)
            )
        )
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

        @testset "centered raw-angle option" begin
            @instrument function gmint_wbc()
                return @jones begin
                    gp ~ ArrayPrior(
                        GaussMarkovSitePrior(ScanSeg(), WrappedBrownian(D = Exponential(1.0)); init = UniformInit(), centered = true);
                        refant = SEFDReference(0.0)
                    )
                    return SingleStokesGain(exp(1im * gp))
                end
            end
            postc = VLBIPosterior(skym, gmint_wbc(), dvis)
            sc = prior_sample(rng, postc)
            @test isfinite(logdensityof(postc, sc))
            tpc = asflat(postc)
            # one coordinate per free phase (no embedding), same D and sky dims
            @test LogDensityProblems.dimension(tpc) == skydim + nsites + nfreew
            xc = prior_sample(rng, tpc)
            @test isfinite(logdensityof(tpc, xc))
            # the centered transport is bijective, so the latent roundtrip is exact
            yc = Comrade.transform(tpc, xc)
            @test Comrade.inverse(tpc, yc) ≈ xc
            # the flat coordinates are the angles themselves: exactly 2π-periodic in
            # every free phase (the identity transport carries the shift through)
            fic = setdiff(1:length(postc.prior.instrument.gp.dists), postc.prior.instrument.gp.dists.fixedinds)
            yc2 = deepcopy(yc)
            parent(yc2.instrument.gp.params)[first(fic)] += 2π
            @test logdensityof(tpc, Comrade.inverse(tpc, yc2)) ≈ logdensityof(tpc, xc)
            @test_throws ArgumentError ascube(postc)

            fc = let tpc = tpc
                x -> logdensityof(tpc, x)
            end
            gzc, = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), Const(fc), xc)
            gfdc, = grad(central_fdm(5, 1), fc, xc)
            @test gzc ≈ gfdc rtol = 1.0e-5
        end
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
