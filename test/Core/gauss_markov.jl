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
        # non-finite hyperparameters put an exact zero in a transition variance
        @test_throws "must be positive and finite" Comrade.OrnsteinUhlenbeck(σ = Inf, τ = 1.0)
        @test_throws "must be positive and finite" Comrade.OrnsteinUhlenbeck(σ = 1.0, τ = Inf)
        @test_throws "must be positive and finite" Comrade.OrnsteinUhlenbeck(σ = 1.0, τ = NaN)
        @test_throws "must be finite" Comrade.OrnsteinUhlenbeck(σ = 1.0, τ = 1.0, μ = Inf)
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
        # tiny Δt/τ must not round Φ to 1, which would make Q = σ²(1-Φ²) zero
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
        spec = Comrade.MarkovChainSpec(p, StationaryInit(), Val(nothing), inds, ts, fixedpos, NonCentered())
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
        spec0 = Comrade.MarkovChainSpec(p, StationaryInit(), Val(nothing), inds, ts, Int[], NonCentered())
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

        # conditional rand when the chain's first point is free
        rng2 = Random.Xoshiro(11)
        spec = Comrade.MarkovChainSpec(p, init, Val(nothing), inds, ts, fp, NonCentered())
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
        @test_throws "must be positive and finite" Comrade.BrownianMotion(D = Inf)
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
        spec = Comrade.MarkovChainSpec(pb, init, Val(nothing), inds, ts, fp, NonCentered())
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
        @test_throws ArgumentError Comrade.WrappedBrownian(τ = -1.0)
        @test_throws "must be positive and finite" Comrade.WrappedBrownian(τ = Inf)
        τw = 2.5                                   # coherence time; diffusion D = 2/τw = 0.8
        pw = Comrade.WrappedBrownian(τ = τw)
        @test !Comrade.isstationary(pw)
        @test Comrade.is_wrapped(pw)
        @test_throws ArgumentError Comrade.stationary_moments(pw)

        # construction rules
        @test_throws ArgumentError GaussMarkovSitePrior(ScanSeg(), pw)                          # StationaryInit invalid
        @test_throws ArgumentError GaussMarkovSitePrior(ScanSeg(), pw; init = GaussianInit(0.0, 1.0))
        @test_throws ArgumentError GaussMarkovSitePrior(ScanSeg(), p; init = UniformInit())     # circular init on a real-line process
        # WrappedBrownian defaults to non-centered coordinates, as the real-line processes do
        sp = GaussMarkovSitePrior(ScanSeg(), pw; init = UniformInit())
        @test sp.param isa Comrade.NonCentered
        @test GaussMarkovSitePrior(ScanSeg(), p).param isa Comrade.NonCentered
        @test GaussMarkovSitePrior(ScanSeg(), pw; init = UniformInit(), centered = true).param isa Comrade.Centered
        @test GaussMarkovSitePrior(ScanSeg(), pw; init = UniformInit(), centered = false).param isa Comrade.NonCentered
        # the angle embedding applies only to wrapped processes and is never the default
        spe = GaussMarkovSitePrior(ScanSeg(), pw; init = UniformInit(), param = AngleEmbedded())
        @test spe.param isa Comrade.AngleEmbedded
        @test_throws ArgumentError GaussMarkovSitePrior(ScanSeg(), p; param = AngleEmbedded())
        # `param` and `centered` cannot both be given
        @test_throws ArgumentError GaussMarkovSitePrior(
            ScanSeg(), pw; init = UniformInit(), param = Centered(), centered = true
        )

        # wrapped-normal density: matches a wide brute-force image sum in both regimes
        for Q in (0.05, 0.5, 2.0, 3.9, 4.1, 8.0, 50.0), Δ in (-3.0, -0.4, 0.0, 1.1, 2.9, 7.0)
            brute = log(sum(k -> exp(-abs2(Δ + 2π * k) / (2Q)), -200:200) / sqrt(2π * Q))
            @test Comrade._wn_logpdf(Δ, Q) ≈ brute rtol = 1.0e-10
        end

        # a zero variance (a fitted hyperparameter at the edge of its support) gives a
        # finite, very negative log density, not NaN
        @test isfinite(Comrade._wn_logpdf(0.5, 0.0))
        @test Comrade._wn_logpdf(0.5, 0.0) < -1.0e10
        @test isfinite(Comrade._normal_logpdf(0.5, 0.0))
        @test Comrade._normal_logpdf(0.5, 0.0) < -1.0e10
        @test isfinite(Comrade._normal_logpdf(0.0, 0.0))
        @test Comrade._increment_var(pw, 0.0) > 0

        # transitions compose exactly: marginalizing the middle point over the circle
        let Q1 = 0.7, Q2 = 1.3, Δ = 1.9
            g(θm) = exp(Comrade._wn_logpdf(θm, Q1) + Comrade._wn_logpdf(Δ - θm, Q2))
            θs = range(-π, π; length = 20_001)[1:(end - 1)]
            quad = log(sum(g, θs) * step(θs))
            @test quad ≈ Comrade._wn_logpdf(Δ, Q1 + Q2) atol = 1.0e-8
        end

        # the chain density is 2π-periodic in each point
        ℓ0 = Comrade._gm_chain_logpdf(pw, UniformInit(), x, inds, ts)
        for i in (1, 4, n)
            xs = copy(x)
            xs[i] += 2π
            @test Comrade._gm_chain_logpdf(pw, UniformInit(), xs, inds, ts) ≈ ℓ0
        end
        @test Comrade._gm_chain_logpdf(pw, UniformInit(), x .+ 2π, inds, ts) ≈ ℓ0

        # unconditional rand: range, uniform start, and the e^{-Δ/τ} coherence decay
        rng4 = Random.Xoshiro(31)
        spec0 = Comrade.MarkovChainSpec(pw, UniformInit(), Val(nothing), inds, ts, Int[], Centered())
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
            @test abs(coh) ≈ exp(-(ts[j] - ts[i]) / τw) atol = 2.0e-2
            @test abs(imag(coh)) < 2.0e-2
        end

        # conditional rand: fixed values hit exactly, in range
        fpw = [2, 6]
        θf = [0.4, -2.9]
        specf = Comrade.MarkovChainSpec(pw, UniformInit(), Val(nothing), inds, ts, fpw, Centered())
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
        spec3 = Comrade.MarkovChainSpec(pw, UniformInit(), Val(nothing), inds3, ts3, [1, 3], Centered())
        d3 = Comrade.GaussMarkovChainDist((AA = spec3,), [1, 3], [θa, θb], 3)
        v = Vector{Float64}(undef, 3)
        dm = map(1:N) do _
            Comrade._rand_chains!(rng4, v, d3, (;))
            v[2]
        end
        Q1, Q2 = 2 * 0.4 / τw, 2 * 0.7 / τw
        h(θm) = exp(Comrade._wn_logpdf(θm - θa, Q1) + Comrade._wn_logpdf(θb - θm, Q2))
        θgrid = range(-π, π; length = 4001)[1:(end - 1)]
        cquad = sum(θ -> cis(θ) * h(θ), θgrid) / sum(h, θgrid)
        @test abs(mean(cis.(dm)) - cquad) < 2.0e-2

        # angle-embedding flat transport: two latents per free phase, angles from the
        # atan pairs, and the per-point logjacs are exactly the AngleTransform's
        specw = Comrade.MarkovChainSpec(pw, UniformInit(), Val(nothing), inds, ts, Int[], AngleEmbedded())
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

        # the sheet weights of any one coordinate sum to one over its 2π sheets, so they
        # leave the circular model unchanged
        for Q in (0.001, 0.05, 0.5, 2.0, 3.9, 4.1, 10.0, 50.0)
            pq = Comrade.WrappedBrownian(τ = 2 / Q)   # Δt = 1, so the transition variance is Q
            for Δ in (0.0, 0.7, -2.5, 3.0, 6.0)
                tot = sum(m -> exp(Comrade._sheet_logweight(pq, Δ + 2π * m, 0.0, 1.0, 0.0)), -40:40)
                @test tot ≈ 1 rtol = 1.0e-12
            end
        end
        for φ in (0.0, 1.0, -2.0, 3.1, 9.0)
            tot = sum(m -> exp(Comrade._init_sheet_logweight(UniformInit(), pw, φ + 2π * m)), -40:40)
            @test tot ≈ 1 rtol = 1.0e-12
        end
        # a 2π step is suppressed by exp(-(2π)²/2Q)
        let pq = Comrade.WrappedBrownian(τ = 20.0)   # Δt = 1 transition variance 2/20 = 0.1
            w(Δ) = Comrade._sheet_logweight(pq, Δ, 0.0, 1.0, 0.0)
            @test w(2π) - w(0.0) ≈ -(2π)^2 / (2 * 0.1) rtol = 1.0e-6
        end

        # centered = true keeps the raw angles as flat coordinates (one coordinate per
        # free phase) and adds the sheet weight
        specc = Comrade.MarkovChainSpec(pw, UniformInit(), Val(nothing), inds, ts, Int[], Centered())
        dc = Comrade.GaussMarkovChainDist((AA = specc,), Int[], Float64[], n)
        nodec = Comrade.PT.transport_node(dc, Comrade.PT.TVFlat())
        @test TV.dimension(nodec) == n
        θin = 2π .* rand(rng4, n) .- π
        θc, ℓcj, _ = TV.transform_with(TV.LogJac(), nodec, θin, 1)
        @test θc ≈ θin
        @test ℓcj ≈ Comrade._init_sheet_logweight(UniformInit(), pw, θin[1]) +
            sum(k -> Comrade._sheet_logweight(pw, θin[k], θin[k - 1], ts[k] - ts[k - 1], 0.0), 2:n)
        # Coordinates outside (−π, π] wrap to angles, and a 2π shift lowers the density.
        # The baseline must be on the favored sheet: `θin` has steps beyond π, so shifting
        # it can raise the weight, while the smooth `θsm` is its own favored unwrapping.
        θsm = collect(0.25 .* (1:n) .- 1.0)
        _, ℓsm, _ = TV.transform_with(TV.LogJac(), nodec, θsm, 1)
        for j in (3, n)
            θshift = copy(θsm)
            θshift[j] += 2π
            θw, ℓw, _ = TV.transform_with(TV.LogJac(), nodec, θshift, 1)
            @test θw ≈ [Comrade._wrap_angle(v) for v in θshift]
            @test ℓw < ℓsm - 10
            θshift[j] -= 4π
            _, ℓw2, _ = TV.transform_with(TV.LogJac(), nodec, θshift, 1)
            @test ℓw2 < ℓsm - 10
        end
        # the inverse unwraps onto the favored sheet: every step lies in (−π, π]
        zc = zeros(n)
        TV.inverse_at!(zc, 1, nodec, θc)
        @test all(k -> abs(zc[k] - zc[k - 1]) ≤ π + 1.0e-12, 2:n)
        θz, _, _ = TV.transform_with(TV.LogJac(), nodec, zc, 1)
        @test θz ≈ θc
        # with and without reference-fixed points: transform ∘ inverse is the identity on
        # angles, and the pushforward onto the circle is the chain density
        for (fpos, fvals) in ((Int[], Float64[]), ([2], [0.4]), ([1, 4], [0.0, -1.2]))
            sp5 = Comrade.MarkovChainSpec(pw, UniformInit(), Val(nothing), 1:5, ts[1:5], fpos, Centered())
            d5 = Comrade.GaussMarkovChainDist((AA = sp5,), collect(fpos), fvals, 5)
            nd5 = Comrade.PT.transport_node(d5, Comrade.PT.TVFlat())
            nf = TV.dimension(nd5)
            @test nf == 5 - length(fpos)
            z5 = [0.3, -0.2, 0.15, -0.1, 0.25][1:nf]
            y5, _, _ = TV.transform_with(TV.LogJac(), nd5, z5, 1)
            zi5 = zeros(nf)
            TV.inverse_at!(zi5, 1, nd5, y5)
            yr5, _, _ = TV.transform_with(TV.LogJac(), nd5, zi5, 1)
            @test yr5 ≈ y5
            lat = map(Iterators.product(ntuple(_ -> (-4):4, nf)...)) do m
                ys, ℓs, _ = TV.transform_with(TV.LogJac(), nd5, z5 .+ 2π .* collect(m), 1)
                return logpdf(d5, ys) + ℓs
            end
            mx = maximum(lat)
            @test mx + log(sum(l -> exp(l - mx), lat)) ≈ logpdf(d5, y5) atol = 1.0e-9
        end

        # non-centered (the default): each coordinate is a scaled step and the phase is the
        # wrapped running sum; the first free point's coordinate is its angle
        specn = Comrade.MarkovChainSpec(pw, UniformInit(), Val(nothing), inds, ts, Int[], NonCentered())
        dn = Comrade.GaussMarkovChainDist((AA = specn,), Int[], Float64[], n)
        noden = Comrade.PT.transport_node(dn, Comrade.PT.TVFlat())
        @test TV.dimension(noden) == n
        zn = randn(rng4, n)
        θn, ℓn, _ = TV.transform_with(TV.LogJac(), noden, zn, 1)
        φn = cumsum(
            vcat(zn[1], [sqrt(Comrade._increment_var(pw, ts[k] - ts[k - 1])) * zn[k] for k in 2:n])
        )
        @test θn ≈ Comrade._wrap_angle.(φn)
        # the flat target is N(0, I) in the steps
        @test logpdf(dn, θn) + ℓn ≈
            -log(2π) + Comrade._init_sheet_logweight(UniformInit(), pw, zn[1]) +
            sum(k -> logpdf(Normal(), zn[k]), 2:n)
        # `_fused_flat`: the transitions cancel against the sheet weights, so `logpdf`
        # holds only the opening init term
        @test logpdf(dn, θn) ≈ -log(2π)
        # transform ∘ inverse is the identity on wrapped angles
        zi = zeros(n)
        TV.inverse_at!(zi, 1, noden, θn)
        θi, _, _ = TV.transform_with(TV.LogJac(), noden, zi, 1)
        @test θi ≈ θn
        # summed over sheets, the lift gives the circular chain density
        let ts2 = [0.0, 0.8], Q = Comrade._increment_var(pw, 0.8), sq = sqrt(Q)
            sp2 = Comrade.MarkovChainSpec(pw, UniformInit(), Val(nothing), [1, 2], ts2, Int[], NonCentered())
            d2 = Comrade.GaussMarkovChainDist((AA = sp2,), Int[], Float64[], 2)
            nd2 = Comrade.PT.transport_node(d2, Comrade.PT.TVFlat())
            θa2, θb2 = 0.7, -2.1
            tot = sum(Iterators.product((-60):60, (-60):60)) do (m1, m2)
                z2 = [θa2 + 2π * m1, (Comrade._wrap_angle(θb2 - θa2) + 2π * m2) / sq]
                y2, ℓ2, _ = TV.transform_with(TV.LogJac(), nd2, z2, 1)
                # |dθ₂/dz₂| = √Q takes the density in the coordinate to one in the angle
                return exp(logpdf(d2, y2) + ℓ2) / sq
            end
            @test tot ≈ exp(Comrade._wn_logpdf(θb2 - θa2, Q)) / (2π) rtol = 1.0e-9
        end

        # A wrapped chain drops reference-fixed stamps: the transition spans the gap and
        # the fixed value has no prior term.
        let indsg = collect(1:6), tsg = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
            spw = GaussMarkovSitePrior(ScanSeg(), pw; init = UniformInit())
            spr = GaussMarkovSitePrior(ScanSeg(), p)
            cw = Comrade._chainspec(spw, :a, indsg, tsg, [2, 5])
            @test cw.inds == [1, 3, 4, 6]
            @test cw.ts == [0.0, 1.0, 1.5, 2.5]
            @test cw.free.dtl[2] ≈ 1.0          # t₃ − t₁, the composed gap
            @test isempty(cw.fsub.inds)          # nothing left to condition on
            # the real-line process keeps them and conditions on them exactly
            cr = Comrade._chainspec(spr, :a, indsg, tsg, [2, 5])
            @test cr.inds == indsg && cr.fsub.inds == [2, 5]
            # the wrapped chain's density does not depend on the value fixed at stamp 2
            lp(sp, v) = let c = Comrade._chainspec(sp, :a, collect(1:4), tsg[1:4], [2])
                Comrade._chain_logpdf(
                    Comrade.GaussMarkovChainDist(
                        (AA = c,), c.fsub.inds, fill(v, length(c.fsub.inds)), 4
                    ), [0.1, v, -0.4, 1.2], (;)
                )
            end
            @test lp(spw, 0.0) ≈ lp(spw, 2.5)
            @test !(lp(spr, 0.0) ≈ lp(spr, 2.5))
            # a FixedInit opening value is kept and conditioned on
            spf = GaussMarkovSitePrior(ScanSeg(), pw; init = FixedInit(0.0))
            cf = Comrade._chainspec(spf, :a, indsg, tsg, [1, 4])
            @test cf.inds == [1, 2, 3, 5, 6] && cf.fsub.inds == [1]
            # fused remainder of a FixedInit chain is exactly zero for any values
            df = Comrade.GaussMarkovChainDist((AA = cf,), cf.fsub.inds, [0.0], 6)
            @test iszero(logpdf(df, [0.0, 0.4, -0.9, 99.0, 1.2, 0.3]))
            # a site fixed at every stamp has an empty chain; `rand` keeps its fixed values
            ca = Comrade._chainspec(spw, :a, indsg, tsg, collect(1:6))
            @test isempty(ca.inds) && Comrade._nfree(ca) == 0
            da = Comrade.GaussMarkovChainDist((AA = ca,), indsg, fill(0.3, 6), 6)
            @test Comrade._flat_dim(da) == 0
            @test iszero(Comrade._chain_logpdf(da, fill(0.3, 6), (;)))
            va = Vector{Float64}(undef, 6)
            Comrade._rand_chains!(Random.Xoshiro(3), va, da, (;))
            @test va == fill(0.3, 6)
        end
    end

    @testset "WrappedOrnsteinUhlenbeck process" begin
        @test_throws ArgumentError Comrade.WrappedOrnsteinUhlenbeck(σ = -1.0, τ = 1.0)
        @test_throws ArgumentError Comrade.WrappedOrnsteinUhlenbeck(σ = 1.0, τ = -1.0)
        @test_throws ArgumentError Comrade.WrappedOrnsteinUhlenbeck(σ = 1.0, τ = 1.0, μ = Normal())
        @test_throws "must be positive and finite" Comrade.WrappedOrnsteinUhlenbeck(σ = 1.0, τ = Inf)

        σw, τw2, μw = 0.3, 2.0, 0.4
        pr = Comrade.WrappedOrnsteinUhlenbeck(σ = σw, τ = τw2, μ = μw)
        @test Comrade.is_wrapped(pr)
        @test Comrade.isstationary(pr)
        @test Comrade._periodic_mean(pr)
        @test Comrade.stationary_moments(pr) == (μw, σw^2)
        @test Comrade._process_mean(pr) == μw
        Φ, Q = Comrade.transition_moments(pr, 0.5)
        @test Φ ≈ exp(-0.5 / τw2)
        @test Q ≈ σw^2 * (1 - Φ^2)
        @test keys(Comrade.hyperprior(Comrade.WrappedOrnsteinUhlenbeck(σ = Exponential(0.2), τ = 2.0))) == (:σ,)

        # with the shortest-arc drift the transition law is a normalized density on the
        # circle, 2π-periodic in both arguments
        θgrid = range(-π, π; length = 20_001)[1:(end - 1)]
        ker(p, xi, xp, Δt) = exp(Comrade._transition_logpdf(p, xi, xp, Δt, Comrade._process_mean(p)))
        for xp in (-3.0, -0.5, 0.9, 3.1), Δt in (0.1, 1.0, 20.0)
            @test sum(θ -> ker(pr, θ, xp, Δt), θgrid) * step(θgrid) ≈ 1 atol = 1.0e-9
            for xi in (-1.0, 0.7, 2.5)
                ℓ = Comrade._transition_logpdf(pr, xi, xp, Δt, μw)
                @test ℓ ≈ Comrade._transition_logpdf(pr, xi + 2π, xp, Δt, μw)
                @test ℓ ≈ Comrade._transition_logpdf(pr, xi, xp - 2π, Δt, μw)
            end
        end

        # the WrappedBrownian limit: σ, τ → ∞ at fixed 2σ²/τ = 2/τ_coh
        let τc = 1.7
            errs = map((30.0, 300.0)) do sc
                q = Comrade.WrappedOrnsteinUhlenbeck(σ = sc, τ = sc^2 * τc)
                return abs(
                    Comrade._transition_logpdf(q, 0.4, 0.1, 0.5, 0.0) -
                        Comrade._transition_logpdf(Comrade.WrappedBrownian(τ = τc), 0.4, 0.1, 0.5, 0.0)
                )
            end
            @test errs[1] < 1.0e-3
            @test errs[2] < 100 * eps() + errs[1] / 50   # falls off as 1/σ²
        end

        # WN(μ, σ²) is the stationary marginal and the kernel composes over a skipped
        # stamp. Both hold up to the WN mass beyond μ ± π, so the error grows with σ;
        # the σ = 2 bound is the one in the WrappedOrnsteinUhlenbeck docstring.
        wn(p, δ) = exp(Comrade._wn_logpdf(δ, p.σ^2))
        push0 = map(θgrid) do xi
            sum(θ -> ker(pr, xi, θ, 0.9) * wn(pr, θ - μw), θgrid) * step(θgrid)
        end
        @test maximum(abs.(push0 .- wn.(Ref(pr), θgrid .- μw))) < 1.0e-8
        function comperr(p, d1, d2, x0)
            μ0 = Comrade._process_mean(p)
            return maximum(θgrid) do xi
                twice = sum(θ -> ker(p, xi, θ, d2) * ker(p, θ, x0, d1), θgrid) * step(θgrid)
                return abs(twice - ker(p, xi, x0, d1 + d2))
            end
        end
        @test comperr(pr, 0.7, 1.3, 0.4) < 1.0e-8
        @test comperr(Comrade.WrappedOrnsteinUhlenbeck(σ = 1.0, τ = τw2), 0.7, 1.3, 0.4) < 1.0e-4
        @test comperr(Comrade.WrappedOrnsteinUhlenbeck(σ = 2.0, τ = τw2), 0.7, 1.3, 0.4) < 5.0e-3

        # sheet weights sum to one around the mean-reverting center, and the first
        # point's weights sum to one over the WN(μ, σ²) window
        for Δt in (0.05, 1.0, 8.0), xp in (-2.0, 0.4, 2.7)
            tot = sum(m -> exp(Comrade._sheet_logweight(pr, 0.37 + 2π * m, xp, Δt, μw)), -40:40)
            @test tot ≈ 1 rtol = 1.0e-12
        end
        for φ in (0.0, 1.0, -2.0, 3.1)
            tot = sum(m -> exp(Comrade._init_sheet_logweight(StationaryInit(), pr, φ + 2π * m)), -40:40)
            @test tot ≈ 1 rtol = 1.0e-12
        end
        # the unused window of a FixedInit opening point is finite, not 0 * NaN
        @test isfinite(Comrade._init_sheet_logweight(FixedInit(0.0), pr, 0.3))

        # construction rules: uniform is not the stationary law here, so it is rejected;
        # StationaryInit is the default
        @test_throws ArgumentError GaussMarkovSitePrior(ScanSeg(), pr; init = UniformInit())
        @test_throws ArgumentError GaussMarkovSitePrior(ScanSeg(), pr; init = GaussianInit(0.0, 1.0))
        @test GaussMarkovSitePrior(ScanSeg(), pr).param isa Comrade.Centered
        @test GaussMarkovSitePrior(ScanSeg(), pr).init isa StationaryInit
        @test GaussMarkovSitePrior(ScanSeg(), pr; init = FixedInit(0.0)) isa GaussMarkovSitePrior

        # the chain density is exactly 2π-periodic in every point, like WrappedBrownian
        nr = 8
        tsr = cumsum([0.0; 0.2 .* ones(nr - 1)])
        indsr = collect(1:nr)
        xr = 0.5 .* randn(Random.Xoshiro(7), nr)
        ℓr = Comrade._gm_chain_logpdf(pr, StationaryInit(), xr, indsr, tsr)
        for i in (1, 4, nr)
            xs = copy(xr)
            xs[i] += 2π
            @test Comrade._gm_chain_logpdf(pr, StationaryInit(), xs, indsr, tsr) ≈ ℓr
        end

        # rand: values on the circle, the WN(μ, σ²) marginal, and the saturating
        # coherence exp(-σ²(1 - e^{-Δ/τ}))
        rngr = Random.Xoshiro(19)
        specr = Comrade.MarkovChainSpec(pr, StationaryInit(), Val(nothing), indsr, tsr, Int[], Centered())
        dr = Comrade.GaussMarkovChainDist((AA = specr,), Int[], Float64[], nr)
        N = 50_000
        dw = Matrix{Float64}(undef, nr, N)
        for k in 1:N
            Comrade._rand_chains!(rngr, view(dw, :, k), dr, (;))
        end
        @test all(abs.(dw) .<= π + 1.0e-8)
        @test abs(mean(cis.(dw[1, :])) - cis(μw) * exp(-σw^2 / 2)) < 2.0e-2
        for (i, j) in ((1, 2), (3, 7), (1, nr))
            coh = mean(cis.(dw[j, :] .- dw[i, :]))
            @test abs(coh) ≈ exp(-σw^2 * (1 - exp(-(tsr[j] - tsr[i]) / τw2))) atol = 2.0e-2
        end
        # conditional rand hits the reference-fixed values exactly
        fpr = [2, 6]
        θfr = [0.4, -0.3]
        specfr = Comrade.MarkovChainSpec(pr, StationaryInit(), Val(nothing), indsr, tsr, fpr, Centered())
        dfr = Comrade.GaussMarkovChainDist((AA = specfr,), fpr, θfr, nr)
        vr = Vector{Float64}(undef, nr)
        for _ in 1:200
            Comrade._rand_chains!(rngr, vr, dfr, (;))
            @test vr[fpr] == θfr
            @test all(abs.(vr) .<= π + 1.0e-8)
        end

        # centered flat lift: sheet-weighted raw angles; the inverse unwraps each point
        # independently, since the transition mean is periodic in the previous angle
        nodew = Comrade.PT.transport_node(dr, Comrade.PT.TVFlat())
        @test TV.dimension(nodew) == nr
        θin = 2π .* rand(rngr, nr) .- π
        θw, ℓw, _ = TV.transform_with(TV.LogJac(), nodew, θin, 1)
        @test θw ≈ θin
        @test ℓw ≈ Comrade._init_sheet_logweight(StationaryInit(), pr, θin[1]) +
            sum(k -> Comrade._sheet_logweight(pr, θin[k], θin[k - 1], tsr[k] - tsr[k - 1], μw), 2:nr)
        zw = zeros(nr)
        TV.inverse_at!(zw, 1, nodew, θw)
        θz, _, _ = TV.transform_with(TV.LogJac(), nodew, zw, 1)
        @test θz ≈ θw
        # every coordinate lands within π of its own transition mean
        @test abs(zw[1] - μw) ≤ π + 1.0e-12
        for k in 2:nr
            Φk, _ = Comrade.transition_moments(pr, tsr[k] - tsr[k - 1])
            @test abs(zw[k] - Comrade._cond_mean(pr, θw[k - 1], Φk, μw)) ≤ π + 1.0e-12
        end
        # the lift is proper and its pushforward onto the circle is the chain density,
        # with and without reference-fixed points
        for (fpos, fvals) in ((Int[], Float64[]), ([2], [0.4]), ([1, 4], [0.0, -1.2]))
            sp5 = Comrade.MarkovChainSpec(pr, StationaryInit(), Val(nothing), 1:5, tsr[1:5], fpos, Centered())
            d5 = Comrade.GaussMarkovChainDist((AA = sp5,), collect(fpos), fvals, 5)
            nd5 = Comrade.PT.transport_node(d5, Comrade.PT.TVFlat())
            nf = TV.dimension(nd5)
            @test nf == 5 - length(fpos)
            z5 = [0.3, -0.2, 0.15, -0.1, 0.25][1:nf]
            y5, _, _ = TV.transform_with(TV.LogJac(), nd5, z5, 1)
            zi5 = zeros(nf)
            TV.inverse_at!(zi5, 1, nd5, y5)
            yr5, _, _ = TV.transform_with(TV.LogJac(), nd5, zi5, 1)
            @test yr5 ≈ y5
            lat = map(Iterators.product(ntuple(_ -> (-4):4, nf)...)) do m
                ys, ℓs, _ = TV.transform_with(TV.LogJac(), nd5, z5 .+ 2π .* collect(m), 1)
                return logpdf(d5, ys) + ℓs
            end
            mx = maximum(lat)
            @test mx + log(sum(l -> exp(l - mx), lat)) ≈ logpdf(d5, y5) atol = 1.0e-9
            # proper: shifting one coordinate far off is strongly suppressed
            _, ℓ5, _ = TV.transform_with(TV.LogJac(), nd5, z5, 1)
            zsh = copy(z5)
            zsh[1] += 2π
            _, ℓsh, _ = TV.transform_with(TV.LogJac(), nd5, zsh, 1)
            @test ℓsh < ℓ5 - 10
        end
    end

    @testset "VonMisesProcess" begin
        @test_throws ArgumentError Comrade.VonMisesProcess(σ = -1.0, τ = 1.0)
        @test_throws ArgumentError Comrade.VonMisesProcess(σ = 1.0, τ = -1.0)
        @test_throws ArgumentError Comrade.VonMisesProcess(σ = 1.0, τ = 1.0, μ = Normal())
        @test_throws "must be positive and finite" Comrade.VonMisesProcess(σ = 1.0, τ = Inf)

        σv, τv, μv = 0.3, 2.0, 0.4
        pv = Comrade.VonMisesProcess(σ = σv, τ = τv, μ = μv)
        @test Comrade.is_wrapped(pv)
        @test Comrade.isstationary(pv)
        @test !Comrade._periodic_mean(pv)
        @test !Comrade.has_affine_lift(pv)
        @test Comrade.has_noncentered_lift(pv)
        @test Comrade.stationary_moments(pv) == (μv, σv^2)
        @test Comrade.transition_moments(pv, 0.5) ==
            Comrade.transition_moments(Comrade.WrappedOrnsteinUhlenbeck(σ = σv, τ = τv), 0.5)
        @test keys(Comrade.hyperprior(Comrade.VonMisesProcess(σ = Exponential(0.2), τ = 2.0))) == (:σ,)

        # the sine drift is 2π-equivariant, monotone, and within (1 − Φ)|ε|³/6 of the
        # shortest-arc drift for small ε
        pw = Comrade.WrappedOrnsteinUhlenbeck(σ = σv, τ = τv, μ = μv)
        for Δt in (0.1, 1.0, 20.0), xp in (-3.0, -0.5, 0.9, 3.1)
            Φ, _ = Comrade.transition_moments(pv, Δt)
            m = Comrade._cond_mean(pv, xp, Φ, μv)
            @test Comrade._cond_mean(pv, xp + 2π, Φ, μv) ≈ m + 2π
            @test Comrade._cond_mean(pv, xp + 1.0e-4, Φ, μv) > m
            ε = xp - μv
            abs(Comrade._wrap_angle(ε)) < 1.5 && @test abs(m - Comrade._cond_mean(pw, xp, Φ, μv)) ≤
                (1 - Φ) * abs(Comrade._wrap_angle(ε))^3 / 6 + 1.0e-12
        end

        # the transition law is a normalized circular kernel, 2π-periodic in both
        # arguments, as for the shortest-arc process
        θgrid = range(-π, π; length = 20_001)[1:(end - 1)]
        ker(p, xi, xp, Δt) = exp(Comrade._transition_logpdf(p, xi, xp, Δt, Comrade._process_mean(p)))
        for xp in (-3.0, 0.9), Δt in (0.1, 1.0)
            @test sum(θ -> ker(pv, θ, xp, Δt), θgrid) * step(θgrid) ≈ 1 atol = 1.0e-9
            ℓ = Comrade._transition_logpdf(pv, 0.7, xp, Δt, μv)
            @test ℓ ≈ Comrade._transition_logpdf(pv, 0.7 + 2π, xp, Δt, μv)
            @test ℓ ≈ Comrade._transition_logpdf(pv, 0.7, xp - 2π, Δt, μv)
        end

        # construction rules: NonCentered is the default; uniform is rejected; the
        # wrapped-only forms behave as for WOU
        @test GaussMarkovSitePrior(ScanSeg(), pv).param isa Comrade.NonCentered
        @test GaussMarkovSitePrior(ScanSeg(), pv).init isa StationaryInit
        @test GaussMarkovSitePrior(ScanSeg(), pv; centered = true).param isa Comrade.Centered
        @test GaussMarkovSitePrior(ScanSeg(), pv; param = AngleEmbedded()).param isa Comrade.AngleEmbedded
        @test_throws ArgumentError GaussMarkovSitePrior(ScanSeg(), pv; init = UniformInit())
        @test_throws ArgumentError GaussMarkovSitePrior(ScanSeg(), pv; init = GaussianInit(0.0, 1.0))
        # the shortest-arc process rejects NonCentered and names VonMisesProcess instead
        @test_throws "VonMisesProcess" GaussMarkovSitePrior(ScanSeg(), pw; centered = false)

        nv = 8
        tsv = cumsum([0.0; 0.2 .* ones(nv - 1)]) .+ 0.1
        indsv = collect(1:nv)
        rngv = Random.Xoshiro(23)

        # the chain density is 2π-periodic in every point and close to the shortest-arc
        # chain density on a prior draw
        specn = Comrade.MarkovChainSpec(pv, StationaryInit(), Val(nothing), indsv, tsv, Int[], NonCentered())
        dn = Comrade.GaussMarkovChainDist((AA = specn,), Int[], Float64[], nv)
        xv = Vector{Float64}(undef, nv)
        Comrade._rand_chains!(rngv, xv, dn, (;))
        ℓv = Comrade._gm_chain_logpdf(pv, StationaryInit(), xv, indsv, tsv)
        for i in (1, 4, nv)
            xs = copy(xv)
            xs[i] += 2π
            @test Comrade._gm_chain_logpdf(pv, StationaryInit(), xs, indsv, tsv) ≈ ℓv
        end
        @test abs(ℓv - Comrade._gm_chain_logpdf(pw, StationaryInit(), xv, indsv, tsv)) < 0.05

        # non-centered flat lift: θₖ = wrap(m(θₖ₋₁) + √Q zₖ), computed sequentially, so the
        # chain takes no affine-scan slots
        @test isempty(Comrade._scan_tgt(specn))
        noden = Comrade.PT.transport_node(dn, Comrade.PT.TVFlat())
        @test TV.dimension(noden) == nv
        zn = randn(rngv, nv)
        θn, ℓn, _ = TV.transform_with(TV.LogJac(), noden, zn, 1)
        φ = zn[1]
        θman = [Comrade._wrap_angle(φ)]
        for k in 2:nv
            Φk, Qk = Comrade.transition_moments(pv, tsv[k] - tsv[k - 1])
            φ = Comrade._wrap_angle(Comrade._cond_mean(pv, θman[k - 1], Φk, μv) + sqrt(Qk) * zn[k])
            push!(θman, φ)
        end
        @test θn ≈ θman
        # flat target: N(0, 1) per step plus the opening point's WN marginal (in `logpdf`)
        # and sheet weight (in the transform)
        @test logpdf(dn, θn) + ℓn ≈
            Comrade._wn_logpdf(zn[1] - μv, σv^2) +
            Comrade._init_sheet_logweight(StationaryInit(), pv, zn[1]) +
            sum(k -> logpdf(Normal(), zn[k]), 2:nv)
        # transform ∘ inverse is the identity on wrapped angles
        zi = zeros(nv)
        TV.inverse_at!(zi, 1, noden, θn)
        θr, _, _ = TV.transform_with(TV.LogJac(), noden, zi, 1)
        @test θr ≈ θn
        # summed over sheets, the lift gives the circular chain density
        let ts2 = [0.0, 0.8], (Φ2, Q2) = Comrade.transition_moments(pv, 0.8), sq = sqrt(Q2)
            sp2 = Comrade.MarkovChainSpec(pv, StationaryInit(), Val(nothing), [1, 2], ts2, Int[], NonCentered())
            d2 = Comrade.GaussMarkovChainDist((AA = sp2,), Int[], Float64[], 2)
            nd2 = Comrade.PT.transport_node(d2, Comrade.PT.TVFlat())
            θa, θb = 0.7, -0.9
            m2 = Comrade._cond_mean(pv, θa, Φ2, μv)
            tot = sum(Iterators.product((-40):40, (-40):40)) do (m1, mm)
                z2 = [θa + 2π * m1, (Comrade._wrap_angle(θb - m2) + 2π * mm) / sq]
                y2, ℓ2, _ = TV.transform_with(TV.LogJac(), nd2, z2, 1)
                return exp(logpdf(d2, y2) + ℓ2) / sq
            end
            @test tot ≈ exp(Comrade._wn_logpdf(θa - μv, σv^2) + Comrade._wn_logpdf(θb - m2, Q2)) rtol = 1.0e-9
        end

        # Enzyme gradient of the flat density matches finite differences
        let node = noden, d = dn
            f = z -> begin
                y, ℓ, _ = TV.transform_with(TV.LogJac(), node, z, 1)
                return logpdf(d, y) + ℓ
            end
            gfd = grad(central_fdm(5, 1), f, zn)[1]
            gz = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), Const(f), zn)[1]
            @test gz ≈ gfd rtol = 1.0e-6
        end

        # centered flat lift: sheet-weighted raw angles; the inverse puts every step
        # within π of its own mean
        specc = Comrade.MarkovChainSpec(pv, StationaryInit(), Val(nothing), indsv, tsv, Int[], Centered())
        dc = Comrade.GaussMarkovChainDist((AA = specc,), Int[], Float64[], nv)
        nodec = Comrade.PT.transport_node(dc, Comrade.PT.TVFlat())
        θin = 2π .* rand(rngv, nv) .- π
        θc, ℓcw, _ = TV.transform_with(TV.LogJac(), nodec, θin, 1)
        @test θc ≈ θin
        @test ℓcw ≈ Comrade._init_sheet_logweight(StationaryInit(), pv, θin[1]) +
            sum(k -> Comrade._sheet_logweight(pv, θin[k], θin[k - 1], tsv[k] - tsv[k - 1], μv), 2:nv)
        zc = zeros(nv)
        TV.inverse_at!(zc, 1, nodec, θc)
        θz, _, _ = TV.transform_with(TV.LogJac(), nodec, zc, 1)
        @test θz ≈ θc
        @test abs(zc[1] - μv) ≤ π + 1.0e-12
        for k in 2:nv
            Φk, _ = Comrade.transition_moments(pv, tsv[k] - tsv[k - 1])
            @test abs(zc[k] - Comrade._cond_mean(pv, zc[k - 1], Φk, μv)) ≤ π + 1.0e-12
        end

        # sheet weights sum to one around the sine-drift mean, so the lift is exact
        for Δt in (0.05, 1.0), xp in (-2.0, 0.4)
            tot = sum(m -> exp(Comrade._sheet_logweight(pv, 0.37 + 2π * m, xp, Δt, μv)), -40:40)
            @test tot ≈ 1 rtol = 1.0e-12
        end

        # a FixedInit chain starts from its fixed opening value; the fused remainder is zero
        let spf = GaussMarkovSitePrior(ScanSeg(), pv; init = FixedInit(0.0))
            cf = Comrade._chainspec(spf, :a, indsv, tsv, [1])
            @test cf.inds == indsv && cf.fsub.inds == [1]
            df = Comrade.GaussMarkovChainDist((AA = cf,), [1], [0.0], nv)
            @test iszero(logpdf(df, xv))
            ndf = Comrade.PT.transport_node(df, Comrade.PT.TVFlat())
            @test TV.dimension(ndf) == nv - 1
            zf = 0.3 .* randn(rngv, nv - 1)
            yf, ℓf, _ = TV.transform_with(TV.LogJac(), ndf, zf, 1)
            @test yf[1] == 0.0
            Φ1, Q1 = Comrade.transition_moments(pv, tsv[2] - tsv[1])
            @test yf[2] ≈ Comrade._wrap_angle(Comrade._cond_mean(pv, 0.0, Φ1, μv) + sqrt(Q1) * zf[1])
            @test ℓf ≈ sum(k -> logpdf(Normal(), zf[k]), 1:(nv - 1))
            zfi = zeros(nv - 1)
            TV.inverse_at!(zfi, 1, ndf, yf)
            yfr, _, _ = TV.transform_with(TV.LogJac(), ndf, zfi, 1)
            @test yfr ≈ yf
        end

        # interior reference-fixed stamps are dropped and the transition spans the gap
        let spv = GaussMarkovSitePrior(ScanSeg(), pv)
            cg = Comrade._chainspec(spv, :a, indsv, tsv, [2, 5])
            @test cg.inds == [1, 3, 4, 6, 7, 8]
            @test isempty(cg.fsub.inds)
            @test cg.free.dtl[2] ≈ tsv[3] - tsv[1]
        end

        # batched group walkers (the path Reactant traces) agree with the per-chain ones,
        # including per-site numeric overrides
        let inds2 = collect((nv + 1):(2nv)), pv2 = Comrade.VonMisesProcess(σ = 0.5, τ = 3.0, μ = μv)
            spa = Comrade.MarkovChainSpec(pv, StationaryInit(), Val(nothing), indsv, tsv, Int[], NonCentered())
            spb = Comrade.MarkovChainSpec(pv2, StationaryInit(), Val(nothing), inds2, tsv, Int[], NonCentered())
            dg = Comrade.GaussMarkovChainDist((AA = spa, AB = spb), Int[], Float64[], 2nv)
            gs = Comrade.chaingroups(dg)
            @test length(gs) == 1 && first(gs) isa Comrade.MarkovGroup
            zg = randn(rngv, 2nv)
            # per-chain flat coloring (the host path)
            yh = zeros(2nv)
            avh, cvh = Comrade._scan_buffers(yh, dg)
            ℓh, _, _ = Comrade._color_specs_flat!(TV.LogJac(), yh, zg, 1, 1, avh, cvh, values(Comrade.chainspecs(dg)), (;))
            Comrade._resolve_scan!(yh, dg, avh, cvh)
            # batched group coloring (the traced path)
            yg = zeros(2nv)
            avg, cvg = Comrade._scan_buffers(yg, dg)
            ℓg, _, _ = Comrade._color_specs_flat!(TV.LogJac(), yg, zg, 1, 1, avg, cvg, gs, (;))
            Comrade._resolve_scan!(yg, dg, avg, cvg)
            @test yg ≈ yh
            @test ℓg ≈ ℓh
            @test Comrade.chain_term(first(gs), yh, (;)) ≈
                Comrade.chain_term(spa, yh, (;)) + Comrade.chain_term(spb, yh, (;))
        end

        # fitted hyperparameters: the whitened coordinates do not depend on σ, so only the
        # opening term changes
        let ph = Comrade.VonMisesProcess(σ = Exponential(0.3), τ = τv, μ = μv)
            sph = Comrade.MarkovChainSpec(ph, StationaryInit(), Val(:AA), indsv, tsv, Int[], NonCentered())
            dh = Comrade.GaussMarkovChainDist((AA = sph,), Int[], Float64[], nv)
            hp = (AA = (σ = σv,),)
            yh = zeros(nv)
            avh, cvh = Comrade._scan_buffers(yh, dh)
            ℓh, _, _ = Comrade._color_specs_flat!(TV.LogJac(), yh, zn, 1, 1, avh, cvh, values(Comrade.chainspecs(dh)), hp)
            Comrade._resolve_scan!(yh, dh, avh, cvh)
            @test yh ≈ θn
            @test ℓh ≈ ℓn
            @test Comrade._chain_logpdf(dh, yh, hp) ≈ logpdf(dn, θn)
        end

        # rand: on the circle, the WN(μ, σ²) marginal, and reference-fixed values hit
        N = 50_000
        dwv = Matrix{Float64}(undef, nv, N)
        for k in 1:N
            Comrade._rand_chains!(rngv, view(dwv, :, k), dn, (;))
        end
        @test all(abs.(dwv) .<= π + 1.0e-8)
        @test abs(mean(cis.(dwv[1, :])) - cis(μv) * exp(-σv^2 / 2)) < 2.0e-2
        for (i, j) in ((1, 2), (3, 7))
            coh = mean(cis.(dwv[j, :] .- dwv[i, :]))
            @test abs(coh) ≈ exp(-σv^2 * (1 - exp(-(tsv[j] - tsv[i]) / τv))) atol = 2.0e-2
        end
        fpv = [2, 6]
        θfv = [0.4, -0.3]
        specfv = Comrade.MarkovChainSpec(pv, StationaryInit(), Val(nothing), indsv, tsv, fpv, NonCentered())
        dfv = Comrade.GaussMarkovChainDist((AA = specfv,), fpv, θfv, nv)
        vv = Vector{Float64}(undef, nv)
        for _ in 1:200
            Comrade._rand_chains!(rngv, vv, dfv, (;))
            @test vv[fpv] == θfv
            @test all(abs.(vv) .<= π + 1.0e-8)
        end
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
        # `(params, hyperparams)` samples pass through PosteriorSamples, `rmap`, `caltable`,
        # and `siteparams`
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
        # combining mean and std: siteparams unwraps both before the elementwise op
        combined = map(
            (x, y) -> Comrade.siteparams(x) .+ Comrade.siteparams(y),
            mchain.instrument, schain.instrument,
        )
        @test instrumentmodel(post, (; instrument = combined)) isa Comrade.SiteArray
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

        # an IID override on the reference site has every point fixed, so its chain term
        # is an empty sum
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
        # Under whitening the flat prior is standard normal in z: changing only z shifts
        # the flat log density by −Δ‖z‖²/2. This checks the coloring map and its Jacobian.
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
        # IID overrides use their own transform, which keeps bounded supports bounded
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
        # Gains simulated from a known OU process: the hyperparameter posterior conditional
        # on the true gains peaks near the truth. A joint MAP over gains and hyperparameters
        # is degenerate for a hierarchical model, so it is not a valid check.
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
        # FixedInit(0.0) fixes every site's chain to zero at its first time, removing the
        # level degeneracy with a separate circular offset term
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

        # ascube works for FixedInit chains (checked without the circular offset term,
        # which refuses Std transports)
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

        # FixedInit + centered + fitted hyperparameters
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

        # a reference value that disagrees with FixedInit at the first point is an error
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
        # a circular phase prior; UniformInit covers the per-track offset, so no separate
        # offset term is needed
        @instrument function gmint_wb()
            return @jones begin
                gp ~ ArrayPrior(
                    GaussMarkovSitePrior(ScanSeg(), WrappedBrownian(τ = InverseGamma(3.0, 6.0)); init = UniformInit());
                    refant = SEFDReference(0.0)
                )
                return SingleStokesGain(exp(1im * gp))
            end
        end
        post = VLBIPosterior(skym, gmint_wb(), dvis)
        s = prior_sample(rng, post)
        @test all(h -> keys(h) == (:τ,), values(s.instrument.gp.hyperparams))
        @test isfinite(logdensityof(post, s))

        # the model is invariant under a 2π shift of any free phase
        freeinds = setdiff(1:length(post.prior.instrument.gp.dists), post.prior.instrument.gp.dists.fixedinds)
        s2 = deepcopy(s)
        parent(s2.instrument.gp.params)[first(freeinds)] += 2π
        @test logdensityof(post, s2) ≈ logdensityof(post, s)

        tp = asflat(post)
        # one coordinate per free phase, plus one τ per site and the sky's dimensions
        dchain = post.prior.instrument.gp.dists
        nfreew = length(dchain) - length(dchain.fixedinds)
        nsites = length(keys(s.instrument.gp.hyperparams))
        skydim = LogDensityProblems.dimension(asflat(VLBIPosterior(skym, dvis)))
        @test LogDensityProblems.dimension(tp) == skydim + nsites + nfreew

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
                        GaussMarkovSitePrior(ScanSeg(), WrappedBrownian(τ = InverseGamma(3.0, 6.0)); init = UniformInit(), param = AngleEmbedded());
                        refant = SEFDReference(0.0)
                    )
                    return SingleStokesGain(exp(1im * gp))
                end
            end
            poste = VLBIPosterior(skym, gmint_wbe(), dvis)
            se = prior_sample(rng, poste)
            @test isfinite(logdensityof(poste, se))
            tpe = asflat(poste)
            # each free phase contributes two latent reals, same D and sky dims
            @test LogDensityProblems.dimension(tpe) == skydim + nsites + 2 * nfreew
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
                    GaussMarkovSitePrior(
                        ScanSeg(),
                        WrappedOrnsteinUhlenbeck(σ = Exponential(0.3), τ = InverseGamma(3.0, 6.0))
                    );
                    refant = SEFDReference(0.0)
                )
                return SingleStokesGain(exp(1im * gp))
            end
        end
        post = VLBIPosterior(skym, gmint_wou(), dvis)
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
        nsites = length(keys(s.instrument.gp.hyperparams))
        skydim = LogDensityProblems.dimension(asflat(VLBIPosterior(skym, dvis)))
        # one coordinate per free phase, plus σ and τ per site
        @test LogDensityProblems.dimension(tp) == skydim + 2 * nsites + nfreew

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
                        GaussMarkovSitePrior(
                            ScanSeg(), WrappedOrnsteinUhlenbeck(σ = 0.3, τ = InverseGamma(3.0, 6.0));
                            param = AngleEmbedded()
                        );
                        refant = SEFDReference(0.0)
                    )
                    return SingleStokesGain(exp(1im * gp))
                end
            end
            poste = VLBIPosterior(skym, gmint_woue(), dvis)
            se = prior_sample(rng, poste)
            @test isfinite(logdensityof(poste, se))
            tpe = asflat(poste)
            @test LogDensityProblems.dimension(tpe) == skydim + nsites + 2 * nfreew
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

    @testset "batched chain groups reproduce the per-chain walkers" begin
        # The batched walkers over `chaingroups(d)` match the per-chain walkers for the
        # chain log density, flat coloring and its log-Jacobian, flat whitening, and both
        # Std transports.
        perchain(dd, xx, hh) = sum(Comrade.chain_term(sp, xx, hh) for sp in values(Comrade.chainspecs(dd)))

        function mkdist(procs, init, par, n; fixed, hpnames)
            K = length(procs)
            chains = NamedTuple{ntuple(i -> Symbol(:c, i), K)}(
                ntuple(K) do i
                    Comrade.MarkovChainSpec(
                        procs[i], init, Val(hpnames === nothing ? nothing : hpnames[i]),
                        collect(((i - 1) * n + 1):(i * n)),
                        collect(range(0.1 * i; step = 0.1, length = n)), fixed[i], par
                    )
                end
            )
            fi = Int[]
            fv = Float64[]
            for (i, sp) in enumerate(values(chains))
                append!(fi, sp.fsub.inds)
                append!(fv, fill(0.05 * i, length(sp.fsub.inds)))
            end
            return Comrade.GaussMarkovChainDist(chains, fi, fv, K * n)
        end

        rngb = Random.Xoshiro(5)
        K, n = 4, 6
        rep(p) = ntuple(_ -> p, K)
        allfix(f) = ntuple(_ -> f, K)
        OU = Comrade.OrnsteinUhlenbeck(σ = 0.2, τ = 0.5)
        OUh = Comrade.OrnsteinUhlenbeck(σ = InverseGamma(3.0, 1.0), τ = InverseGamma(3.0, 1.0))
        WB = Comrade.WrappedBrownian(τ = 0.5)
        WBh = Comrade.WrappedBrownian(τ = InverseGamma(3.0, 1.0))
        WOU = Comrade.WrappedOrnsteinUhlenbeck(σ = 0.3, τ = 0.7)
        BM = Comrade.BrownianMotion(D = 0.3)
        sites = (:a, :b, :c, :d)

        cases = (
            ("OU whitened", rep(OU), StationaryInit(), NonCentered(), allfix(Int[]), nothing),
            ("OU centered", rep(OU), StationaryInit(), Centered(), allfix(Int[]), nothing),
            ("OU refant", rep(OU), StationaryInit(), NonCentered(), allfix([2, 5]), nothing),
            ("OU fitted", rep(OUh), StationaryInit(), NonCentered(), allfix(Int[]), sites),
            ("OU fitted refant", rep(OUh), StationaryInit(), Centered(), allfix([1, 4]), sites),
            ("BM GaussianInit", rep(BM), GaussianInit(0.1, 2.0), NonCentered(), allfix(Int[]), nothing),
            ("BM FixedInit refant", rep(BM), FixedInit(0.0), NonCentered(), allfix([1]), nothing),
            ("BM GaussianInit refant", rep(BM), GaussianInit(0.1, 2.0), NonCentered(), allfix([2, 5]), nothing),
            ("WB centered", rep(WB), UniformInit(), Centered(), allfix(Int[]), nothing),
            ("WB centered refant", rep(WB), UniformInit(), Centered(), allfix([2, 5]), nothing),
            ("WB noncentered", rep(WB), UniformInit(), NonCentered(), allfix(Int[]), nothing),
            ("WB noncentered refant", rep(WB), UniformInit(), NonCentered(), allfix([2, 5]), nothing),
            ("WB embedding", rep(WB), UniformInit(), AngleEmbedded(), allfix([3]), nothing),
            ("WB fitted", rep(WBh), UniformInit(), Centered(), allfix([2]), sites),
            ("WB fitted noncentered", rep(WBh), UniformInit(), NonCentered(), allfix([2]), sites),
            ("WOU stationary", rep(WOU), StationaryInit(), Centered(), allfix(Int[]), nothing),
            ("WOU refant", rep(WOU), StationaryInit(), Centered(), allfix([1, 4]), nothing),
            # chains with different numeric σ or τ share one group; each reads its own values
            # from `numvals`
            (
                "mixed numeric fields",
                (
                    OU, Comrade.OrnsteinUhlenbeck(σ = 0.7, τ = 0.5),
                    Comrade.OrnsteinUhlenbeck(σ = 0.2, τ = 1.9), OU,
                ),
                StationaryInit(), NonCentered(), allfix(Int[]), nothing,
            ),
            # chains with and without reference-fixed points share a group: with `mskf = 0`,
            # `_bridge_weights(Q₁, 0, 1)` reduces to the one-sided conditional
            ("mixed refant", rep(OU), StationaryInit(), NonCentered(), ([2], Int[], [1, 5], Int[]), nothing),
            ("mixed refant wrapped", rep(WB), UniformInit(), Centered(), (Int[], [3], Int[], [2, 6]), nothing),
            ("mixed refant wrapped nc", rep(WB), UniformInit(), NonCentered(), (Int[], [3], Int[], [2, 6]), nothing),
        )

        for (nm, procs, init, par, fixed, hpn) in cases
            d = mkdist(procs, init, par, n; fixed, hpnames = hpn)
            @testset "$nm" begin
                # every chain lands in a single group
                @test length(Comrade.chaingroups(d)) == 1
                @test Comrade._nfree(only(Comrade.chaingroups(d))) ==
                    sum(Comrade._nfree, values(Comrade.chainspecs(d)))

                hp = hpn === nothing ? (;) :
                    NamedTuple{hpn}(
                        ntuple(K) do i
                            procs[i] isa Comrade.WrappedBrownian ? (τ = 0.4 + 0.1i,) :
                            (σ = 0.1 + 0.05i, τ = 0.4 + 0.1i)
                    end
                    )
                y = Comrade.is_wrapped(first(procs)) ? (2π .* rand(rngb, K * n) .- π) :
                    randn(rngb, K * n)
                for (i, sp) in enumerate(values(Comrade.chainspecs(d)))
                    y[sp.fsub.inds] .= 0.05 * i
                end
                @test Comrade._chain_logpdf(d, y, hp) ≈ perchain(d, y, hp) rtol = 1.0e-13

                nflat = sum(Comrade._flat_dim, Comrade.chaingroups(d))
                z = randn(rngb, nflat)
                function color(units, zz)
                    yy = zeros(length(d))
                    Comrade._fill_fixed!(yy, d)
                    av, cv = Comrade._scan_buffers(yy, d)
                    ll, _, _ = Comrade._color_specs_flat!(TV.LogJac(), yy, zz, 1, 1, av, cv, units, hp)
                    Comrade._resolve_scan!(yy, d, av, cv)
                    return yy, ll
                end
                yb, ℓb = color(Comrade.chaingroups(d), z)
                yr, ℓr = color(values(Comrade.chainspecs(d)), z)
                @test yb ≈ yr rtol = 1.0e-13
                @test ℓb ≈ ℓr rtol = 1.0e-13
                @test isfinite(Comrade._chain_logpdf(d, yb, hp))

                wb = zeros(nflat)
                Comrade._whiten_specs_flat!(wb, 1, yb, Comrade.chaingroups(d), hp)
                wr = zeros(nflat)
                Comrade._whiten_specs_flat!(wr, 1, yb, values(Comrade.chainspecs(d)), hp)
                @test wb ≈ wr rtol = 1.0e-13
                @test first(color(Comrade.chaingroups(d), wb)) ≈ yb rtol = 1.0e-12

                if all(Comrade._std_transportable, Comrade.chaingroups(d))
                    sp = Comrade.PT.StdNormal()
                    nstd = sum(u -> Comrade._std_dim(u, sp), Comrade.chaingroups(d))
                    zs = randn(rngb, nstd)
                    function colorstd(units)
                        yy = zeros(length(d))
                        Comrade._fill_fixed!(yy, d)
                        av, cv = Comrade._scan_buffers(yy, d)
                        Comrade._color_specs_std!(yy, zs, 1, 1, av, cv, units, hp, sp)
                        Comrade._resolve_scan!(yy, d, av, cv)
                        return yy
                    end
                    ysg = colorstd(Comrade.chaingroups(d))
                    @test ysg ≈ colorstd(values(Comrade.chainspecs(d))) rtol = 1.0e-13
                    sb = zeros(nstd)
                    Comrade._whiten_specs_std!(sb, 1, ysg, Comrade.chaingroups(d), hp, sp)
                    sr = zeros(nstd)
                    Comrade._whiten_specs_std!(sr, 1, ysg, values(Comrade.chainspecs(d)), hp, sp)
                    @test sb ≈ sr rtol = 1.0e-13
                end
            end
        end

        # Chains that differ in process type, init, or `centered` are not merged, and only
        # adjacent chains are merged, so the flat coordinate layout is unchanged.
        @testset "grouping boundaries" begin
            mixed = (
                c1 = Comrade.MarkovChainSpec(OU, StationaryInit(), Val(nothing), 1:4, collect(0.1:0.1:0.4), Int[], NonCentered()),
                c2 = Comrade.MarkovChainSpec(OU, StationaryInit(), Val(nothing), 5:8, collect(0.1:0.1:0.4), Int[], Centered()),
                c3 = Comrade.MarkovChainSpec(BM, FixedInit(0.0), Val(nothing), 9:12, collect(0.1:0.1:0.4), [1], NonCentered()),
                c4 = Comrade.MarkovChainSpec(BM, FixedInit(1.0), Val(nothing), 13:16, collect(0.1:0.1:0.4), [1], NonCentered()),
            )
            dm = Comrade.GaussMarkovChainDist(mixed, [9, 13], [0.0, 1.0], 16)
            @test length(Comrade.chaingroups(dm)) == 4
            # same key but non-adjacent: three groups, not two, so the layout is preserved
            order = (a = mixed.c1, b = mixed.c2, c = mixed.c1)
            da = Comrade.GaussMarkovChainDist(order, Int[], Float64[], 16)
            @test length(Comrade.chaingroups(da)) == 3
        end
    end

    @testset "chain times are absolute and on the data's clock" begin
        # `spec.ts` uses the data's `Ti` clock, with one origin shared by every site
        dtab = Comrade.datatable(arrayconfig(dvis))
        function _specs(seg, centered)
            sp = GaussMarkovSitePrior(seg, Comrade.WrappedBrownian(τ = 2.0); init = UniformInit(), centered = centered)
            op = Comrade.ObservedArrayPrior(ArrayPrior(sp), arrayconfig(dvis))
            return Comrade.chainspecs(op.dists)
        end
        for (seg, exact) in ((IntegSeg(), true), (ScanSeg(), false))
            specst = _specs(seg, true)
            for (site, sp) in pairs(specst)
                st = sort(unique(dtab.Ti[findall(b -> site in b, dtab.sites)]))
                # the data start at ~0.58 h, so no chain starts at 0; every chain time lies
                # within the site's observed span
                @test !iszero(first(sp.ts))
                @test all(t -> minimum(st) - 1.0e-8 ≤ t ≤ maximum(st) + 1.0e-8, sp.ts)
                @test issorted(sp.ts) && allunique(sp.ts)
                if exact
                    # one parameter per integration: chain times are datum times
                    @test all(t -> any(u -> isapprox(u, t; atol = 1.0e-8), st), sp.ts)
                    @test first(sp.ts) ≈ first(st)
                end
            end
            # one shared clock: sites observing together agree on a scan's time
            allts = [Set(round.(collect(sp.ts); digits = 8)) for sp in values(specst)]
            @test any(!isempty, [intersect(allts[a], allts[b]) for a in 1:length(allts) for b in (a + 1):length(allts)])
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
