using LinearAlgebra
using Random
using Statistics
import TransformVariables as TV

@testset "low-rank preconditioner" begin
    rng = Random.Xoshiro(42)
    n, m = 24, 3
    V = Matrix(qr(randn(rng, n, m)).Q)[:, 1:m]
    b = randn(rng, n)
    d = exp.(0.3 .* randn(rng, n))
    s = [8.0, 4.0, 2.5]
    p = LowRankPreconditioner(b, d, V, s)
    A = Diagonal(d) * (I + V * Diagonal(s .- 1) * V')

    z = randn(rng, n)
    @test Comrade._affine_fwd(p, z) ≈ b .+ A * z
    @test Comrade._affine_inv(p, Comrade._affine_fwd(p, z)) ≈ z
    @test Comrade._affine_logdet(p) ≈ first(logabsdet(A))

    @test_throws ArgumentError LowRankPreconditioner(b, d, randn(rng, n, m), s)
    @test_throws ArgumentError LowRankPreconditioner(b, -d, V, s)
    @test_throws DimensionMismatch LowRankPreconditioner(b[1:3], d, V, s)

    # rank-0: pure diagonal standardization
    p0 = LowRankPreconditioner(b, d, zeros(n, 0), Float64[])
    @test Comrade._affine_fwd(p0, z) ≈ b .+ d .* z
    @test Comrade._affine_logdet(p0) ≈ sum(log, d)

    # TV node over an identity inner transform: x = b + A z, constant log-Jacobian
    t = Comrade.PreconditionedFlat(p, TV.as(Array, n))
    x, ℓ, ix = TV.transform_with(TV.LogJac(), t, z, 1)
    @test x ≈ b .+ A * z
    @test ℓ ≈ first(logabsdet(A))
    @test ix == n + 1
    @test TV.inverse(t, x) ≈ z

    # bounded inner transform: log-Jacobian is the inner's at (b + A z) plus the constant
    tb = Comrade.PreconditionedFlat(p, TV.as(Array, TV.as𝕀, n))
    xb, ℓb, _ = TV.transform_with(TV.LogJac(), tb, z, 1)
    xi, ℓi, _ = TV.transform_with(TV.LogJac(), TV.as(Array, TV.as𝕀, n), b .+ A * z, 1)
    @test xb ≈ xi
    @test ℓb ≈ ℓi + first(logabsdet(A))
    @test TV.inverse(tb, xb) ≈ z

    # Fisher-divergence estimator, N > n: draws plus their exact scores pin the Gaussian
    # down completely, so the fitted transform whitens it.
    nfd = 12
    dd = [100.0, 25.0, 9.0, 1.0e-4, 1.0e-2, 0.04, 1, 1, 1, 1, 1, 1]
    Qf = Matrix(qr(randn(rng, nfd, nfd)).Q)
    Σf = Symmetric(Qf * Diagonal(dd) * Qf')
    Xf = sqrt(Σf) * randn(rng, nfd, 40) .+ randn(rng, nfd)
    Gf = -(Σf \ (Xf .- mean(Xf; dims = 2)))
    pf = Comrade._fisher_lowrank(Xf, Gf; rank = 12, cutoff = 1.3)
    Af = Diagonal(pf.d) * (I + pf.V * Diagonal(pf.s .- 1) * pf.V')
    evf = eigvals(Symmetric(Af \ Matrix(Σf) / Af'))
    @test maximum(evf) / minimum(evf) < 1.5
    zf = randn(rng, nfd)
    @test Comrade._affine_inv(pf, Comrade._affine_fwd(pf, zf)) ≈ zf
end

@testset "Fisher fit: ridge and scale floor at N << n" begin
    # The ridge γ and the stiff-scale floor `s_floor` are the noise control that makes
    # the estimator usable when far fewer draws than dimensions are available. With N
    # draws the joint draw+score subspace has dimension ≤ 2(N-1); sampling noise inside
    # it passes the two-sided eigenvalue filter as spurious wide/stiff directions unless
    # the covariances are ridged, and a genuine stiff direction's fitted scale can round
    # toward zero and drive the sampler's step size with it unless it is floored.
    rng = Random.Xoshiro(1234)
    n, N = 400, 60
    u = normalize(randn(rng, n))                       # planted wide direction
    v = normalize(randn(rng, n)); v .-= dot(v, u) * u; normalize!(v)   # planted stiff one
    Σ = Symmetric(I + 35.0 * u * u' - 0.99 * v * v')   # sd 6 along u, 0.1 along v
    Z = sqrt(Σ) * randn(rng, n, N)
    G = -(Σ \ (Z .- mean(Z; dims = 2)))

    pre = Comrade._fisher_lowrank(Z, G; rank = 64)
    # Both planted directions are found ...
    @test maximum(abs.(pre.V' * normalize(u ./ pre.d))) > 0.8
    @test maximum(abs.(pre.V' * normalize(v ./ pre.d))) > 0.8
    # ... and essentially nothing else: the fit tracks real anisotropy, not noise.
    @test length(pre.s) <= 6
    # No fitted scale sits below the floor.
    @test minimum(pre.s) >= 0.02

    # Weakening the ridge floods the fit with noise directions. This is the regression
    # guard: the ridge is load-bearing, not decoration.
    weak = Comrade._fisher_lowrank(Z, G; rank = 64, γ = 1.0e-9)
    @test length(weak.s) > 4 * length(pre.s)

    # The floor binds from below only; raising it lifts the stiff scales to it.
    lifted = Comrade._fisher_lowrank(Z, G; rank = 64, s_floor = 0.5)
    @test minimum(lifted.s) >= 0.5
end

@testset "metric adaptors" begin
    @test Comrade.adapts_welford(WelfordDiagonal())
    @test !Comrade.adapts_welford(FixedMetric())
    @test !Comrade.adapts_welford(FisherLowRank())
    @test isempty(Comrade.metric_refit_steps(WelfordDiagonal(), 1000, 100))
    @test isnothing(Comrade.init_metric_adaptation(WelfordDiagonal()))

    @testset "construction is validated" begin
        @test_throws "rank must be positive" FisherLowRank(; rank = 0)
        @test_throws "cutoff must be greater than 1" FisherLowRank(; cutoff = 1.0)
        @test_throws "ridge γ must be positive" FisherLowRank(; γ = 0.0)
        @test_throws "s_floor must be positive" FisherLowRank(; s_floor = -1.0)
        @test_throws "min_draws must be at least 4" FisherLowRank(; min_draws = 3)
        @test_throws "unknown refit schedule" FisherLowRank(; schedule = :welford)
        @test_throws "must lie strictly in (0, 1)" FisherLowRank(; schedule = [0.5, 1.5])
        @test_throws "must be :stan, :nutpie" FisherLowRank(; schedule = "stan")
    end

    @testset "refit schedules" begin
        na, chunk = 2000, 100
        stan = Comrade.metric_refit_steps(FisherLowRank(), na, chunk)
        @test first(stan) == 100
        @test issorted(stan) && allunique(stan)
        @test all(s -> 0 < s < na, stan)
        @test last(stan) <= 0.85 * na
        # Gaps grow, so late segments are the long ones dual averaging needs.
        gaps = diff(stan)
        @test issorted(gaps)
        @test length(gaps) < 2 || last(gaps) > first(gaps)

        nutpie = Comrade.metric_refit_steps(FisherLowRank(; schedule = :nutpie), na, chunk)
        @test issorted(nutpie) && allunique(nutpie)
        @test length(nutpie) > length(stan)          # nutpie updates far more often

        manual = Comrade.metric_refit_steps(
            FisherLowRank(; schedule = [0.25, 0.5, 0.75]), na, chunk
        )
        @test manual == [500, 1000, 1500]

        # A warmup too short for the schedule simply gets fewer refits, never an
        # out-of-range step (which the sampler would reject).
        @test all(s -> 0 < s < 120, Comrade.metric_refit_steps(FisherLowRank(), 120, 100))
    end

    @testset "base-flat capture through a transform" begin
        rng = Random.Xoshiro(99)
        n, m = 30, 2
        V = Matrix(qr(randn(rng, n, m)).Q)[:, 1:m]
        pre = LowRankPreconditioner(
            randn(rng, n), exp.(0.3 .* randn(rng, n)), V, [5.0, 0.3]
        )
        A = Diagonal(pre.d) * (I + V * Diagonal(pre.s .- 1) * V')

        # A draw maps forward through the transform, a gradient back through A⁻ᵀ, so a
        # base-flat score recovered from a sampled-space one round-trips exactly.
        z = randn(rng, n)
        @test Comrade._baseflat_draw(pre, z) ≈ A * z .+ pre.b
        gx = randn(rng, n)
        gz = A' * gx
        @test Comrade._baseflat_score(pre, gz) ≈ gx
        @test Comrade._affine_invT(pre, gz) ≈ gx

        # rank-0 and no-transform cases
        p0 = LowRankPreconditioner(pre.b, pre.d, zeros(n, 0), Float64[])
        @test Comrade._affine_invT(p0, gx) ≈ gx ./ pre.d
        @test Comrade._baseflat_draw(nothing, z) == z
        @test Comrade._baseflat_score(nothing, gz) == gz
    end

    @testset "accumulate and refit" begin
        # A Gaussian target with one wide and one stiff direction, observed the way the
        # sampler would report it: positions and gradients in the CURRENT latent space.
        rng = Random.Xoshiro(2024)
        n, N = 60, 40
        u = normalize(randn(rng, n))
        v = normalize(randn(rng, n)); v .-= dot(v, u) * u; normalize!(v)
        Σ = Symmetric(I + 35.0 * u * u' - 0.99 * v * v')
        L = sqrt(Σ)

        pre = LowRankPreconditioner(
            zeros(n), fill(2.0, n), reshape(normalize(randn(rng, n)), n, 1), [3.0]
        )
        A = Diagonal(pre.d) * (I + pre.V * Diagonal(pre.s .- 1) * pre.V')

        ad = FisherLowRank(; rank = 8, min_draws = 12)
        st = Comrade.init_metric_adaptation(ad)
        @test isnothing(Comrade.metric_refit(ad, st))   # nothing recorded yet

        for _ in 1:N
            x = L * randn(rng, n)                        # a base-flat draw
            z = A \ (x .- pre.b)                         # as the sampler sees it
            gz = A' * (-(Σ \ x))                         # and its gradient there
            xb = Comrade.observe_draw!(st, pre, z, gz)
            @test xb ≈ x
        end
        @test length(st.draws) == N
        @test length(st.scores) == N

        fit = Comrade.metric_refit(ad, st)
        @test fit isa LowRankPreconditioner
        @test length(fit.b) == n
        # The planted directions are recovered from sampler-reported quantities alone.
        @test maximum(abs.(fit.V' * normalize(u ./ fit.d))) > 0.8
        @test maximum(abs.(fit.V' * normalize(v ./ fit.d))) > 0.8
        @test minimum(fit.s) >= ad.s_floor
    end
end
