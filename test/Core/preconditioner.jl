using LinearAlgebra
using Serialization
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

@testset "Fisher fit: regularization at N << n" begin
    # The projected covariances are formed as XXᵀ/γ + I, so the shrinkage target is the
    # identity. That is what keeps the estimator usable when far fewer draws than dimensions
    # are available: a direction the window carries no information about lands at eigenvalue
    # 1, where the two-sided filter ignores it, rather than at 0 or ∞ where the filter would
    # admit it as an extreme stiff or wide direction. It also keeps both covariances positive
    # definite, which is what makes the geometric mean of the two well posed.
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
    # Every fitted scale is a geometry estimate, not a numerical zero. Shrinking the
    # covariances toward the identity bounds how extreme a correction a finite window can
    # propose, so the stiff scale stays near the planted 0.1 instead of underflowing.
    @test minimum(pre.s) > 0.03

    @test all(isfinite, pre.s)
end

@testset "Fisher fit: score-informed center along the fitted directions" begin
    # The center carries the mean-score correction on the low-rank part as well as on the
    # diagonal. Dropping the low-rank half leaves the center biased along exactly the
    # directions the fit corrects, and unlike sampling noise that bias does not shrink as
    # draws accumulate. Measured as where a Gaussian's known center lands in the fitted
    # latent space, which is the origin for a perfect fit.
    rng = Random.Xoshiro(3)
    n, N = 400, 320
    u = normalize(randn(rng, n))
    v = normalize(randn(rng, n)); v .-= dot(v, u) * u; normalize!(v)
    Σ = Symmetric(I + 35.0 * u * u' - 0.99 * v * v')
    m = 2.0 .* randn(rng, n)                        # true center, unknown to the fit
    Z = m .+ cholesky(Σ).L * randn(rng, n, N)
    G = -(Σ \ (Z .- m))

    pre = Comrade._fisher_lowrank(Z, G; rank = 64)
    @test norm(Comrade._affine_inv(pre, m)) < 1.0

    # The diagonal-only center, for contrast: the same fit with the low-rank half of the
    # correction removed.
    diagonly = Comrade.LowRankPreconditioner(
        vec(mean(Z; dims = 2)) .+ pre.d .^ 2 .* vec(mean(G; dims = 2)),
        pre.d, pre.V, pre.s
    )
    @test norm(Comrade._affine_inv(pre, m)) < 0.2 * norm(Comrade._affine_inv(diagonly, m))
end

@testset "SPD geometric mean under a graded spectrum" begin
    # `_fisher_lowrank` hands `_spdm` projected covariances whose directions span orders of
    # magnitude, because a real posterior's do, and whose draw and score power grade in
    # opposite senses. Forming the congruence A^½ B A^½ explicitly squares that dynamic
    # range and the small eigenvalues of the result come back negative; the caller's
    # two-sided filter then reads each of those as an infinitely stiff direction and the
    # transform compresses that axis by an unbounded 1/√λ. Routing through the Cholesky
    # factors holds the result positive definite.
    rng = Random.Xoshiro(11)
    m, k, γ = 382, 192, 1.0e-5
    g = exp10.(range(-4, 4, length = k))
    Px = randn(rng, m, k) * Diagonal(g)
    Pa = randn(rng, m, k) * Diagonal(reverse(g))
    Cx = Symmetric(Px * Px' ./ γ + I)
    Ca = Symmetric(Pa * Pa' ./ γ + I)
    Σ = Comrade._spdm(Ca, Cx)
    @test minimum(eigvals(Σ)) > 0
    # Σ is what it is defined to be: the solution of Σ Ca Σ = Cx.
    @test norm(Σ * Ca * Σ - Cx) / norm(Cx) < 1.0e-8
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
        @test_throws "weight γ must be positive" FisherLowRank(; γ = 0.0)
        @test_throws "min_draws must be at least 4" FisherLowRank(; min_draws = 3)
        @test_throws "unknown refit schedule" FisherLowRank(; schedule = :welford)
        @test_throws "must lie strictly in (0, 1)" FisherLowRank(; schedule = [0.5, 1.5])
        @test_throws "must be :stan, :nutpie" FisherLowRank(; schedule = "stan")
        @test_throws "discard must lie in [0, 1)" FisherLowRank(; discard = 1.0)
        @test_throws "discard must lie in [0, 1)" FisherLowRank(; discard = -0.1)
    end

    @testset "refit draw selection" begin
        # without discard only the first 3 draws are dropped
        @test Comrade._refit_selection(FisherLowRank(), 100) == 4:100
        # discard drops the leading fraction of the recorded draws
        @test Comrade._refit_selection(FisherLowRank(; discard = 0.5), 100) == 51:100
        # but never below min_draws
        @test Comrade._refit_selection(FisherLowRank(; discard = 0.9, min_draws = 12), 20) == 9:20
        # and the kept draws are thinned to at most max_fit_draws
        sel = Comrade._refit_selection(FisherLowRank(; discard = 0.5, max_fit_draws = 20), 400)
        @test length(sel) == 20 && first(sel) == 201 && last(sel) == 400
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
        @test minimum(fit.s) > 0.03

        # A carrying refit keeps the current transform's directions (here one the window
        # does not single out, next to a zero-padded device column) and adds the planted ones.
        w = normalize(randn(rng, n)); w .-= dot(w, u) * u .+ dot(w, v) * v; normalize!(w)
        cur = LowRankPreconditioner(zeros(n), ones(n), hcat(w, zeros(n)), [0.01, 1.0])
        cad = FisherLowRank(; rank = 8, min_draws = 12, carry = :rescale)
        cfit = Comrade.metric_refit(cad, st; current = cur)
        @test length(cfit.s) <= 8
        @test maximum(abs.(cfit.V' * normalize(w ./ cfit.d))) > 0.99
        @test maximum(abs.(cfit.V' * normalize(u ./ cfit.d))) > 0.8
        @test maximum(abs.(fit.V' * normalize(w ./ fit.d))) < 0.9
        @test Comrade.metric_refit(ad, st; current = cur).V ≈ fit.V
        # `:keep` carries the starting transform exactly: its direction, scale and diagonal
        kad = FisherLowRank(; rank = 8, min_draws = 12, carry = :keep)
        kcur = LowRankPreconditioner(zeros(n), fill(2.0, n), hcat(w, zeros(n)), [0.01, 1.0])
        kfit = Comrade.metric_refit(kad, st; current = cur, initial = kcur)
        @test kfit.d == kcur.d
        @test kfit.V[:, 1] == w && kfit.s[1] == 0.01
        @test length(kfit.s) <= 8
        @test kfit.V' * kfit.V ≈ I atol = 1.0e-8
        @test maximum(abs.(kfit.V' * normalize(u ./ kfit.d))) > 0.8
        @test Comrade.metric_refit(kad, st; current = kcur).V ≈ fit.V
        @test_throws "carry must be :none, :rescale or :keep" FisherLowRank(; carry = :yes)
    end
end

@testset "low-rank preconditioner in the StdNormal space" begin
    rng = Random.Xoshiro(3)
    gprior = (
        f1 = VLBIGaussian(1.0, 0.1), σ1 = VLBIGaussian(μas2rad(20.0), μas2rad(2.0)),
        τ1 = VLBIGaussian(0.5, 0.05), ξ1 = VLBIGaussian(0.0, 0.3),
        f2 = VLBIGaussian(0.5, 0.1), σ2 = VLBIGaussian(μas2rad(20.0), μas2rad(2.0)),
        τ2 = VLBIGaussian(0.5, 0.05), ξ2 = VLBIGaussian(0.0, 0.3),
        x = VLBIGaussian(0.0, μas2rad(20.0)), y = VLBIGaussian(0.0, μas2rad(20.0)),
    )
    _, vis, _, _, _ = load_data()
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 12, 12)
    post = VLBIPosterior(SkyModel(test_model, gprior, g), vis)
    sp = Comrade.PT.StdNormal()
    tstd = Comrade.transport_to(post, sp)
    n = dimension(tstd)
    V = Matrix(qr(randn(rng, n, 2)).Q)[:, 1:2]
    pre = LowRankPreconditioner(randn(rng, n), exp.(0.3 .* randn(rng, n)), V, [3.0, 0.5])
    tpre = Comrade.transport_to(post, Preconditioned(sp, pre))
    @test dimension(tpre) == n

    # z ↦ w = b + A z, then the StdNormal transport; the density picks up log|det A|
    z = randn(rng, n)
    w = Comrade._affine_fwd(pre, z)
    x = transform(tpre, z)
    @test inverse(tstd, x) ≈ w
    @test inverse(tpre, x) ≈ z
    @test logdensityof(tpre, z) ≈ logdensityof(tstd, w) + Comrade._affine_logdet(pre)

    # the reference draws z = A⁻¹(w − b) with w ~ N(0, I)
    stop = getfield(tpre.transform, :stop)
    @test rand(Random.Xoshiro(5), stop) ≈
        Comrade._affine_inv(pre, rand(Random.Xoshiro(5), getfield(tstd.transform, :stop)))

    @test Comrade._base_space(tpre) === sp
    @test Comrade._base_space(tstd) === sp
    @test isnothing(Comrade._base_space(asflat(post)))
    @test Comrade._transport_pre(tpre) === pre
    @test isnothing(Comrade._transport_pre(tstd))
    @test Comrade._in_space(sp, pre) isa Preconditioned
    @test Comrade._in_space(nothing, pre) === pre
    bad = LowRankPreconditioner(zeros(n + 1), ones(n + 1), zeros(n + 1, 0), Float64[])
    @test_throws DimensionMismatch Comrade.transport_to(post, Preconditioned(sp, bad))

    # score init centers the transform at the start point in the base space; the score comes
    # from the cache (here the standard-normal score -x)
    θ0 = prior_sample(rng, post)
    cache = Comrade._new_refit_cache()
    cache.fn[] = x -> -x
    tm = Comrade._score_init_pre(post, θ0, cache; reactant = false, space = sp)
    @test tm isa Preconditioned && tm.space === sp
    @test tm.pre.b ≈ inverse(tstd, θ0)
    @test Comrade._score_init_pre(post, θ0, cache; reactant = false) isa LowRankPreconditioner

    # a pilot's recorded base-space draws fit only a preconditioner for the same space
    mktempdir() do dir
        draws = [randn(rng, n) for _ in 1:12]
        serialize(joinpath(dir, "metric_adaptation.jls"), Comrade.FisherAdaptation(draws, [-d for d in draws]))
        serialize(joinpath(dir, "transport.jls"), sp)
        @test fit_preconditioner(dir, post; rank = 2, space = sp) isa Preconditioned
        @test_throws "sampled in the StdNormal space" fit_preconditioner(dir, post; rank = 2)
        serialize(joinpath(dir, "transport.jls"), nothing)
        @test fit_preconditioner(dir, post; rank = 2) isa LowRankPreconditioner
        @test_throws "sampled in the flat space" fit_preconditioner(dir, post; rank = 2, space = sp)
    end
end
