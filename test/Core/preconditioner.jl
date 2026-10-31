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

    # Stiff only: the planted stiff direction with the same scale, and no wide direction.
    stiff = Comrade._fisher_lowrank(Z, G; rank = 64, max_wide = 0)
    @test all(<(1), stiff.s)
    @test maximum(abs.(stiff.V' * normalize(v ./ stiff.d))) > 0.8
    @test maximum(abs.(stiff.V' * normalize(u ./ stiff.d))) < 0.5
    @test sort(stiff.s) ≈ sort(filter(<(1), pre.s))
    refit = Comrade.metric_refit(
        FisherLowRank(; rank = 64, max_wide = 0),
        Comrade.FisherAdaptation(collect.(eachcol(Z)), collect.(eachcol(G)))
    )
    @test !isempty(refit.s) && all(<(1), refit.s)

    # One wide direction allowed: the planted wide one, beside the same stiff ones.
    mixed = Comrade._fisher_lowrank(Z, G; rank = 64, max_wide = 1)
    @test count(>(1), mixed.s) == 1
    @test maximum(abs.(mixed.V' * normalize(u ./ mixed.d))) > 0.8
    @test sort(filter(<(1), mixed.s)) ≈ sort(stiff.s)
    @test_throws "max_wide must be nonnegative" FisherLowRank(; max_wide = -1)
end

@testset "Fisher fit: held-out rank selection" begin
    # Stiff directions of variance 1e-4, 1e-3 and 1e-2 dominate the leapfrog stiffness;
    # 0.1 and 0.3 sit near the noise floor of the held-out scores' Gram matrix, so the
    # selected rank stops at the three that set the step.
    rng = Random.Xoshiro(7)
    n, N = 300, 200
    B = Matrix(qr(randn(rng, n, 7)).Q)[:, 1:7]
    Σ = Symmetric(I + B * Diagonal([1.0e-4, 1.0e-3, 1.0e-2, 0.1, 0.3, 30.0, 10.0] .- 1) * B')
    Z = cholesky(Σ).L * randn(rng, n, N)
    G = -(Σ \ Z)
    st = Comrade.FisherAdaptation(collect.(eachcol(Z)), collect.(eachcol(G)))
    a = FisherLowRank(; rank = 64, max_wide = 0, holdout = 0.25, discard = 0.0)
    r = Comrade._holdout_rank(a, st, 1:N)
    @test 3 <= r <= 5
    pre = Comrade.metric_refit(a, st)
    @test length(pre.s) == r && all(<(1), pre.s)
    @test Comrade._holdout_rank(a, Comrade.FisherAdaptation(st.draws[1:20], st.scores[1:20]), 1:20) == 64

    @test_throws "holdout must lie in [0, 1)" FisherLowRank(; holdout = 1.0)
    @test_throws "rank_tol must lie in (0, 1)" FisherLowRank(; holdout = 0.25, rank_tol = 0.0)
    @test_throws "needs carry = :none" FisherLowRank(; holdout = 0.25, carry = :rescale)
end

@testset "nutpie schedule adaptation" begin
    nut, stan = FisherLowRank(; rank = 50, schedule = :nutpie), FisherLowRank(; rank = 50)
    # step-size adaptation carries across nutpie refits but the last; Stan restarts at each
    @test !Comrade.restart_step_size(nut, false) && Comrade.restart_step_size(nut, true)
    @test Comrade.restart_step_size(stan, false) && Comrade.restart_step_size(stan, true)
    @test Comrade.restart_step_size(FixedMetric(), false)
    # the nutpie device block holds the full rank; other schedules grow it with the fits
    @test Comrade.fixed_rank_cap(nut) == 50 && isnothing(Comrade.fixed_rank_cap(stan))
    @test Comrade.fixed_rank_cap(FisherLowRank(; rank = 4, schedule = :nutpie)) == 16
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
        # the window starts after the first 3 records, or after the `discard` fraction
        @test Comrade._window_start(FisherLowRank(; discard = 0.0), 100) == 4
        @test Comrade._window_start(FisherLowRank(), 100) == 51
        # and a fit thins the stored draws evenly to at most max_fit_draws
        @test Comrade._refit_selection(FisherLowRank(), 50) == 1:50
        sel = Comrade._refit_selection(FisherLowRank(; max_fit_draws = 20), 400)
        @test length(sel) == 20 && first(sel) == 1 && last(sel) == 400
    end

    @testset "draw store" begin
        a = FisherLowRank(; max_fit_draws = 12, min_draws = 6)
        observe!(st, i) = Comrade.observe_draw!(a, st, nothing, fill(Float64(i), 3), fill(-Float64(i), 3))
        st = Comrade.FisherAdaptation()
        for i in 1:1000
            observe!(st, i)
        end
        @test st.nrecords == 1000
        # bounded, inside the window, evenly spaced, and the stored values are the records'
        @test Comrade._store_capacity(a) ÷ 2 <= length(st.draws) <= Comrade._store_capacity(a)
        @test first(st.record) >= Comrade._window_start(a, 1000)
        @test st.stride > 1 && all(d -> d <= 2 * st.stride, diff(st.record))
        @test all(k -> st.draws[k] == fill(st.record[k], 3) && st.scores[k] == fill(-st.record[k], 3), eachindex(st.record))
        # a large discard still keeps min_draws
        b = FisherLowRank(; discard = 0.9, min_draws = 6)
        sb = Comrade.FisherAdaptation()
        for i in 1:20
            Comrade.observe_draw!(b, sb, nothing, [Float64(i)], [0.0])
        end
        @test length(sb.draws) == 6

        # The checkpoint replays to the live store at the step it was last saved.
        path = joinpath(mktempdir(), "metric_adaptation.jls")
        ck = Comrade.adaptation_checkpoint(path, a)
        live = Comrade.FisherAdaptation()
        saved = nothing
        for i in 1:995
            observe!(live, i)
            if i % 10 == 0
                Comrade.save_adaptation!(ck, live)
                saved = deepcopy(live)
            end
        end
        loaded = Comrade.resume_adaptation!(a, Comrade.load_adaptation(path, a), 990 * a.record_every)
        @test loaded.record == saved.record && loaded.draws == saved.draws &&
            loaded.scores == saved.scores && loaded.stride == saved.stride &&
            loaded.nrecords == saved.nrecords
        # an interrupted append leaves a partial draw at the end, which is skipped
        open(io -> write(io, Int64(991), 1.0), path * ".log", "a")
        @test_logs (:warn, r"partially written draw") Comrade.load_adaptation(path, a)
        # other adaptors checkpoint their state whole
        @test Comrade.adaptation_checkpoint(path, FixedMetric()) == path
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
            xb = Comrade.observe_draw!(ad, st, pre, z, gz)
            @test xb ≈ x
        end
        # the store keeps the next fit's window
        @test st.nrecords == N
        @test length(st.draws) == N - Comrade._window_start(ad, N) + 1
        @test length(st.scores) == length(st.draws)

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

    @testset "fixed directions" begin
        # The target is wide along w1 and w2 (sd 20 and 8), but the window's draws barely
        # move along them, so a fit from the window alone leaves them narrow.
        rng = Random.Xoshiro(9)
        n, N = 50, 40
        Q = Matrix(qr(randn(rng, n, 2)).Q)[:, 1:2]
        w1, w2 = Q[:, 1], Q[:, 2]
        τ = [20.0, 8.0]
        Σwin = Symmetric(I + Q * Diagonal([0.25, 0.25] .- 1) * Q')
        Σ = Symmetric(I + Q * Diagonal(τ .^ 2 .- 1) * Q')
        Z = cholesky(Σwin).L * randn(rng, n, N)
        G = -(Σ \ Z)
        fst = Comrade.FisherAdaptation(collect(eachcol(Z)), collect(eachcol(G)))
        metric_cov(p) = Diagonal(p.d) * (I + p.V * Diagonal(p.s .^ 2 .- 1) * p.V') * Diagonal(p.d)

        plain = Comrade.metric_refit(FisherLowRank(; rank = 8), fst)
        @test w1' * metric_cov(plain) * w1 < 10

        fit = Comrade.metric_refit(FisherLowRank(; rank = 8, fixed = (reshape(3w1, n, 1), [τ[1]])), fst)
        @test w1' * metric_cov(fit) * w1 ≈ τ[1]^2 rtol = 0.02
        @test abs(dot(fit.V[:, 1], normalize(w1 ./ fit.d))) ≈ 1
        @test fit.V' * fit.V ≈ I atol = 1.0e-8
        @test length(fit.s) <= 8

        two = Comrade.metric_refit(FisherLowRank(; rank = 8, fixed = (hcat(w1, w2), τ)), fst)
        @test Q' * metric_cov(two) * Q ≈ Diagonal(τ .^ 2) rtol = 0.05
        @test two.V' * two.V ≈ I atol = 1.0e-8

        fx = (reshape(w1, n, 1), [τ[1]])
        @test_throws "fixed directions need carry = :none" FisherLowRank(; carry = :rescale, fixed = fx)
        @test_throws DimensionMismatch FisherLowRank(; fixed = (hcat(w1, w2), [1.0]))
        @test_throws "standard deviations must be positive" FisherLowRank(; fixed = (reshape(w1, n, 1), [0.0]))
        @test_throws "between 1 and rank - 1" FisherLowRank(; rank = 2, fixed = (hcat(w1, w2), τ))
        @test_throws "linearly independent" Comrade.metric_refit(FisherLowRank(; fixed = (hcat(w1, 2w1), τ)), fst)
        @test_throws DimensionMismatch Comrade.metric_refit(
            FisherLowRank(; fixed = (reshape(w1[1:10], 10, 1), [1.0])), fst
        )
        @test occursin("1 fixed directions", sprint(show, FisherLowRank(; fixed = fx)))
    end

    @testset "carried stiff directions outnumbering the draws" begin
        # The target's stiffness lies along the carried directions, with the carried scales
        # off by up to a factor 3 and a target diagonal away from the carried unit one.
        rng = Random.Xoshiro(5)
        n, k, N = 200, 60, 20
        dt = exp.(0.3 .* randn(rng, n))
        Vt = Matrix(qr(randn(rng, n, k)).Q)[:, 1:k]
        st = exp.(range(log(1.0e-3), log(10.0), length = k))
        P = Symmetric(Diagonal(inv.(dt)) * (I + Vt * Diagonal(inv.(st .^ 2) .- 1) * Vt') * Diagonal(inv.(dt)))
        Vc = Matrix(qr(Vt ./ dt).Q)[:, 1:k]
        E = eigen(Symmetric(Vc' * P * Vc))
        Vc = Vc * E.vectors
        cur = LowRankPreconditioner(zeros(n), ones(n), Vc, inv.(sqrt.(E.values)) .* exp.(log(1 / 3) .* rand(rng, k)))
        Z = cholesky(P).U \ randn(rng, n, N)
        G = -(P * Z)
        fst = Comrade.FisherAdaptation(collect(eachcol(Z)), collect(eachcol(G)))
        cfit = Comrade.metric_refit(FisherLowRank(; rank = 80, min_draws = 12, carry = :rescale), fst; current = cur)
        @test cfit.d == cur.d
        @test cfit.V[:, 1:k] == cur.V
        # The smallest target-to-metric sd ratio sets the step size; the refit keeps it near
        # the carried transform's rather than re-estimating a diagonal the window cannot support.
        metric_cov(p) = Diagonal(p.d) * (I + p.V * Diagonal(p.s .^ 2 .- 1) * p.V') * Diagonal(p.d)
        worst(p) = minimum(inv.(sqrt.(eigvals(Matrix(P), inv(metric_cov(p))))))
        @test worst(cfit) > 0.5 * worst(cur)
    end

    @testset "stiff-core correction" begin
        @test_throws "directions must lie in 1:steps" StiffCore(; directions = 50, steps = 40)
        @test_throws "level must be positive" StiffCore(; level = 0.0)
        @test_throws "power must be at least 1" StiffCore(; power = 0.5)
        @test_throws "from must lie in [0, 1)" StiffCore(; from = 1.0)
        @test occursin("stiff core", repr(FisherLowRank(; stiff = StiffCore())))
        @test isnothing(Comrade.stiff_core(FisherLowRank()))
        due = FisherLowRank(; schedule = :nutpie, stiff = StiffCore(; from = 0.3))
        @test filter(t -> Comrade.stiff_core_due(due, t, 3000), Comrade.metric_refit_steps(due, 3000, 100)) == [900, 1700, 2500]
        @test !Comrade.stiff_core_due(FisherLowRank(), 2500, 3000)
        # The stiff direction rotates with the draw among four orthonormal directions over a
        # bulk of curvature 3: H(x) = 3I + λ u(x) u(x)ᵀ.
        rng = Random.Xoshiro(11)
        n, N, level, κ = 40, 20, 2.0, 3.0
        Us = Matrix(qr(randn(rng, n, 4)).Q)[:, 1:4]
        λ = 3.0e4
        u(x) = Us[:, mod(round(Int, 1000 * abs(x[1])), 4) + 1]
        H(x) = κ * I + λ * u(x) * u(x)'
        xs = [randn(rng, n) for _ in 1:N]
        st = Comrade.FisherAdaptation(xs, [zeros(n) for _ in 1:N])
        a = FisherLowRank(; stiff = StiffCore(; draws = 20, directions = 2, steps = 10, level))
        pre = LowRankPreconditioner(zeros(n), ones(n), zeros(n, 0), Float64[])
        fixed = Comrade.stiff_core_refit(a, st, pre, (x, Y) -> H(x) * Y)
        A = I + fixed.V * Diagonal(fixed.s .- 1) * fixed.V'
        worst = maximum(x -> eigmax(Symmetric(A * H(x) * A)), xs)
        # A direction stiff at a fraction f of the draws has power-16 mean ≈ (κ + λ) f^(1/16),
        # so the correction leaves the least frequent direction at level κ f^(-1/16).
        f = minimum(i -> count(x -> u(x) == Us[:, i], xs), 1:4) / N
        @test worst ≈ level * κ * f^(-1 / 16) rtol = 1.0e-3
        @test length(fixed.s) == 4
        v = nullspace(Us')[:, 1]
        @test A * v ≈ v
    end
end

@testset "Gauss–Newton low-rank adaptor" begin
    rng = Random.Xoshiro(17)
    n = 40
    lowrank(λ) = (Q = Matrix(qr(randn(rng, n, length(λ))).Q)[:, eachindex(λ)]; (Q * Diagonal(λ) * Q', Q))
    sp = Comrade.PT.StdNormal()

    @testset "construction is validated" begin
        H(x, W) = W
        @test_throws "rank must be positive" GaussNewtonLowRank(H; rank = 0)
        @test_throws "oversample must be at least 1" GaussNewtonLowRank(H; rank = 4, oversample = 0)
        @test_throws "probes_per_draw must lie in 1:14" GaussNewtonLowRank(H; rank = 4, probes_per_draw = 15)
        @test_throws "threshold must be positive" GaussNewtonLowRank(H; rank = 4, threshold = 0)
        @test_throws "min_draws must be at least 1" GaussNewtonLowRank(H; rank = 4, min_draws = 0)
        @test_throws "refit schedule" GaussNewtonLowRank(H; rank = 4, schedule = :often)
        a = GaussNewtonLowRank(H; rank = 4)
        @test isnothing(Comrade.check_metric_space(a, sp))
        @test_throws "GaussNewtonLowRank needs the StdNormal latent space" Comrade.check_metric_space(a, nothing)
        @test Comrade.metric_refit_steps(a, 1000, 10) == Comrade.metric_refit_steps(FisherLowRank(), 1000, 10)
    end

    @testset "constant curvature: exact eigenpairs" begin
        λ0 = [5.0e4, 900.0, 150.0, 20.0, 0.5]
        A, Q = lowrank(λ0)
        a = GaussNewtonLowRank((x, W) -> A * W; rank = 4, oversample = 3, threshold = 100.0, min_draws = 2)
        st = Comrade.init_metric_adaptation(a)
        @test isnothing(Comrade.metric_refit(a, st))
        # draws arrive in the sampler's coordinates and are recorded in the base space
        pre = LowRankPreconditioner(zeros(n), fill(2.0, n), reshape(normalize(randn(rng, n)), n, 1), [3.0])
        for _ in 1:2
            z = randn(rng, n)
            @test Comrade.observe_draw!(a, st, pre, z, zero(z)) ≈ Comrade._affine_fwd(pre, z)
        end
        fit = Comrade.metric_refit(a, st)
        @test fit isa LowRankPreconditioner
        @test fit.b == zeros(n) && fit.d == ones(n)
        @test fit.s ≈ inv.(sqrt.(1 .+ λ0[1:3])) rtol = 1.0e-8
        @test abs.(fit.V' * Q[:, 1:3]) ≈ I atol = 1.0e-8
        # the transform whitens the posterior Hessian I + A on its span
        M = I + fit.V * Diagonal(fit.s .- 1) * fit.V'
        @test fit.V' * (M * (I + A) * M) * fit.V ≈ I atol = 1.0e-8
        # the sketch is reset by a fit
        @test st.ndraws == 0 && all(iszero, st.counts) && iszero(st.Y)
        @test isnothing(Comrade.metric_refit(a, st))

        # a complete probe set (rank + oversample = n) of a low-rank curvature is exact too
        full = GaussNewtonLowRank((x, W) -> A * W; rank = n - 4, oversample = 4, threshold = 100.0, min_draws = 1)
        sf = Comrade.init_metric_adaptation(full)
        Comrade.observe_draw!(full, sf, nothing, randn(rng, n), zeros(n))
        @test Comrade.metric_refit(full, sf).s ≈ inv.(sqrt.(1 .+ λ0[1:3])) rtol = 1.0e-8
    end

    @testset "varying curvature: the sketch averages the draws" begin
        As = [lowrank([400.0 * k, 30.0])[1] for k in 1:3]
        xs = [randn(rng, n) for _ in As]
        curv(x, W) = As[findfirst(==(x), xs)] * W
        a = GaussNewtonLowRank(curv; rank = 8, oversample = 2, threshold = 1.0, min_draws = 3)
        st = Comrade.init_metric_adaptation(a)
        foreach(x -> Comrade.observe_draw!(a, st, nothing, x, zero(x)), xs)
        Ω = Comrade._probe_matrix(a, n)
        # Ω is not stored: its columns are generated the same way each time
        @test Ω == Comrade._probe_matrix(a, n)
        @test Comrade._probe_columns(a, n, [3, 1]) == Ω[:, [3, 1]]
        @test all(c -> 0.5 < norm(c) < 1.5, eachcol(Ω))
        @test st.Y ./ st.counts' ≈ sum(As) / 3 * Ω
        E = eigen(Symmetric(sum(As) / 3); sortby = -)
        fit = Comrade.metric_refit(a, st)
        @test fit.s ≈ inv.(sqrt.(1 .+ filter(>=(1.0), E.values))) rtol = 1.0e-6

        # probes cycle through the columns; a refit waits until each has been applied
        b = GaussNewtonLowRank(curv; rank = 8, oversample = 2, probes_per_draw = 4, threshold = 1.0, min_draws = 1)
        sb = Comrade.init_metric_adaptation(b)
        Comrade.observe_draw!(b, sb, nothing, xs[1], zero(xs[1]))
        Comrade.observe_draw!(b, sb, nothing, xs[2], zero(xs[2]))
        @test sb.counts == [1, 1, 1, 1, 1, 1, 1, 1, 0, 0] && sb.next == 9
        @test isnothing(Comrade.metric_refit(b, sb))
        Comrade.observe_draw!(b, sb, nothing, xs[3], zero(xs[3]))
        @test sb.counts == [2, 2, 1, 1, 1, 1, 1, 1, 1, 1] && sb.next == 3
        Ωb = Comrade._probe_matrix(b, n)
        @test sb.Y[:, 1] ≈ (As[1] + As[3]) * Ωb[:, 1]
        @test sb.Y[:, 5] ≈ As[2] * Ωb[:, 5]

        # products formed a few columns at a time add up to the same sketch
        s1, s3 = Comrade.init_metric_adaptation(a), Comrade.init_metric_adaptation(a)
        foreach(st -> (st.Y = zeros(n, 10); st.counts = zeros(Int, 10)), (s1, s3))
        Comrade._apply_probes!(a, s1, xs[2], 1:10)
        Comrade._apply_probes!(a, s3, xs[2], 1:10; block = 3)
        @test s3.Y ≈ s1.Y rtol = 1.0e-12
        @test s1.Y ≈ As[2] * Comrade._probe_matrix(a, n) rtol = 1.0e-12

        # the sketch round-trips through serialization and keeps accumulating
        path = tempname()
        serialize(path, st)
        foreach(x -> Comrade.observe_draw!(a, st, nothing, x, zero(x)), xs)
        st2 = deserialize(path)
        rm(path)
        foreach(x -> Comrade.observe_draw!(a, st2, nothing, x, zero(x)), xs)
        @test st2.Y == st.Y && st2.counts == st.counts
    end

    @testset "directions confined to rows" begin
        # RowSupportedMatrix behaves as the dense matrix it stands for
        rows = sort(randperm(rng, n)[1:25])
        M = Matrix(qr(randn(rng, 25, 4)).Q)[:, 1:4]
        R = RowSupportedMatrix(n, rows, M)
        D = Array(R)
        @test size(R) == (n, 4) && D[rows, :] == M && iszero(D[setdiff(1:n, rows), :])
        @test R[rows[3], 2] == M[3, 2] && R[first(setdiff(1:n, rows)), 2] == 0
        z, c = randn(rng, n), randn(rng, 4)
        @test R' * z ≈ D' * z && R * c ≈ D * c
        @test R' * [z z] ≈ D' * [z z] && R * [c c] ≈ D * [c c]
        @test Array(R[:, [2, 4]]) == D[:, [2, 4]]
        @test Comrade.Adapt.parent_type(typeof(R)) === typeof(R)
        @test_throws "strictly increasing" RowSupportedMatrix(n, [3, 2], M[1:2, :])
        @test_throws DimensionMismatch RowSupportedMatrix(n, rows[1:3], M)
        # a preconditioner on a RowSupportedMatrix is the same map as on its dense matrix
        s = [0.1, 0.5, 2.0, 7.0]
        pr = LowRankPreconditioner(randn(rng, n), exp.(randn(rng, n)), R, s)
        pd = LowRankPreconditioner(pr.b, pr.d, D, s)
        @test Comrade._affine_fwd(pr, z) ≈ Comrade._affine_fwd(pd, z)
        @test Comrade._affine_inv(pr, z) ≈ Comrade._affine_inv(pd, z)
        @test Comrade._affine_invT(pr, z) ≈ Comrade._affine_invT(pd, z)
        @test Comrade._affine_inv(pr, Comrade._affine_fwd(pr, z)) ≈ z
        @test Comrade._hostify(pr) === pr
        @test_throws "orthonormal columns" LowRankPreconditioner(pr.b, pr.d, RowSupportedMatrix(n, rows, 2 .* M), s)

        # the Gauss–Newton fit on rows is the eigendecomposition of the curvature restricted to them
        λ0 = [5.0e4, 900.0, 150.0, 20.0, 0.5]
        A, _ = lowrank(λ0)
        g = GaussNewtonLowRank((x, W) -> A * W; rank = 6, oversample = 4, threshold = 100.0, min_draws = 1, rows)
        st = Comrade.init_metric_adaptation(g)
        Comrade.observe_draw!(g, st, nothing, randn(rng, n), zeros(n))
        @test size(st.Y) == (length(rows), 10) && st.n == n
        fit = Comrade.metric_refit(g, st)
        Er = eigen(Symmetric(A[rows, rows]); sortby = -)
        keep = findall(>=(100.0), Er.values)
        @test fit.V isa RowSupportedMatrix && fit.V.rows == rows && size(fit.V) == (n, length(keep))
        @test fit.s ≈ inv.(sqrt.(1 .+ Er.values[keep])) rtol = 1.0e-8
        @test abs.(fit.V.M' * Er.vectors[:, keep]) ≈ I atol = 1.0e-8
        @test_throws "strictly increasing" GaussNewtonLowRank((x, W) -> W; rank = 2, rows = [4, 2])
        @test_throws "exceeds the 5 latent coordinates" Comrade.observe_draw!(
            GaussNewtonLowRank((x, W) -> W; rank = 4, oversample = 2, rows = rows[1:5]),
            Comrade.GaussNewtonSketch(), nothing, randn(rng, n), zeros(n)
        )
    end

    @testset "a basis built from the draw" begin
        # a curvature whose range is the span of the first coordinates the draw marks, and a
        # basis that spans those coordinates plus a few more
        support(x) = findall(>(0.5), x)
        function curv(x, W)
            S = support(x)
            B = zeros(n, length(S))
            for (j, i) in enumerate(S)
                B[i, j] = 1.0
            end
            G = B * Diagonal(range(5.0e3, 50.0; length = length(S))) * B'
            return G * W
        end
        basis(x) = Matrix(qr(hcat([Float64.(1:n .== i) for i in support(x)]..., randn(Random.Xoshiro(9), n, 3))).Q)
        a = GaussNewtonLowRank(curv; rank = n, threshold = 100.0, min_draws = 2, basis)
        st = Comrade.init_metric_adaptation(a)
        x0 = vcat(ones(6), zeros(n - 6))
        Comrade.observe_draw!(a, st, nothing, x0, zeros(n))
        @test st.anchor == x0 && size(st.Y) == (n, 9)
        # the window keeps its anchor's basis; the products are taken at the new draw
        Comrade.observe_draw!(a, st, nothing, x0, zeros(n))
        @test st.counts == fill(2, 9) && st.ndraws == 2
        fit = Comrade.metric_refit(a, st)
        λ0 = collect(range(5.0e3, 50.0; length = 6))
        keep = findall(>=(100.0), λ0)
        @test size(fit.V, 2) == length(keep)
        @test sort(fit.s) ≈ sort(inv.(sqrt.(1 .+ λ0[keep]))) rtol = 1.0e-8
        @test isempty(st.anchor) && st.ndraws == 0
        # the next window anchors at its first draw, with a basis of a different size
        x1 = vcat(ones(3), zeros(n - 3))
        Comrade.observe_draw!(a, st, nothing, x1, zeros(n))
        @test st.anchor == x1 && size(st.Y) == (n, 6)
        @test_throws "rows and basis are exclusive" GaussNewtonLowRank(curv; rank = 2, rows = [1, 2], basis)
        @test_throws DimensionMismatch Comrade.observe_draw!(
            GaussNewtonLowRank(curv; rank = 2, basis = x -> zeros(n - 1, 2)), Comrade.GaussNewtonSketch(), nothing, x0, zeros(n)
        )
    end

    @testset "more eigenvalues than rank keeps the largest" begin
        A, _ = lowrank([1.0e4, 5.0e3, 2.0e3, 1.0e3])
        a = GaussNewtonLowRank((x, W) -> A * W; rank = 2, threshold = 10.0, min_draws = 1)
        st = Comrade.init_metric_adaptation(a)
        Comrade.observe_draw!(a, st, nothing, randn(rng, n), zeros(n))
        fit = @test_logs (:info, r"4 eigenvalues exceed threshold = 10.0; keeping the largest rank = 2, the largest dropped is λ = 2000") Comrade.metric_refit(a, st)
        E = eigen(Symmetric(A); sortby = -)
        @test size(fit.V) == (n, 2)
        @test fit.s ≈ inv.(sqrt.(1 .+ [1.0e4, 5.0e3])) rtol = 1.0e-8
        @test abs.(fit.V' * E.vectors[:, 1:2]) ≈ I atol = 1.0e-8
    end

    @testset "residual refits" begin
        H(x, W) = W
        @test_throws "residual_draws must be positive" GaussNewtonLowRank(H; rank = 2, residual_draws = 0)
        @test_throws "residual_at must hold fractions of warmup in (0, 1)" GaussNewtonLowRank(H; rank = 2, residual_at = [1.0])
        @test_throws "residual_directions must lie in 1:residual_steps = 1:5" GaussNewtonLowRank(
            H; rank = 2, residual_steps = 5, residual_directions = 6
        )
        @test_throws "residual_threshold must exceed 1" GaussNewtonLowRank(H; rank = 2, residual_threshold = 1.0)
        @test occursin("residual_at = [0.5]", repr(GaussNewtonLowRank(H; rank = 2, residual_at = [0.5])))
        @test !occursin("residual_at", repr(GaussNewtonLowRank(H; rank = 2)))

        # residual refits run at their own steps, which must not coincide with fits
        na, chunk = 100, 10
        steps = (; n_adapts = na, chunk)
        @test Comrade.metric_refit_steps(GaussNewtonLowRank(H; rank = 2, schedule = [0.2], residual_at = [0.5]), na, chunk) == [20, 50]
        @test_throws "residual_at and schedule share the warmup steps [50]" Comrade.metric_refit_steps(
            GaussNewtonLowRank(H; rank = 2, schedule = [0.2, 0.5], residual_at = [0.5]), na, chunk
        )

        # the kept draws stay evenly spread over the window, between residual_draws and twice that
        a = GaussNewtonLowRank(H; rank = 2, oversample = 2, min_draws = 1, schedule = [0.2], residual_at = [0.5], residual_draws = 3)
        st = Comrade.init_metric_adaptation(a)
        xs = [fill(Float64(i), n) for i in 1:13]
        foreach(x -> Comrade.observe_draw!(a, st, nothing, x, zero(x)), xs)
        @test st.stride == 4 && first.(st.kept) == [1.0, 5.0, 9.0, 13.0]
        @test_throws "needs the refit's step, n_adapts and chunk" Comrade.metric_refit(a, st)
        Comrade.metric_refit(a, st; step = 20, steps...)
        @test isempty(st.kept) && st.stride == 1
        @test isempty(Comrade.init_metric_adaptation(GaussNewtonLowRank(H; rank = 2)).kept)

        whitened(fit, A) = (L = I + fit.V * Diagonal(fit.s .- 1) * fit.V'; Symmetric(L * (I + A) * L))
        # a Nyström fit at step 20 from `fitdraws`, then a residual refit at step 50 from `resdraws`
        function fit_then_residual(g, fitdraws, resdraws)
            sg = Comrade.init_metric_adaptation(g)
            foreach(x -> Comrade.observe_draw!(g, sg, nothing, x, zero(x)), fitdraws)
            fit = @test_logs min_level = Base.CoreLogging.Warn Comrade.metric_refit(g, sg; step = 20, steps...)
            foreach(x -> Comrade.observe_draw!(g, sg, fit, Comrade._affine_inv(fit, x), zero(x)), resdraws)
            return fit, sg
        end

        # directions a rank-truncated fit drops are whitened by a residual refit
        A, _ = lowrank([1.0e4, 5.0e3, 2.0e3])
        g = GaussNewtonLowRank(
            (x, W) -> A * W; rank = 1, threshold = 10.0, min_draws = 1, schedule = [0.2],
            residual_at = [0.5], residual_draws = 2, residual_steps = 10, residual_directions = 3
        )
        fit, sg = fit_then_residual(g, [randn(rng, n)], [randn(rng, n) for _ in 1:2])
        res = @test_logs (:info, r"residual refit: largest curvature left at the draws .*; \d+ directions from 2 draws.*whitening 2 ") Comrade.metric_refit(
            g, sg; current = fit, step = 50, steps...
        )
        @test opnorm(whitened(fit, A)) > 4.0e3
        @test opnorm(whitened(res, A)) ≈ 1 atol = 1.0e-6
        @test res.V' * res.V ≈ I atol = 1.0e-10
        @test size(res.V, 2) == 3
        @test isempty(sg.kept)
        @test isnothing(Comrade.metric_refit(g, sg; current = res, step = 50, steps...))   # no draws kept since

        # a stiff direction that turns slightly from draw to draw, within a few directions,
        # leaves curvature outside the averaged span at draws the fit did not see; a residual
        # refit at such draws takes most of it out at others
        q = normalize(randn(rng, n))
        P = Matrix(qr((I - q * q') * randn(rng, n, 2)).Q)[:, 1:2]
        stiff() = (v = normalize(q + 0.02 * normalize(P * randn(rng, 2))); 1.0e6 * v * v')
        dr = [randn(rng, n) for _ in 1:9]
        Ad = [stiff() for _ in dr]
        # draws come back through the preconditioner, so match them to the nearest one
        curv(x, W) = Ad[argmin(map(d -> norm(d - x), dr))] * W
        g = GaussNewtonLowRank(curv; rank = 1, oversample = 4, threshold = 10.0, min_draws = 3, schedule = [0.2], residual_at = [0.5])
        fit, sg = fit_then_residual(g, dr[1:3], dr[4:6])
        res = @test_logs min_level = Base.CoreLogging.Warn Comrade.metric_refit(g, sg; current = fit, step = 50, steps...)
        worst(p) = maximum(Ak -> opnorm(whitened(p, Ak)), Ad[7:9])
        @test worst(fit) > 100
        @test worst(res) < worst(fit) / 10

        # with rows, the added directions stay on them
        rows = sort(randperm(rng, n)[1:25])
        Ar, _ = lowrank([1.0e4, 5.0e3, 2.0e3])
        g = GaussNewtonLowRank(
            (x, W) -> Ar * W; rank = 1, oversample = 4, threshold = 10.0, min_draws = 1, rows,
            schedule = [0.2], residual_at = [0.5], residual_draws = 1, residual_steps = 10, residual_directions = 3
        )
        fit, sg = fit_then_residual(g, [randn(rng, n)], [randn(rng, n)])
        res = @test_logs min_level = Base.CoreLogging.Warn Comrade.metric_refit(g, sg; current = fit, step = 50, steps...)
        @test res.V isa RowSupportedMatrix && res.V.rows == rows
        L = I + res.V.M * Diagonal(res.s .- 1) * res.V.M'
        @test opnorm(Symmetric(L * (I + Ar[rows, rows]) * L)) ≈ 1 atol = 1.0e-6

        # a residual refit corrects a current preconditioner with zero shift and unit diagonal
        g = GaussNewtonLowRank(H; rank = 2, min_draws = 1, schedule = [0.2], residual_at = [0.5])
        sg = Comrade.init_metric_adaptation(g)
        Comrade.observe_draw!(g, sg, nothing, randn(rng, n), zeros(n))
        @test_throws "the chain samples without one" Comrade.metric_refit(g, sg; step = 50, steps...)
        scaled = LowRankPreconditioner(zeros(n), fill(2.0, n), zeros(n, 0), Float64[])
        @test_throws "zero shift and unit diagonal" Comrade.metric_refit(g, sg; current = scaled, step = 50, steps...)
    end

    @testset "rescaled directions" begin
        H(x, W) = W
        @test_throws "rescale must be nothing or a NutpieDiagonal" GaussNewtonLowRank(H; rank = 2, rescale = 1)
        @test_throws "needs min_draws ≥ 2" GaussNewtonLowRank(H; rank = 2, min_draws = 1, rescale = NutpieDiagonal())
        @test occursin("rescale = ", repr(GaussNewtonLowRank(H; rank = 2, rescale = NutpieDiagonal())))

        # a Gaussian posterior 4× wider than its curvature along the stiffest direction: the
        # curvature scale is too narrow there, and the draws and gradients recover the width
        A, Qa = lowrank([1.0e4; 1.0e2; ones(n - 2)])
        λ, Q = [1.0e4, 1.0e2], Qa[:, 1:2]
        widen = [16.0, 1.0]
        C = I + Q * Diagonal(widen ./ (1 .+ λ) .- 1) * Q'
        L = cholesky(Symmetric(C)).L
        xs = [L * randn(rng, n) for _ in 1:1500]
        fits = map((nothing, NutpieDiagonal(; n0 = 0, bound = 1.0e6))) do r
            g = GaussNewtonLowRank((x, W) -> A * W; rank = 2, oversample = 4, threshold = 10.0, min_draws = 2, residual_draws = 1000, rescale = r)
            sg = Comrade.init_metric_adaptation(g)
            foreach(x -> Comrade.observe_draw!(g, sg, nothing, x, -(C \ x)), xs)
            isnothing(r) ? Comrade.metric_refit(g, sg) :
                (@test_logs (:info, r"rescale from 1500 draws: scale factors over 2 directions") Comrade.metric_refit(g, sg))
        end
        @test first(fits[1].s) ≈ 1 / sqrt(1 + λ[1]) rtol = 1.0e-3
        @test fits[2].s ≈ sqrt.(widen ./ (1 .+ λ)) rtol = 0.05
        @test fits[2].V ≈ fits[1].V

        # the factor is bounded, and with few draws shrunk toward the curvature scale
        C2 = I + Q * Diagonal([1.0e8, 1.0] ./ (1 .+ λ) .- 1) * Q'
        L2 = cholesky(Symmetric(C2)).L
        x2 = [L2 * randn(rng, n) for _ in 1:200]
        factor(nd) = begin
            g = GaussNewtonLowRank((x, W) -> A * W; rank = 2, oversample = 4, threshold = 10.0, min_draws = 2, residual_draws = 1000, rescale = nd)
            sg = Comrade.init_metric_adaptation(g)
            foreach(x -> Comrade.observe_draw!(g, sg, nothing, x, -(C2 \ x)), x2)
            first(Comrade.metric_refit(g, sg).s) * sqrt(1 + λ[1])
        end
        @test 9 < factor(NutpieDiagonal(; n0 = 0, bound = 10)) < 10
        @test factor(NutpieDiagonal(; n0 = 200, bound = 10)) < factor(NutpieDiagonal(; n0 = 0, bound = 10))

        # a residual refit rescales the combined directions too
        g = GaussNewtonLowRank(
            (x, W) -> A * W; rank = 1, oversample = 4, threshold = 10.0, min_draws = 2, schedule = [0.2], residual_at = [0.5],
            residual_draws = 1000, residual_steps = 10, residual_directions = 2, rescale = NutpieDiagonal(; n0 = 0, bound = 1.0e6)
        )
        sg = Comrade.init_metric_adaptation(g)
        foreach(x -> Comrade.observe_draw!(g, sg, nothing, x, -(C \ x)), xs[1:200])
        fit = Comrade.metric_refit(g, sg; step = 20, n_adapts = 100, chunk = 10)
        foreach(x -> Comrade.observe_draw!(g, sg, fit, Comrade._affine_inv(fit, x), Comrade._affine_fwd(fit, -(C \ x))), xs[201:1500])
        res = Comrade.metric_refit(g, sg; current = fit, step = 50, n_adapts = 100, chunk = 10)
        Minv = I + res.V * Diagonal(inv.(res.s) .- 1) * res.V'
        @test opnorm(Symmetric(Minv * C * Minv) - I) < 0.1
    end

    @testset "failures are errors" begin
        bad = GaussNewtonLowRank((x, W) -> fill(NaN, size(W)); rank = 2, min_draws = 1)
        @test_throws "curvature product at a warmup draw is not finite" Comrade.observe_draw!(
            bad, Comrade.init_metric_adaptation(bad), nothing, randn(rng, n), zeros(n)
        )
        short = GaussNewtonLowRank((x, W) -> W[1:3, :]; rank = 2, min_draws = 1)
        @test_throws DimensionMismatch Comrade.observe_draw!(
            short, Comrade.init_metric_adaptation(short), nothing, randn(rng, n), zeros(n)
        )
        @test_throws "exceeds the $n latent coordinates" Comrade._probe_matrix(GaussNewtonLowRank((x, W) -> W; rank = n), n)
        indefinite = GaussNewtonLowRank((x, W) -> -W; rank = 2, min_draws = 1)
        si = Comrade.init_metric_adaptation(indefinite)
        Comrade.observe_draw!(indefinite, si, nothing, randn(rng, n), zeros(n))
        @test_throws "not positive definite on its probes" Comrade.metric_refit(indefinite, si)
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

@testset "NutpieDiagonal" begin
    s = [0.05, 0.5, 1.0, 3.0, 40.0]
    rng = Random.Xoshiro(3)
    d = NutpieDiagonal(; n0 = 0)
    m = Comrade._DiagonalMoments(length(s))
    for _ in 1:4000
        x = s .* randn(rng, length(s))
        Comrade._observe!(m, x, -x ./ s .^ 2)
    end
    # For a diagonal Gaussian the raw scale is exact; the bound maps it through
    # c tanh(log σ̂ / c) and keeps every scale inside [1/bound, bound].
    c = log(10)
    σ = sqrt.(Comrade.inverse_mass(d, m))
    @test σ ≈ exp.(c .* tanh.(log.(s) ./ c)) rtol = 1.0e-12
    @test all(x -> 0.1 < x < 10, σ)
    @test σ[3] ≈ 1
    @test σ[2] ≈ 0.5 rtol = 0.05

    # With few draws the estimate is shrunk toward the unit scale.
    few = Comrade._DiagonalMoments(length(s))
    for _ in 1:3
        x = s .* randn(rng, length(s))
        Comrade._observe!(few, x, -x ./ s .^ 2)
    end
    σ3 = sqrt.(Comrade.inverse_mass(NutpieDiagonal(; n0 = 5), few))
    @test σ3 ≈ exp.(c .* tanh.((3 / 8) .* log.(s) ./ c)) rtol = 1.0e-12

    # A coordinate that never moves keeps the unit scale.
    still = Comrade._DiagonalMoments(2)
    Comrade._observe!(still, [1.0, 0.0], [2.0, 0.0])
    Comrade._observe!(still, [1.0, 1.0], [2.0, -1.0])
    @test Comrade.inverse_mass(d, still)[1] == 1

    Comrade._reset!(m)
    @test m.n == 0 && all(iszero, m.sx)
    @test_throws "at least two draws" Comrade.inverse_mass(d, m)
    @test_throws "bound must exceed 1" NutpieDiagonal(; bound = 1)
    @test_throws "n0 must be nonnegative" NutpieDiagonal(; n0 = -1)
    @test_throws DimensionMismatch Comrade._observe!(m, zeros(4), zeros(5))
end
