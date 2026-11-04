using LinearAlgebra
using Random
import Distributions as Dists
using Distributions: Exponential
using VLBIFiles

# Toy skies: the flux depends on `a * b` (and `a * c`) only, and the shift on `φ` only through
# `cos φ`, `sin φ`.
function _moves_toy_sky(θ, _, _)
    g = stretched(Gaussian(), μas2rad(20.0), μas2rad(20.0))
    return (θ.a * θ.b * θ.c) * shifted(g, μas2rad(5.0) * cos(θ.φ), μas2rad(5.0) * sin(θ.φ))
end

struct _SheetToy <: Comrade.AbstractMove end
Comrade.step_distribution(::_SheetToy) = SymmetricDiscreteStep(rng -> rand(rng, (-1, 1)), -)
Comrade.move_name(::_SheetToy) = "sheet"
Comrade.likelihood_invariant(::_SheetToy) = true
function Comrade.propose(::_SheetToy, x, s, ctx)
    path = (:sky, :φ)
    x′ = copy(x)
    x′[Comrade.coords(ctx.view, path)] = Comrade.latent(ctx.view, path, Comrade.value(ctx.view, x, path) + 2π * s)
    return x′, 0.0
end

struct _BadStep <: Comrade.AbstractMove end
Comrade.step_distribution(::_BadStep) = :gaussian
Comrade.move_name(::_BadStep) = "badstep"

@testset "moves" begin
    _, vis, _, _, _ = load_data()
    prior = (
        a = Dists.LogNormal(0.0, 0.3), b = Dists.LogNormal(0.0, 0.3),
        c = Dists.Uniform(0.2, 5.0), φ = Dists.Normal(0.0, 3.0),
    )
    g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 8, 8)
    post = VLBIPosterior(vis; sky = SkyModel(_moves_toy_sky, prior, g))
    sp = Comrade.PT.StdNormal()
    rng = Random.Xoshiro(3)
    θs = [prior_sample(rng, post) for _ in 1:3]

    @testset "coordinate view: $(Comrade._space_name(space))" for space in (nothing, sp)
        view = CoordinateView(post, space)
        @test Comrade.space(view) === space
        @test [Comrade.coords(view, (:sky, k)) for k in (:a, :b, :c, :φ)] == [1:1, 2:2, 3:3, 4:4]
        x = Comrade.inverse(view.tbase, θs[1])
        for k in (:a, :b, :c, :φ)
            v = Comrade.value(view, x, (:sky, k))
            @test v ≈ θs[1].sky[k]
            @test Comrade.latent(view, (:sky, k), v) ≈ x[Comrade.coords(view, (:sky, k))]
        end
        @test Comrade.value(view, x, (:sky,)).b ≈ θs[1].sky.b
        @test_throws "no parameter d at (:sky,); found [:a, :b, :c, :φ]" Comrade.coords(view, (:sky, :d))
        @test_throws "is a" Comrade.coords(view, (:sky, :a, :x))
    end
    @test_throws "moves act in the flat or StdNormal latent space" CoordinateView(post, Comrade.PT.StdUniform())

    # b ↦ b a / a′ keeps a b fixed; a log-normal b has an affine latent in both spaces, so
    # the latent map is a shift with unit Jacobian
    ab(view) = CompensatedMove(
        "a_b", view, (:sky, :a), (:sky, :b),
        (vb, x, x′, ctx) -> vb * Comrade.value(ctx.view, x, (:sky, :a)) / Comrade.value(ctx.view, x′, (:sky, :a));
        likelihood_invariant = true
    )
    @testset "compensated move: $(Comrade._space_name(space))" for space in (nothing, sp)
        m = ab(CoordinateView(post, space))
        @test Comrade.step_distribution(m) == NormalStep()
        r = check_move(m, post, θs; space, rng, τ = 0.3, nprior = 400)
        @test r.reversal < 1.0e-12 && r.logdet < 1.0e-6
        @test length(r.ks) == 2 && all(<=(r.ks_critical), r.ks)
    end

    # c ↦ c a / a′ for a uniform c: in the flat (logit) space the Jacobian is not one
    view = CoordinateView(post)
    lo, hi = 0.2, 5.0
    ac(logdet) = CompensatedMove(
        "a_c", view, (:sky, :a), (:sky, :c),
        (vc, x, x′, ctx) -> vc * Comrade.value(ctx.view, x, (:sky, :a)) / Comrade.value(ctx.view, x′, (:sky, :a));
        logdet, likelihood_invariant = true
    )
    function ac_logdet(vc, x, x′, ctx)
        r = Comrade.value(ctx.view, x, (:sky, :a)) / Comrade.value(ctx.view, x′, (:sky, :a))
        c′ = vc * r
        return log(r) + log((vc - lo) * (hi - vc)) - log((c′ - lo) * (hi - c′))
    end
    θc = [merge(θ, (; sky = merge(θ.sky, (; c = 1.0 + 0.5 * k)))) for (k, θ) in enumerate(θs)]
    @test check_move(ac(ac_logdet), post, θc; rng, τ = 0.01).logdet < 1.0e-6
    @test_throws "reports logdet = 0.0, but the finite-difference Jacobian" check_move(ac((vc, x, x′, ctx) -> 0.0), post, θc; rng, τ = 0.01)

    @testset "discrete move: $(Comrade._space_name(space))" for space in (nothing, sp)
        r = check_move(_SheetToy(), post, θs; space, rng, nprior = 400)
        @test r.logdet < 1.0e-6
        @test only(r.ks) <= r.ks_critical
    end

    # each property fails with its own message
    a_only = CompensatedMove("a_only", view, (:sky, :a), (:sky, :b), (vb, x, x′, ctx) -> vb; likelihood_invariant = true)
    @test Comrade.likelihood_invariant(a_only)
    @test_throws "move a_only changed the log-likelihood" check_move(a_only, post, θs; rng)
    rwa = CompensatedMove("rw", view, (:sky, :a), (:sky, :b), (vb, x, x′, ctx) -> vb)
    @test !Comrade.likelihood_invariant(rwa)
    @test check_move(rwa, post, θs; rng, nprior = 400).loglikelihood == 0
    skew = CompensatedMove(
        "skew", view, (:sky, :a), (:sky, :b),
        (vb, x, x′, ctx) -> vb * exp(Comrade.value(ctx.view, x, (:sky, :a)) - 1)
    )
    @test_throws "move skew is not reversed by its reverse step" check_move(skew, post, θs; rng)
    # a random walk that does not leave the prior invariant: it moves `a` but ignores the
    # log-determinant of the b rescaling it makes
    biased = CompensatedMove(
        "biased", view, (:sky, :a), (:sky, :b), (vb, x, x′, ctx) -> vb .* 1.0;
        logdet = (vb, x, x′, ctx) -> 2 * (x′[1] - x[1])
    )
    @test_throws "move biased reports logdet" check_move(biased, post, θs; rng)
    @test_throws "does not leave the prior invariant" check_move(
        biased, post, θs; rng, τ = 0.3, logdet_atol = Inf, nprior = 400, nsteps = 50
    )

    @test_throws "lies in the block" CompensatedMove("x", view, (:sky, :a), (:sky,), identity)
    @test_throws "index 2 is outside" CompensatedMove("x", view, (:sky, :a), (:sky, :b), identity; index = 2)
    @test_throws "the initial_scale of move x must be positive, got 0" MoveSet(
        post, (CompensatedMove("x", view, (:sky, :a), (:sky, :b), identity),); initial_scale = 0
    )
    @test_throws "the step distribution of move badstep is a Symbol" MoveSet(post, (_BadStep(),))
    @test_throws "the step distribution of move badstep is a Symbol" check_move(_BadStep(), post, θs)
    @test_throws "at least one component" ComponentwiseNormalStep(0)
end

# A small model with a gain log-amplitude chain `lg` (fitted OU hyperparameters) and a
# real-line phase chain `gp` with a reference site: V = f1 · gᵢ conj(gⱼ) · Gaussian.
_moves_gm_data() = load_test_table()
function _moves_gm_sky()
    tm(θ, _, _) = θ.f1 * stretched(Gaussian(), θ.σ1, θ.σ1)
    prior = (f1 = VLBIUniform(0.5, 1.5), σ1 = VLBIUniform(μas2rad(10.0), μas2rad(40.0)))
    return SkyModel(tm, prior, spatialgrid(μas2rad(150.0), μas2rad(150.0), 16, 16))
end
@instrument function _moves_gm_instrument()
    return @jones begin
        lg ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = Exponential(0.1), τ = Exponential(2.0))))
        gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 0.5, τ = 2.0)); refant = SensitiveReference(0.0))
        return Jones(exp(complex(lg, gp)))
    end
end
@instrument function _moves_wrapped_instrument()
    return @jones begin
        gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), WrappedOrnsteinUhlenbeck(σ = 0.5, τ = 2.0)); refant = SensitiveReference(0.0))
        return Jones(exp(complex(0.0, gp)))
    end
end

# `move` without its `likelihood_invariant` flag.
struct _NotInvariant{M} <: Comrade.AbstractMove
    move::M
end
Comrade.step_distribution(w::_NotInvariant) = Comrade.step_distribution(w.move)
Comrade.move_name(::_NotInvariant) = "not_invariant"
Comrade.propose(w::_NotInvariant, x, u, ctx) = Comrade.propose(w.move, x, u, ctx)
Comrade.component_coords(w::_NotInvariant) = Comrade.component_coords(w.move)

struct _TracedDiscrete <: Comrade.AbstractMove end
Comrade.step_distribution(::_TracedDiscrete) = SymmetricDiscreteStep(rng -> 1, -)
Comrade.move_name(::_TracedDiscrete) = "traced_discrete"
Comrade.traceable(::_TracedDiscrete) = true

mutable struct _MovesState
    position::Vector{Float64}
end

@testset "chain moves and MoveSet" begin
    post = VLBIPosterior(_moves_gm_data(); sky = _moves_gm_sky(), instrument = _moves_gm_instrument())
    sp = Comrade.PT.StdNormal()
    rng = Random.Xoshiro(11)
    θs = [prior_sample(rng, post) for _ in 1:2]

    @testset "chain moves: $(Comrade._space_name(space))" for space in (nothing, sp)
        view = CoordinateView(post, space)
        sheet = PhaseSheetMove(post, (:gp,); space)
        @test Comrade.move_name(sheet) == "phase_sheet"
        # the reference site's points are fixed; every other site has free points
        @test !isempty(sheet.points) && any(sheet.fixed[1]) && !all(sheet.fixed[1])
        r = check_move(sheet, post, θs; space, rng, nprior = 300)
        @test r.logdet < 1.0e-6 && only(unique(length(r.ks))) >= 1
        x = Comrade.inverse(view.tbase, θs[1])
        ctx = Comrade.move_context(view, (sheet,))
        step = (1, first(sheet.points[1][4]), 1)
        x′, _ = Comrade.propose(sheet, x, step, ctx)
        d = parent(Comrade.value(view, x′, (:instrument, :gp))) .- parent(θs[1].instrument.gp)
        @test all(v -> isapprox(v, 0; atol = 1.0e-9) || isapprox(v, 2π; atol = 1.0e-9), d)
        @test count(v -> isapprox(v, 2π; atol = 1.0e-9), d) >= 1
        @test Comrade.step_label(sheet, step) == (:gp, sheet.points[1][2])

        fg = flux_gain_move(view; flux = (:sky, :f1), gains = (:instrument, :lg))
        r = check_move(fg, post, θs; space, rng, nprior = 300)
        @test r.logdet < 1.0e-6
        x′, _ = Comrade.propose(fg, x, 0.2, ctx)
        f, f′ = θs[1].sky.f1, Comrade.value(view, x′, (:sky, :f1))
        @test parent(Comrade.value(view, x′, (:instrument, :lg)).params) ≈
            parent(θs[1].instrument.lg.params) .- log(f′ / f) / 2

        hm = chain_hyper_moves(view, (:instrument, :lg))
        @test sort(Comrade.move_name.(hm)) == ["chain_hyper[lg.σ]", "chain_hyper[lg.τ]"]
        for m in hm
            @test length(m.sites) > 1 && all(!isempty, m.innovations)
            @test Comrade.step_distribution(m) == ComponentwiseNormalStep(length(m.sites))
            r = check_move(m, post, θs; space, rng, nprior = m === first(hm) ? 300 : 0)
            @test r.logdet < 1.0e-6
            x′, ld = Comrade.propose(m, x, fill(0.3, length(m.sites)), ctx)
            v′ = Comrade.value(view, x′, (:instrument, :lg))
            @test parent(v′.params) ≈ parent(θs[1].instrument.lg.params) rtol = 1.0e-10
            @test all(s -> v′.hyperparams[s] != θs[1].instrument.lg.hyperparams[s], m.sites)
            @test length(ld) == length(m.sites)
        end
    end
    @test_throws "is not a Gauss–Markov chain term with fitted hyperparameters" chain_hyper_moves(CoordinateView(post), (:instrument, :gp))
    @test MoveSet(post, chain_hyper_moves(CoordinateView(post), (:instrument, :lg))) isa MoveSet
    @test_throws "is componentwise, which needs a move that leaves the likelihood unchanged" MoveSet(post, (_NotInvariant(first(chain_hyper_moves(CoordinateView(post), (:instrument, :lg)))),))
    @test_throws "may be wrapped" PhaseSheetMove(VLBIPosterior(_moves_gm_data(); sky = _moves_gm_sky(), instrument = _moves_wrapped_instrument()), (:gp,))
    @test_throws "no instrument parameter gq" PhaseSheetMove(post, (:gq,))
    @test_throws "power must be positive" flux_gain_move(CoordinateView(post); flux = (:sky, :f1), gains = (:instrument, :lg), power = 0)

    @testset "MoveSet" begin
        view = CoordinateView(post)
        fg = flux_gain_move(view; flux = (:sky, :f1), gains = (:instrument, :lg))
        sheet = PhaseSheetMove(post, (:gp,))
        dir = mktempdir()
        ms = MoveSet(
            post, (fg, sheet); θ0 = θs[1], rounds = [5, 3], initial_scale = [0.05, nothing],
            output = joinpath(dir, "moves.jls")
        )
        @test move_summary(ms)[1].τ ≈ 0.05
        tflat = asflat(post)
        x0 = Comrade.inverse(tflat, θs[1])
        @test move_seconds(ms) == 0
        x1 = ms(copy(x0), tflat, :warmup, Random.Xoshiro(1); step = 10, total = 100)
        @test move_seconds(ms) > 0
        @test x1 != x0
        s = move_summary(ms)
        @test [x.warmup.proposed for x in s] == [5, 3]
        @test s[1].τ != 0.05 && isnothing(s[2].τ)
        @test sum(v[1] for v in values(s[2].bystep)) == 3
        @test Comrade.deserialize(joinpath(dir, "moves.jls"))[1].warmup == s[1].warmup
        # the likelihood is invariant under both moves
        @test Comrade.loglikelihood(post, transform(tflat, x1)) ≈ Comrade.loglikelihood(post, θs[1]) rtol = 1.0e-8
        # scales are frozen during sampling
        τ = s[1].τ
        x2 = ms(x1, tflat, :sampling, Random.Xoshiro(2))
        @test move_summary(ms)[1].τ == τ && move_summary(ms)[1].sampling.proposed == 5
        x32 = ms(Float32.(x1), tflat, :sampling, Random.Xoshiro(2))
        @test eltype(x32) == Float32 && x32 != Float32.(x1)
        # the hook form moves `state.position`; a warmup call logs only when `info.report` is
        # true or absent, with the acceptances of every call since the last table
        warm(k) = (; phase = :warmup, step = 10k, total = 100)
        st = _MovesState(copy(x2))
        @test_logs ms(st, tflat, merge(warm(2), (; report = false)), Random.Xoshiro(3))
        @test st.position != x2
        @test_logs (:info, r"moves after warmup step 30/100 \(2 calls from step 20, .* s\)(.|\n)*flux_gain\s*│\s*\d+\s*│\s*10\s*│") ms(
            st, tflat, merge(warm(3), (; report = true)), Random.Xoshiro(4)
        )
        @test_logs (:info, r"moves after warmup step 40/100 \(1 call, ") ms(st, tflat, warm(4), Random.Xoshiro(5))
        @test_logs (:info, r"moves after warmup \(1 call, ") ms(st.position, tflat, :warmup, Random.Xoshiro(6))
        # a call that accepts nothing returns its input
        never = CompensatedMove(
            "never", view, (:sky, :f1), (:instrument, :lg), (v, x, x′, ctx) -> v;
            logdet = (vb, x, x′, ctx) -> -Inf
        )
        @test MoveSet(post, (never,))(x0, tflat, :sampling, Random.Xoshiro(7)) === x0

        # StdNormal: invariant moves are accepted without the likelihood
        vs = CoordinateView(post, sp)
        hm = chain_hyper_moves(vs, (:instrument, :lg))
        fgs = flux_gain_move(vs; flux = (:sky, :f1), gains = (:instrument, :lg))
        mss = MoveSet(post, (fgs, hm...); space = sp, θ0 = θs[1], rounds = [4, 3, 3])
        @test all(mss.free)
        tstd = Comrade.transport_to(post, sp)
        ss = _MovesState(Comrade.inverse(tstd, θs[1]))
        for k in 1:5
            mss(ss, tstd, warm(k), Random.Xoshiro(20 + k))
        end
        sm = move_summary(mss)
        nsite = [length(m.sites) for m in hm]
        @test [x.warmup.proposed for x in sm] == 5 .* [4, 3 * nsite[1], 3 * nsite[2]]
        @test all(x -> x.warmup.accepted > 0, sm)
        @test sm[2].τ isa Vector && length(sm[2].τ) == nsite[1] && sm[2].τ != fill(0.1, nsite[1])
        @test sum(v[1] for v in values(sm[2].bystep)) == 15 * nsite[1]
        θe = transform(tstd, ss.position)
        @test Comrade.loglikelihood(post, θe) ≈ Comrade.loglikelihood(post, θs[1]) rtol = 1.0e-8
        @test θe.instrument.lg.hyperparams != θs[1].instrument.lg.hyperparams
        # The acceptance never reads the likelihood, so prior draws passed through the
        # MoveSet stay prior draws on the stepped hyperparameter coordinates.
        msp = MoveSet(post, hm; space = sp, rounds = 10, initial_scale = 1.0)
        idx = reduce(vcat, [m.hcoords for m in hm])
        prng = Random.Xoshiro(31)
        nd = 300
        moved = map(1:nd) do _
            x = randn(prng, Comrade.dimension(tstd))
            for _ in 1:3
                x = msp(x, tstd, :sampling, prng)
            end
            x[idx]
        end
        @test all(x -> x.sampling.accepted > 0, move_summary(msp))
        fresh = [randn(prng, length(idx)) for _ in 1:nd]
        crit = Comrade._ks_critical(1.0e-3, length(idx), nd)
        @test all(i -> Comrade._ks_statistic(getindex.(moved, i), getindex.(fresh, i)) <= crit, eachindex(idx))
        # a move flagged invariant that changes the likelihood fails the audit
        bad = CompensatedMove("bad", vs, (:sky, :f1), (:instrument, :lg), (vb, x, x′, ctx) -> vb; likelihood_invariant = true)
        msb = MoveSet(post, (bad,); space = sp, rounds = 20)
        @test_throws "changed the likelihood" msb(Comrade.inverse(tstd, θs[1]), tstd, :warmup, Random.Xoshiro(1))

        @test_throws "the moves act in the StdNormal space but the sampler samples the flat space" MoveSet(post, (PhaseSheetMove(post, (:gp,); space = sp),); space = sp)(st, tflat, warm(1), rng)
        @test_throws "rounds has 1 entries for 2 moves" MoveSet(post, (fg, sheet); rounds = [1])
        @test_throws "initial_scale has 3 entries for 2 moves" MoveSet(post, (fg, sheet); initial_scale = [0.1, 0.1, 0.1])
        @test_throws "the initial_scale of move flux_gain must be positive, got nothing" MoveSet(post, (fg, sheet); initial_scale = [nothing, 0.1])
        @test_throws "target_accept must lie in (0, 1)" MoveSet(post, (fg,); target_accept = 1.0)
        @test_throws "move names must be distinct" MoveSet(post, (fg, fg))
        @test_throws "MoveSet needs at least one move" MoveSet(post, ())
        @test_throws "unknown sampler phase" ms(st, tflat, (; phase = :burnin, step = 1, total = 1), rng)
        @test_throws "unknown sampler phase" ms(x0, tflat, :burnin, rng)
        @test_throws "is traceable but does not have a NormalStep" MoveSet(post, (_TracedDiscrete(),))

        # with the scale fixed (sampling), a traceable move gives the same chain as
        # proposing one step at a time
        rw(traceable) = CompensatedMove("f1_rw", view, (:sky, :f1), (:instrument, :lg), (v, x, x′, ctx) -> v; traceable)
        msf, msq = MoveSet(post, (rw(true),); rounds = 7), MoveSet(post, (rw(false),); rounds = 7)
        xf0 = Comrade.inverse(tflat, θs[2])
        xf = msf(copy(xf0), tflat, :sampling, Random.Xoshiro(4))
        xq = msq(copy(xf0), tflat, :sampling, Random.Xoshiro(4))
        @test xf == xq && xf != xf0
        @test only(move_summary(msf)).sampling == only(move_summary(msq)).sampling

        lgonly = CompensatedMove("lg_only", view, (:sky, :f1), (:instrument, :lg), (v, x, x′, ctx) -> v; likelihood_invariant = true)
        @test_throws "move lg_only changed the log-likelihood" MoveSet(post, (lgonly,); θ0 = θs[1])
    end
end
