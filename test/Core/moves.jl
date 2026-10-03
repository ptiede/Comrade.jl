using LinearAlgebra
using Random
import Distributions as Dists
using Distributions: Exponential
using VLBIFiles

# Toy skies: the flux depends on `a * b` (and `a * c`) only, and the shift on `φ` only through
# `cos φ`, `sin φ`.
function _moves_toy_sky(θ, meta)
    g = stretched(Gaussian(), μas2rad(20.0), μas2rad(20.0))
    return (θ.a * θ.b * θ.c) * shifted(g, μas2rad(5.0) * cos(θ.φ), μas2rad(5.0) * sin(θ.φ))
end

struct _SheetToy <: Comrade.AbstractMove end
Comrade.step_kind(::_SheetToy) = DiscreteSymmetric()
Comrade.move_name(::_SheetToy) = "sheet"
Comrade.draw_step(::_SheetToy, rng::AbstractRNG, τ) = rand(rng, (-1, 1))
Comrade.reverse_step(::_SheetToy, s) = -s
function Comrade.propose(::_SheetToy, x, s, ctx)
    path = (:sky, :φ)
    x′ = copy(x)
    x′[Comrade.coords(ctx.view, path)] = Comrade.latent(ctx.view, path, Comrade.value(ctx.view, x, path) + 2π * s)
    return x′, 0.0
end

struct _NoStep <: Comrade.AbstractMove end
Comrade.step_kind(::_NoStep) = DiscreteSymmetric()
Comrade.move_name(::_NoStep) = "nostep"

@testset "moves" begin
    _, vis, _, _, _ = load_data()
    prior = (
        a = Dists.LogNormal(0.0, 0.3), b = Dists.LogNormal(0.0, 0.3),
        c = Dists.Uniform(0.2, 5.0), φ = Dists.Normal(0.0, 3.0),
    )
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 8, 8)
    post = VLBIPosterior(SkyModel(_moves_toy_sky, prior, g), vis)
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
        initial_scale = 0.3
    )
    @testset "compensated move: $(Comrade._space_name(space))" for space in (nothing, sp)
        m = ab(CoordinateView(post, space))
        @test Comrade.step_kind(m) == RandomWalk(0.3)
        r = check_move(m, post, θs; space, rng, nprior = 400)
        @test r.reversal < 1.0e-12 && r.logdet < 1.0e-6
        @test length(r.ks) == 2 && all(<=(r.ks_critical), r.ks)
    end

    # c ↦ c a / a′ for a uniform c: in the flat (logit) space the Jacobian is not one
    view = CoordinateView(post)
    lo, hi = 0.2, 5.0
    ac(logdet) = CompensatedMove(
        "a_c", view, (:sky, :a), (:sky, :c),
        (vc, x, x′, ctx) -> vc * Comrade.value(ctx.view, x, (:sky, :a)) / Comrade.value(ctx.view, x′, (:sky, :a));
        logdet, initial_scale = 0.01
    )
    function ac_logdet(vc, x, x′, ctx)
        r = Comrade.value(ctx.view, x, (:sky, :a)) / Comrade.value(ctx.view, x′, (:sky, :a))
        c′ = vc * r
        return log(r) + log((vc - lo) * (hi - vc)) - log((c′ - lo) * (hi - c′))
    end
    θc = [merge(θ, (; sky = merge(θ.sky, (; c = 1.0 + 0.5 * k)))) for (k, θ) in enumerate(θs)]
    @test check_move(ac(ac_logdet), post, θc; rng).logdet < 1.0e-6
    @test_throws "reports logdet = 0.0, but the finite-difference Jacobian" check_move(ac((vc, x, x′, ctx) -> 0.0), post, θc; rng)

    @testset "discrete move: $(Comrade._space_name(space))" for space in (nothing, sp)
        r = check_move(_SheetToy(), post, θs; space, rng, nprior = 400)
        @test r.logdet < 1.0e-6
        @test only(r.ks) <= r.ks_critical
    end

    # each property fails with its own message
    a_only = CompensatedMove("a_only", view, (:sky, :a), (:sky, :b), (vb, x, x′, ctx) -> vb)
    @test_throws "move a_only changed the log-likelihood" check_move(a_only, post, θs; rng)
    @test !Comrade.is_invariant(CompensatedMove("rw", view, (:sky, :a), (:sky, :b), (vb, x, x′, ctx) -> vb; invariant = false))
    @test check_move(
        CompensatedMove("rw", view, (:sky, :a), (:sky, :b), (vb, x, x′, ctx) -> vb; invariant = false),
        post, θs; rng, nprior = 400
    ).loglikelihood == 0
    skew = CompensatedMove(
        "skew", view, (:sky, :a), (:sky, :b),
        (vb, x, x′, ctx) -> vb * exp(Comrade.value(ctx.view, x, (:sky, :a)) - 1);
        invariant = false
    )
    @test_throws "move skew is not reversed by reverse_step" check_move(skew, post, θs; rng)
    # a random walk that does not leave the prior invariant: it moves `a` but ignores the
    # log-determinant of the b rescaling it makes
    biased = CompensatedMove(
        "biased", view, (:sky, :a), (:sky, :b), (vb, x, x′, ctx) -> vb .* 1.0;
        logdet = (vb, x, x′, ctx) -> 2 * (x′[1] - x[1]), invariant = false
    )
    @test_throws "move biased reports logdet" check_move(biased, post, θs; rng)
    @test_throws "does not leave the prior invariant" check_move(
        biased, post, θs; rng, τ = 0.3, logdet_atol = Inf, nprior = 400, nsteps = 50
    )

    @test_throws "lies in the block" CompensatedMove("x", view, (:sky, :a), (:sky,), identity)
    @test_throws "index 2 is outside" CompensatedMove("x", view, (:sky, :a), (:sky, :b), identity; index = 2)
    @test_throws "initial_scale must be positive" CompensatedMove("x", view, (:sky, :a), (:sky, :b), identity; initial_scale = 0)
    @test_throws MethodError Comrade.draw_step(_NoStep(), rng, 1.0)
    @test_throws MethodError Comrade.reverse_step(_NoStep(), 1)
end

# A small model with a gain log-amplitude chain `lg` (fitted OU hyperparameters) and a
# real-line phase chain `gp` with a reference site: V = f1 · gᵢ conj(gⱼ) · Gaussian.
_moves_gm_data() = extract_table(
    VLBIFiles.load(VLBIFiles.UVData, joinpath(@__DIR__, "..", "test_data.uvfits")),
    Visibilities(; time_average = VLBI.GapBasedScans())
)
function _moves_gm_sky()
    tm(θ, meta) = θ.f1 * stretched(Gaussian(), θ.σ1, θ.σ1)
    prior = (f1 = VLBIUniform(0.5, 1.5), σ1 = VLBIUniform(μas2rad(10.0), μas2rad(40.0)))
    return SkyModel(tm, prior, imagepixels(μas2rad(150.0), μas2rad(150.0), 16, 16))
end
@instrument function _moves_gm_instrument()
    return @jones begin
        lg ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = Exponential(0.1), τ = Exponential(2.0))))
        gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 0.5, τ = 2.0)); refant = SEFDReference(0.0))
        return SingleStokesGain(exp(complex(lg, gp)))
    end
end

mutable struct _MovesState
    position::Vector{Float64}
end

@testset "chain moves and MoveSet" begin
    post = VLBIPosterior(_moves_gm_sky(), _moves_gm_instrument(), _moves_gm_data())
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
    end
    @test_throws "no instrument parameter gq" PhaseSheetMove(post, (:gq,))
    @test_throws "power must be positive" flux_gain_move(CoordinateView(post); flux = (:sky, :f1), gains = (:instrument, :lg), power = 0)

    @testset "MoveSet" begin
        view = CoordinateView(post)
        fg = flux_gain_move(view; flux = (:sky, :f1), gains = (:instrument, :lg))
        sheet = PhaseSheetMove(post, (:gp,))
        dir = mktempdir()
        ms = MoveSet(post, (fg, sheet); θ0 = θs[1], rounds = [5, 3], output = joinpath(dir, "moves.jls"))
        tflat = asflat(post)
        st = _MovesState(Comrade.inverse(tflat, θs[1]))
        warm(k) = (; phase = :warmup, step = 10k, total = 100)
        ms(st, tflat, warm(1), Random.Xoshiro(1))
        s = move_summary(ms)
        @test [x.warmup.proposed for x in s] == [5, 3]
        @test s[1].τ != 0.05 && isnothing(s[2].τ)
        @test sum(v[1] for v in values(s[2].bystep)) == 3
        @test Comrade.deserialize(joinpath(dir, "moves.jls"))[1].warmup == s[1].warmup
        # the likelihood is invariant under both moves
        @test loglikelihood(post, transform(tflat, st.position)) ≈ loglikelihood(post, θs[1]) rtol = 1.0e-8
        # scales are frozen during sampling
        τ = s[1].τ
        ms(st, tflat, (; phase = :sampling, step = 1, total = 10), Random.Xoshiro(2))
        @test move_summary(ms)[1].τ == τ && move_summary(ms)[1].sampling.proposed == 5

        # through a preconditioner the base-space chain and the decisions are the same
        n = dimension(tflat)
        V = Matrix(qr(randn(Random.Xoshiro(5), n, 3)).Q)[:, 1:3]
        pre = LowRankPreconditioner(0.1 .* randn(Random.Xoshiro(6), n), exp.(0.2 .* randn(Random.Xoshiro(7), n)), V, [0.5, 2.0, 0.8])
        tpre = Comrade.transport_to(post, pre)
        msa = MoveSet(post, (fg, sheet); rounds = [8, 8])
        msb = MoveSet(post, (fg, sheet); rounds = [8, 8])
        x0 = Comrade.inverse(tflat, θs[2])
        sa = _MovesState(copy(x0))
        sb = _MovesState(Comrade._affine_inv(pre, x0))
        msa(sa, tflat, warm(1), Random.Xoshiro(9))
        msb(sb, tpre, warm(1), Random.Xoshiro(9))
        @test Comrade._affine_fwd(pre, sb.position) ≈ sa.position rtol = 1.0e-8
        @test [x.warmup.accepted for x in move_summary(msa)] == [x.warmup.accepted for x in move_summary(msb)]
        @test sa.position != x0

        @test_throws "the moves act in the StdNormal space but the sampler samples the flat space" MoveSet(post, (PhaseSheetMove(post, (:gp,); space = sp),); space = sp)(st, tflat, warm(1), rng)
        @test_throws "rounds has 1 entries for 2 moves" MoveSet(post, (fg, sheet); rounds = [1])
        @test_throws "target_accept must lie in (0, 1)" MoveSet(post, (fg,); target_accept = 1.0)
        @test_throws "move names must be distinct" MoveSet(post, (fg, fg))
        @test_throws "MoveSet needs at least one move" MoveSet(post, ())
        @test_throws "unknown sampler phase" ms(st, tflat, (; phase = :burnin, step = 1, total = 1), rng)
        lgonly = CompensatedMove("lg_only", view, (:sky, :f1), (:instrument, :lg), (v, x, x′, ctx) -> v)
        @test_throws "move lg_only changed the log-likelihood" MoveSet(post, (lgonly,); θ0 = θs[1])
    end
end
