using Random
import Distributions as Dists

# Geweke joint-distribution tests (`check_move` with `ngeweke`, and `check_moveset`) of every
# move kind.

# `move` with the sign of its log-determinant flipped.
struct _FlippedLogdet{M} <: Comrade.AbstractMove
    move::M
end
Comrade.step_distribution(w::_FlippedLogdet) = Comrade.step_distribution(w.move)
Comrade.move_name(w::_FlippedLogdet) = "flipped_" * Comrade.move_name(w.move)
Comrade.likelihood_invariant(w::_FlippedLogdet) = Comrade.likelihood_invariant(w.move)
Comrade.component_coords(w::_FlippedLogdet) = Comrade.component_coords(w.move)
function Comrade.propose(w::_FlippedLogdet, x, step, ctx)
    x′, ld = Comrade.propose(w.move, x, step, ctx)
    return x′, -ld
end

# The random-walk `move` with its step shifted by `δ`, so its proposal is not symmetric.
struct _Drift{M} <: Comrade.AbstractMove
    move::M
    δ::Float64
end
Comrade.step_distribution(w::_Drift) = Comrade.step_distribution(w.move)
Comrade.move_name(w::_Drift) = "drift_" * Comrade.move_name(w.move)
Comrade.likelihood_invariant(w::_Drift) = Comrade.likelihood_invariant(w.move)
Comrade.propose(w::_Drift, x, u, ctx) = Comrade.propose(w.move, x, u + w.δ, ctx)

# A wide phase chain, so that its 2π sheet moves are accepted often enough to test.
@instrument function _geweke_instrument()
    return @jones begin
        lg ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = Exponential(0.1), τ = Exponential(2.0))))
        gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 4.0, τ = 2.0)); refant = SEFDReference(0.0))
        return SingleStokesGain(exp(complex(lg, gp)))
    end
end

# A phase chain with fitted hyperparameters and SEFD referencing. The τ prior spans
# correlation times both shorter and longer than the scan gaps.
@instrument function _geweke_sefd_instrument()
    return @jones begin
        gp ~ ArrayPrior(
            GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = Exponential(1.0), τ = Dists.LogNormal(log(0.02), 1.0)));
            refant = SEFDReference(0.0)
        )
        return SingleStokesGain(exp(complex(0.0, gp)))
    end
end

_noisier(d, k) = Comrade.EHTObservationTable{Comrade.datumtype(d)}(
    Comrade.measurement(d), k .* Comrade.noise(d), Comrade.arrayconfig(d)
)

# Every statistic is under the critical value and every move was accepted at least sometimes.
_passes(r) = all(<=(r.geweke_critical), r.geweke) && all(>(0.05), r.geweke_accept)

@testset "joint-distribution test of moves" begin
    sp = Comrade.PT.StdNormal()
    kw = (; ngeweke = 1000, geweke_steps = 20)

    @testset "toy sky: $(Comrade._space_name(space))" for space in (nothing, sp)
        _, vis, _, _, _ = load_data()
        prior = (
            a = Dists.LogNormal(0.0, 0.3), b = Dists.LogNormal(0.0, 0.3),
            c = Dists.Uniform(0.2, 5.0), φ = Dists.Normal(0.0, 3.0),
        )
        # total S/N of about 12: the likelihood shapes the posterior but still allows
        # large steps for the non-invariant moves
        g = imagepixels(μas2rad(150.0), μas2rad(150.0), 8, 8)
        post = VLBIPosterior(SkyModel(_moves_toy_sky, prior, g), _noisier(vis, 1000.0))
        rng = Random.Xoshiro(41)
        θs = [prior_sample(rng, post) for _ in 1:2]
        view = CoordinateView(post, space)
        ab = CompensatedMove(
            "a_b", view, (:sky, :a), (:sky, :b),
            (vb, x, x′, ctx) -> vb * Comrade.value(ctx.view, x, (:sky, :a)) / Comrade.value(ctx.view, x′, (:sky, :a));
            likelihood_invariant = true
        )
        rw(traceable) = CompensatedMove("a_rw", view, (:sky, :a), (:sky, :b), (vb, x, x′, ctx) -> vb; traceable)
        @test _passes(check_move(ab, post, θs; space, rng, τ = 0.3, kw...))
        @test _passes(check_move(rw(false), post, θs; space, rng, kw...))
        @test _passes(check_move(rw(true), post, θs; space, rng, kw...))
        @test _passes(check_move(_SheetToy(), post, θs; space, rng, kw...))
        @test_throws "is not reversed by its reverse step" check_move(_Drift(rw(false), 0.05), post, θs; space, rng)
        @test_throws "does not leave the posterior invariant" check_move(_Drift(rw(false), 0.05), post, θs; space, rng, rtol = Inf, kw...)
    end

    post = VLBIPosterior(_moves_gm_sky(), _geweke_instrument(), _moves_gm_data())
    rng = Random.Xoshiro(43)
    θs = [prior_sample(rng, post) for _ in 1:2]
    @testset "chain moves: $(Comrade._space_name(space))" for space in (nothing, sp)
        view = CoordinateView(post, space)
        fg = flux_gain_move(view; flux = (:sky, :f1), gains = (:instrument, :lg))
        @test _passes(check_move(fg, post, θs; space, rng, kw...))
        @test _passes(check_move(PhaseSheetMove(post, (:gp,); space), post, θs; space, rng, kw...))
    end

    @testset "componentwise and MoveSet: $(Comrade._space_name(space))" for space in (nothing, sp)
        view = CoordinateView(post, space)
        fg = flux_gain_move(view; flux = (:sky, :f1), gains = (:instrument, :lg))
        hm = chain_hyper_moves(view, (:instrument, :lg))
        @testset "componentwise $(Comrade.move_name(m))" for m in hm
            @test _passes(check_move(m, post, θs; space, rng, τ = 0.5, kw...))
        end
        sheet = PhaseSheetMove(post, (:gp,); space)
        @test _passes(
            check_moveset(
                MoveSet(
                    post, (fg, sheet, hm...); space, rounds = [2, 1, 1, 1],
                    initial_scale = [0.05, nothing, 0.5, 0.5]
                ); rng, kw...
            )
        )
    end
    vs = CoordinateView(post, sp)
    fgs = flux_gain_move(vs; flux = (:sky, :f1), gains = (:instrument, :lg))
    hm = chain_hyper_moves(vs, (:instrument, :lg))
    @test_throws "ngeweke must be positive" check_moveset(MoveSet(post, (fgs,); space = sp); ngeweke = 0)

    rpost = VLBIPosterior(_moves_gm_sky(), _geweke_sefd_instrument(), _moves_gm_data())
    rview = CoordinateView(rpost, sp)
    rθs = [prior_sample(rng, rpost) for _ in 1:2]
    @testset "SEFD-referenced chain: $(Comrade.move_name(m))" for m in chain_hyper_moves(rview, (:instrument, :gp))
        ctx = Comrade.move_context(rview, (m,))
        for _ in 1:20
            x = Comrade.inverse(rview.tbase, prior_sample(rng, rpost))
            @test isnothing(Comrade._check_components(m, x, 0.5 .* randn(rng, length(m.sites)), ctx))
        end
        @test _passes(check_move(m, rpost, rθs; space = sp, rng, τ = 0.5, kw...))
    end

    @testset "broken moves fail" begin
        @test_throws "does not leave the posterior invariant" check_move(
            _FlippedLogdet(first(hm)), post, θs; space = sp, rng, τ = 0.5, logdet_atol = Inf, kw...
        )
        @test_throws "does not leave the posterior invariant" check_move(
            _Drift(fgs, 0.025), post, θs; space = sp, rng, τ = 0.05, rtol = Inf, kw...
        )
        liar = CompensatedMove("liar", vs, (:sky, :f1), (:instrument, :lg), (vb, x, x′, ctx) -> vb; likelihood_invariant = true)
        @test_throws "changed the likelihood" check_move(liar, post, θs; space = sp, rng, rtol = Inf, kw...)
    end
end
