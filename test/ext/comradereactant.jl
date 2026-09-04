using Reactant
using Random
using Serialization
using Distributions
import TransformVariables as TV

const ReactantEx = Comrade.ComradeBase.ReactantEx

# Reference-antenna gains fix some sites to a constant value. Rebuilding the full
# parameter vector used to scatter those constants into a freshly-allocated array
# (`yfv[fixed_index] .= fixed_values`), which forces scalar indexing and fails to
# trace under Reactant. This checks the gather-based path traces and matches the CPU
# result for both the flat and cube transforms.
@testset "PartiallyFixedTransform under Reactant" begin
    dist = product_distribution([Normal(0.0, 1.0), Normal(0.0, 1.0), Normal(0.0, 1.0)])
    variate_index = [1, 2, 4]
    fixed_index = [3, 5]
    fixed_values = [7.0, 9.0]
    pcd = Comrade.PartiallyConditionedDist(dist, variate_index, fixed_index, fixed_values)

    # flat path (the gradient path NUTS uses): the raw TransformVariables node, on which
    # `TV.transform_with` runs directly. `asflat` now returns a `TransportedDistribution`
    # wrapper (no `transform_with` method), so grab the node via `transport_node`.
    let t = Comrade.transport_node(pcd, Comrade.TVFlat())
        x = rand(TV.dimension(t))
        y, _, _ = TV.transform_with(TV.LogJac(), t, x, firstindex(x))
        @test y[fixed_index] == fixed_values
        f(xx) = sum(first(TV.transform_with(TV.LogJac(), t, xx, 1)))
        @test convert(Float64, @jit f(Reactant.to_rarray(x))) ≈ sum(y)
    end

    # cube path: the forward transport must still place the fixed values correctly
    # on the CPU. (The inner cube transform itself is not yet Reactant-traceable,
    # independent of this fix, so only the flat path is jit'd.)
    let t = ascube(pcd)
        u = rand(TV.dimension(t))
        y = Comrade.latent_pfwd(t, u)
        @test y[fixed_index] == fixed_values
    end
end

@testset "ComradeReactantExt" begin

    # NB: closures (lcamp/cphase) carry SparseMatrixCSC design matrices that
    # Reactant.to_rarray cannot currently convert; use complex visibilities here
    # so the data tuple round-trips cleanly onto the device (matches the NeuralFields
    # example pattern).
    _, vis, _, _, _ = load_data()
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 32, 32)
    skym = SkyModel(test_model, test_prior(), g)
    post_cpu = VLBIPosterior(skym, vis; admode = nothing)
    post = Comrade.prepare_device(post_cpu, ReactantEx())

    tpost = asflat(post)

    # logdensity round-trips through Reactant
    x0_r = Reactant.to_rarray(prior_sample(Random.default_rng(), tpost))
    ld = @jit logdensityof(tpost, x0_r)
    @test isfinite(convert(Float64, ld))

    # Small warmup so the test stays in CI budget. The Stan windowed schedule is
    # run internally by the ProbProg NUTS pass; only n_adapts is configurable here.
    s = ReactantNUTS(; n_adapts = 50, max_tree_depth = 6)

    # MemoryStore
    res = sample(post, s, 100; chunk_size = 25)
    chain = res.out
    @test length(Comrade.postsamples(chain)) == 100
    @test hasproperty(samplerstats(chain), :numerical_error)
    # cross-backend contract: per-draw Bool divergence flags (`count`-able by the
    # default disk callback)
    @test eltype(samplerstats(chain).numerical_error) == Bool
    @test haskey(samplerinfo(chain), :warmup_history)
    @test haskey(samplerinfo(chain), :sample_history)

    # `sample` also returns the live final MCMC state.
    @test hasproperty(res.state, :position)
    @test length(Array(res.state.position)) == dimension(tpost)

    show(IOBuffer(), MIME"text/plain"(), chain)

    # DiskStore + restart
    dir = mktempdir()
    out = sample(post, s, 100; saveto = DiskStore(name = dir, stride = 25)).out
    @test out isa Comrade.DiskOutput
    @test out.nsamples == 100
    @test isfile(joinpath(dir, "state.jls"))
    @test isfile(joinpath(dir, "metadata.jls"))

    # Metadata.jls captured the sample + warmup history.
    let meta = open(deserialize, joinpath(dir, "metadata.jls"))
        @test meta[:sampler] == :ReactantNUTS
        @test haskey(meta, :warmup_history)
        @test haskey(meta, :sample_history)
    end

    # restart skips (already-completed) warmup, reloads state.jls, and appends
    # sampling chunks to reach the new total.
    out2 = sample(post, s, 200; saveto = DiskStore(name = dir, stride = 25), restart = true).out
    @test out2.nsamples == 200
    c = load_samples(out2)
    @test length(Comrade.postsamples(c)) == 200

    # a fresh run into a directory that already holds a chain warns, clears it, and
    # starts over; `restart = true` is the way to continue
    out3 = (@test_logs (:warn, r"already contains a sampled chain") match_mode = :any sample(
        post, s, 100; saveto = DiskStore(name = dir, stride = 25)
    )).out
    @test out3.nsamples == 100
    @test length(Comrade.postsamples(load_samples(out3))) == 100

    rm(dir, recursive = true)
end

@testset "ReactantNUTS warmup checkpoint/resume" begin
    ext = Base.get_extension(Comrade, :ComradeReactantExt)
    ProbProg = Reactant.ProbProg

    # Gaussian (StdNormal) priors so the prior logdensity traces under Reactant — the bounded
    # StdUniform path is not relevant here. Warmup correctness is independent of the prior.
    gprior = (
        f1 = VLBIGaussian(1.0, 0.1), σ1 = VLBIGaussian(μas2rad(20.0), μas2rad(2.0)),
        τ1 = VLBIGaussian(0.5, 0.05), ξ1 = VLBIGaussian(0.0, 0.3),
        f2 = VLBIGaussian(0.5, 0.1), σ2 = VLBIGaussian(μas2rad(20.0), μas2rad(2.0)),
        τ2 = VLBIGaussian(0.5, 0.05), ξ2 = VLBIGaussian(0.0, 0.3),
        x = VLBIGaussian(0.0, μas2rad(20.0)), y = VLBIGaussian(0.0, μas2rad(20.0)),
    )

    _, vis, _, _, _ = load_data()
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 12, 12)
    skym = SkyModel(test_model, gprior, g)
    post = Comrade.prepare_device(VLBIPosterior(skym, vis; admode = nothing), ReactantEx())
    tpost = asflat(post)
    ldf = ext._default_ldf

    na = 30
    sampler = ReactantNUTS(; n_adapts = na, max_tree_depth = 4, init_step_size = 0.01)
    x0 = Reactant.to_rarray(prior_sample(Random.default_rng(), tpost))
    freshrng() = Reactant.ReactantRNG(Reactant.to_rarray(UInt64[1, 5]))
    quiet = _ -> nothing

    # Chunked warmup is bit-identical to a single fused warmup (windows are anchored to the
    # global warmup length via total_warmup/warmup_offset).
    s_single, _ = ext.warmup_chunked(freshrng(), ldf, x0, tpost, sampler; chunk = na, callback = quiet)
    s_multi, _ = ext.warmup_chunked(freshrng(), ldf, x0, tpost, sampler; chunk = 10, callback = quiet)
    @test Array(s_single.step_size) == Array(s_multi.step_size)
    @test Array(s_single.inverse_mass_matrix) == Array(s_multi.inverse_mass_matrix)
    @test Array(s_single.position) == Array(s_multi.position)

    # Mid-warmup checkpoint (captured through the warmup callback) round-trips through
    # save_state/load_state carrying the adaptation accumulators, and resuming from it
    # reproduces the full-warmup state exactly.
    ckpt = tempname()
    capture = info -> (info.step == 20 && ProbProg.save_state(ckpt, info.state); nothing)
    ext.warmup_chunked(freshrng(), ldf, x0, tpost, sampler; chunk = 10, callback = capture)
    @test isfile(ckpt)
    loaded = ProbProg.load_state(ckpt)
    rm(ckpt; force = true)
    @test loaded.adaptation !== nothing

    s_resume, _ = ext.warmup_chunked(
        freshrng(), ldf, nothing, tpost, sampler;
        chunk = 10, resume_state = loaded, warmup_done = 20, callback = quiet,
    )
    @test Array(s_resume.step_size) == Array(s_single.step_size)
    @test Array(s_resume.inverse_mass_matrix) == Array(s_single.inverse_mass_matrix)
    @test Array(s_resume.position) == Array(s_single.position)
end

@testset "ReactantNUTS FisherLowRank warmup refits" begin
    ext = Base.get_extension(Comrade, :ComradeReactantExt)

    gprior = (
        f1 = VLBIGaussian(1.0, 0.1), σ1 = VLBIGaussian(μas2rad(20.0), μas2rad(2.0)),
        τ1 = VLBIGaussian(0.5, 0.05), ξ1 = VLBIGaussian(0.0, 0.3),
        f2 = VLBIGaussian(0.5, 0.1), σ2 = VLBIGaussian(μas2rad(20.0), μas2rad(2.0)),
        τ2 = VLBIGaussian(0.5, 0.05), ξ2 = VLBIGaussian(0.0, 0.3),
        x = VLBIGaussian(0.0, μas2rad(20.0)), y = VLBIGaussian(0.0, μas2rad(20.0)),
    )
    _, vis, _, _, _ = load_data()
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 12, 12)
    post = Comrade.prepare_device(
        VLBIPosterior(SkyModel(test_model, gprior, g), vis; admode = nothing), ReactantEx()
    )
    tpost = asflat(post)
    ldf = ext._default_ldf
    x0 = Reactant.to_rarray(prior_sample(Random.default_rng(), tpost))
    freshrng() = Reactant.ReactantRNG(Reactant.to_rarray(UInt64[1, 5]))
    quiet = _ -> nothing

    na, chunk = 40, 2
    adaptor = FisherLowRank(; rank = 4, schedule = [0.5], min_draws = 6)
    sampler = ReactantNUTS(;
        n_adapts = na, max_tree_depth = 4, init_step_size = 0.01,
        metric_adaptor = adaptor
    )
    @test Comrade.metric_refit_steps(adaptor, na, chunk) == [20]

    # The adaptation window is only checkpointed alongside the sampler state: a window
    # persisted without the state it belongs to would let a restart refit from draws that
    # do not match the position it resumes at.
    sckpt, ckpt = tempname(), tempname()
    state, _, newt = ext.warmup_chunked(
        freshrng(), ldf, x0, tpost, sampler;
        chunk, callback = quiet, checkpoint = sckpt, adaptation_checkpoint = ckpt,
    )

    # Warmup ran to completion in a latent space the refit installed: the transform is a
    # different object and now carries a preconditioner, while base-flat dimension is
    # unchanged (the refit reparameterizes, it does not project).
    @test newt !== tpost
    pre = Comrade._transport_pre(newt)
    @test pre isa LowRankPreconditioner
    @test length(pre.b) == dimension(tpost)
    # The installed transform's low-rank block is padded to a fixed device rank cap so
    # later refits can overwrite it in place; only the non-unit columns are corrections.
    scales = Array(pre.s)
    @test 0 < count(!=(1.0), scales) <= adaptor.rank
    @test length(scales) >= count(!=(1.0), scales)
    @test all(isfinite, Array(state.position))

    # Draws and scores were accumulated from the sampler itself, one per chunk, and
    # checkpointed for restart.
    @test isfile(ckpt)
    st = deserialize(ckpt)
    rm(ckpt; force = true)
    rm(sckpt; force = true)
    @test st isa Comrade.FisherAdaptation
    @test length(st.draws) == na ÷ chunk
    @test length(st.scores) == length(st.draws)
    @test all(d -> length(d) == dimension(tpost), st.draws)
    @test all(g -> all(isfinite, g), st.scores)

    # Resuming hands the accumulated window back, so a post-restart refit sees the whole
    # run rather than only the steps since the restart.
    _, _, resumed = ext.warmup_chunked(
        freshrng(), ldf, nothing, newt, sampler;
        chunk, callback = quiet, resume_state = state,
        warmup_done = na - 4, adaptation = st, segment_start = 20,
    )
    @test length(st.draws) > na ÷ chunk
    @test Comrade._transport_pre(resumed) isa LowRankPreconditioner

    # WelfordDiagonal keeps the sampler's own adaptation and never touches the transform.
    wsampler = ReactantNUTS(; n_adapts = na, max_tree_depth = 4, init_step_size = 0.01)
    _, _, samet = ext.warmup_chunked(
        freshrng(), ldf, x0, tpost, wsampler; chunk = na, callback = quiet
    )
    @test samet === tpost
end

# GaussMarkov (time-correlated) instrument priors trace through the flat path with no
# Reactant-specific code: the `@trace`d chain logpdf and the branchless whitened coloring
# only need `rgetindex`/`rsetindex!` for the scalar-indexing opt-in (Reactant promotes
# the captured static chain tables itself). Pointwise-evaluated distributions
# (hyperpriors, IID overrides) must accept `::Number`, so use the VLBI*/PT distributions.
# Each process stresses a different traced path: OU the fitted-hyperparameter coloring,
# BrownianMotion+FixedInit the `scatter_values!` fixed-value fill, WrappedBrownian the
# log-sum-exp wrapped-normal image sum plus, in its default non-centered form, the masked
# wrap of the batched affine scan and, in its centered form, the sheet-weight loop (whose
# anchor read is a table-driven index into the latent segment); `AngleEmbedded` the
# angle-embedding transform, and WrappedOrnsteinUhlenbeck the shortest-arc `_cond_mean`
# and the pointwise (carry-free) unwrap of the centered lift.
@testset "GaussMarkov priors under Reactant" begin
    using Enzyme

    _, vis, _, _, _ = load_data()
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 12, 12)
    skym = SkyModel(test_model, test_prior(), g)

    @instrument function gmint_reactant()
        return @jones begin
            lg ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = VLBIExponential(0.1), τ = VLBIExponential(2.0))))
            gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 2.0, τ = VLBIInverseGamma(3.0, 6.0))); refant = SEFDReference(0.0))
            return SingleStokesGain(exp(complex(lg, gp)))
        end
    end

    @instrument function gmint_reactant_bm()
        return @jones begin
            lg ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), BrownianMotion(D = VLBIExponential(0.1)); init = FixedInit(0.0)))
            return SingleStokesGain(exp(lg))
        end
    end

    @instrument function gmint_reactant_wb()
        return @jones begin
            gp ~ ArrayPrior(
                GaussMarkovSitePrior(ScanSeg(), WrappedBrownian(τ = VLBIInverseGamma(3.0, 6.0)); init = UniformInit());
                refant = SEFDReference(0.0)
            )
            return SingleStokesGain(exp(1im * gp))
        end
    end

    @instrument function gmint_reactant_wbc()
        return @jones begin
            gp ~ ArrayPrior(
                GaussMarkovSitePrior(
                    ScanSeg(), WrappedBrownian(τ = VLBIInverseGamma(3.0, 6.0));
                    init = UniformInit(), param = Centered()
                );
                refant = SEFDReference(0.0)
            )
            return SingleStokesGain(exp(1im * gp))
        end
    end

    @instrument function gmint_reactant_wbe()
        return @jones begin
            gp ~ ArrayPrior(
                GaussMarkovSitePrior(
                    ScanSeg(), WrappedBrownian(τ = VLBIInverseGamma(3.0, 6.0));
                    init = UniformInit(), param = AngleEmbedded()
                );
                refant = SEFDReference(0.0)
            )
            return SingleStokesGain(exp(1im * gp))
        end
    end

    @instrument function gmint_reactant_wou()
        return @jones begin
            gp ~ ArrayPrior(
                GaussMarkovSitePrior(
                    ScanSeg(),
                    WrappedOrnsteinUhlenbeck(σ = VLBIExponential(0.3), τ = VLBIInverseGamma(3.0, 6.0))
                );
                refant = SEFDReference(0.0)
            )
            return SingleStokesGain(exp(1im * gp))
        end
    end

    # gradient (what the compiled NUTS kernels differentiate) matches CPU Enzyme
    fgrad(x, tp) = Enzyme.gradient(
        Enzyme.set_runtime_activity(Enzyme.Reverse), Const(Base.Fix1(logdensityof, tp)), x
    )[1]

    nwhile(hlo) = length(collect(eachmatch(r"stablehlo\.while", repr(hlo))))

    function check_matches_cpu(intm)
        post_cpu = VLBIPosterior(skym, intm, vis; admode = nothing)
        tpost_cpu = asflat(post_cpu)
        x0 = prior_sample(Random.Xoshiro(31), tpost_cpu)

        post = Comrade.prepare_device(post_cpu, ReactantEx())
        tpost = asflat(post)
        x0_r = Reactant.to_rarray(x0)

        # regression guard: the whole flat path — coloring (staged writes + affine
        # scan + scatter), chain logpdf (re-read-style gathers), fixed-value fill —
        # must stay fully raised. A serialized `stablehlo.while` costs ~100-200 us of
        # device latency per leapfrog step, doubled in the gradient.
        @test nwhile(Reactant.@code_hlo optimize = true logdensityof(tpost, x0_r)) == 0
        @test nwhile(Reactant.@code_hlo optimize = true fgrad(x0_r, tpost)) == 0

        # value matches the CPU path (whitened hierarchical transform + chain logpdf +
        # refant conditioning all traced)
        ld = @jit logdensityof(tpost, x0_r)
        @test convert(Float64, ld) ≈ logdensityof(tpost_cpu, x0) rtol = 1.0e-10

        g_r = convert(Vector{Float64}, @jit fgrad(x0_r, tpost))
        g_c = fgrad(x0, tpost_cpu)
        @test g_r ≈ g_c rtol = 1.0e-8
        return nothing
    end

    @testset "OrnsteinUhlenbeck" begin
        check_matches_cpu(gmint_reactant())
    end
    @testset "BrownianMotion FixedInit" begin
        check_matches_cpu(gmint_reactant_bm())
    end
    @testset "WrappedBrownian UniformInit refant (non-centered default)" begin
        check_matches_cpu(gmint_reactant_wb())
    end
    @testset "WrappedBrownian centered" begin
        check_matches_cpu(gmint_reactant_wbc())
    end
    @testset "WrappedBrownian angle embedding" begin
        check_matches_cpu(gmint_reactant_wbe())
    end
    @testset "WrappedOrnsteinUhlenbeck StationaryInit refant (centered default)" begin
        check_matches_cpu(gmint_reactant_wou())
    end
end

# Device staging for the low-rank preconditioner. The buffers are padded to a fixed rank
# so that a warmup refit can overwrite them in place: the compiled program sees the same
# shapes and is not rebuilt. Padding columns carry `s = 1` and are exact no-ops.
@testset "LowRankPreconditioner device buffers" begin
    rng = Random.Xoshiro(7)
    n, m = 40, 3
    V = Matrix(qr(randn(rng, n, m)).Q)[:, 1:m]
    p = LowRankPreconditioner(randn(rng, n), exp.(randn(rng, n)), V, [4.0, 2.0, 0.5])

    dev = Comrade._device_pre(p; rank_cap = 8)
    @test size(dev.V) == (n, 8)
    @test Array(dev.V)[:, 1:m] ≈ V
    @test all(iszero, Array(dev.V)[:, (m + 1):end])
    @test Array(dev.s) == [4.0, 2.0, 0.5, 1, 1, 1, 1, 1]
    @test_throws ArgumentError Comrade._device_pre(p; rank_cap = 2)

    # A host vector goes through the host mirror of the device buffers, so the padded
    # transform still agrees with the host fit it was built from.
    z = randn(rng, n)
    @test Comrade._affine_fwd(Comrade._pre_for(dev, z), z) ≈ Comrade._affine_fwd(p, z)

    # In-place update at a lower rank: the tail returns to the no-op padding.
    p2 = LowRankPreconditioner(randn(rng, n), exp.(randn(rng, n)), V[:, 1:2], [3.0, 0.4])
    Comrade._update_device_pre!(dev, p2)
    @test Array(dev.s) == [3.0, 0.4, 1, 1, 1, 1, 1, 1]
    @test Array(dev.b) ≈ p2.b
    @test Comrade._affine_fwd(Comrade._pre_for(dev, z), z) ≈ Comrade._affine_fwd(p2, z)

    p3 = LowRankPreconditioner(randn(rng, n), exp.(randn(rng, n)), Matrix(qr(randn(rng, n, 9)).Q)[:, 1:9], fill(2.0, 9))
    @test_throws ArgumentError Comrade._update_device_pre!(dev, p3)
end
