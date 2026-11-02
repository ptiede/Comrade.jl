using Reactant
using Random
using Serialization
using Distributions
import TransformVariables as TV

const ReactantEx = Comrade.ComradeBase.ReactantEx

mutable struct _DeviceMovesState
    position::Any
end

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
    # the remaining per-draw stats: the stored log density is that of the flat posterior
    # at the draw, the step size is frozen after warmup, and the wall time is positive
    let st = samplerstats(chain), tcpu = asflat(post_cpu), ps = Comrade.postsamples(chain)
        @test length(st.log_density) == length(st.step_size) == length(st.time) == 100
        for i in (1, 25, 26, 100)
            @test st.log_density[i] ≈ logdensityof(tcpu, Comrade.inverse(tcpu, ps[i])) rtol = 1.0e-8
        end
        @test allequal(st.step_size) && first(st.step_size) > 0
        @test all(>(0), st.time)
    end
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
    # the scan files carry the same per-draw stats as the in-memory chain
    @test propertynames(samplerstats(load_samples(out))) ==
        (:numerical_error, :log_density, :step_size, :time)

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
    out3 = (
        @test_logs (:warn, r"already contains a sampled chain") match_mode = :any sample(
            post, s, 100; saveto = DiskStore(name = dir, stride = 25)
        )
    ).out
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

    # Checkpoints land every `checkpoint_every` steps, independently of the chunk length: a
    # warmup interrupted at step 24 resumes from the checkpoint at step 20 and reproduces the
    # full-warmup state.
    sckpt, pckpt = tempname(), tempname()
    interrupt = info -> (info.step == 24 && error("interrupted"); nothing)
    @test_throws "interrupted" ext.warmup_chunked(
        freshrng(), ldf, x0, tpost, sampler;
        chunk = 4, callback = interrupt, checkpoint = sckpt, progress_checkpoint = pckpt,
        checkpoint_every = 10,
    )
    @test deserialize(pckpt).warmup_done == 20
    s_ckpt, _ = ext.warmup_chunked(
        freshrng(), ldf, nothing, tpost, sampler;
        chunk = 4, callback = quiet, resume_state = ProbProg.load_state(sckpt), warmup_done = 20,
    )
    rm(sckpt; force = true)
    rm(pckpt; force = true)
    @test Array(s_ckpt.step_size) == Array(s_single.step_size)
    @test Array(s_ckpt.inverse_mass_matrix) == Array(s_single.inverse_mass_matrix)
    @test Array(s_ckpt.position) == Array(s_single.position)
end

@testset "ReactantNUTS in the StdNormal space" begin
    gprior = (
        f1 = VLBIGaussian(1.0, 0.1), σ1 = VLBIGaussian(μas2rad(20.0), μas2rad(2.0)),
        τ1 = VLBIGaussian(0.5, 0.05), ξ1 = VLBIGaussian(0.0, 0.3),
        f2 = VLBIGaussian(0.5, 0.1), σ2 = VLBIGaussian(μas2rad(20.0), μas2rad(2.0)),
        τ2 = VLBIGaussian(0.5, 0.05), ξ2 = VLBIGaussian(0.0, 0.3),
        x = VLBIGaussian(0.0, μas2rad(20.0)), y = VLBIGaussian(0.0, μas2rad(20.0)),
    )
    _, vis, _, _, _ = load_data()
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 12, 12)
    post_cpu = VLBIPosterior(SkyModel(test_model, gprior, g), vis; admode = nothing)
    post = Comrade.prepare_device(post_cpu, ReactantEx())
    sp = Comrade.PT.StdNormal()

    # MoveSet between NUTS chunks on a device position: the log density is compiled once and
    # every round is proposed
    rwmove = CompensatedMove(
        "f1_rw", CoordinateView(post_cpu, sp), (:sky, :f1), (:sky, :f2), (v, x, x′, ctx) -> v
    )
    ms = MoveSet(post_cpu, (rwmove,); space = sp, rounds = 4)
    mdir = mktempdir()
    mout = sample(
        post, ReactantNUTS(; n_adapts = 20, max_tree_depth = 4, init_step_size = 0.01), 20;
        saveto = DiskStore(name = mdir, stride = 10), transport_method = sp, between_chunks = ms
    ).out
    msum = only(move_summary(ms))
    @test msum.warmup.proposed == 4 * 2 && msum.sampling.proposed == 4 * 1
    @test 0 < msum.warmup.accepted + msum.sampling.accepted
    @test all(p -> all(isfinite, values(p.sky)), Comrade.postsamples(load_samples(mout)))
    rm(mdir; recursive = true)

    # On a device position a traceable move proposes in a compiled program and a host move
    # from a copy of the point; both make the same proposals and decisions as the host
    # kernels on the same random numbers.
    mview = CoordinateView(post_cpu, sp)
    scaled(v, x, x′, ctx) = v * sum(ctx.c)
    devmove = CompensatedMove(
        "f1_dev", mview, (:sky, :f1), (:sky, :f2), scaled;
        context = (; c = [1.0]), traceable = true
    )
    devmove2 = CompensatedMove(
        "τ1_dev", mview, (:sky, :τ1), (:sky, :f2), scaled;
        context = (; c = [1.0]), traceable = true
    )
    hostmove = CompensatedMove("ξ1_host", mview, (:sky, :ξ1), (:sky, :f2), scaled; context = (; c = [1.0]))
    @test Comrade.traceable(devmove) && !Comrade.traceable(hostmove)
    dmoves, drounds = (devmove, devmove2, hostmove), [6, 3, 4]
    thost = Comrade.transport_to(post_cpu, sp)
    tdev = Comrade.transport_to(post, sp)
    xstart = Comrade.inverse(thost, prior_sample(Random.Xoshiro(3), post_cpu))
    msd = MoveSet(post_cpu, dmoves; space = sp, rounds = drounds)
    msh = MoveSet(post_cpu, dmoves; space = sp, rounds = drounds)
    warm = (; phase = :warmup, step = 2, total = 20)
    sh = msh(_DeviceMovesState(copy(xstart)), thost, warm, Random.Xoshiro(30))
    sd = msd(_DeviceMovesState(Reactant.to_rarray(xstart)), tdev, warm, Random.Xoshiro(30))
    @test sd.position isa Reactant.AbstractConcreteArray
    @test Array(sd.position) ≈ sh.position rtol = 1.0e-9
    @test [m.warmup for m in move_summary(msd)] == [m.warmup for m in move_summary(msh)]
    @test [m.warmup.proposed for m in move_summary(msd)] == drounds
    @test [m.τ for m in move_summary(msd)] ≈ [m.τ for m in move_summary(msh)]
    @test sh.position != xstart

    # A context update reaches the compiled program without recompiling it: with the scale
    # `c` doubled, the device makes the host's proposals and decisions. The host context holds
    # the moves' own arrays, so a MoveSet built from the same moves sees the update too.
    k0 = msd.compiled[]
    Comrade.update_context!(msd, :c, [2.0])
    msh2 = MoveSet(post_cpu, dmoves; space = sp, rounds = drounds)
    @test only(msh2.ctx.c) == 2.0
    msd.stats .= [Comrade._MoveStats(d, τ) for (d, τ) in zip(msd.dists, msd.initial_scale)]
    sh2 = msh2(_DeviceMovesState(copy(xstart)), thost, warm, Random.Xoshiro(33))
    sd2 = msd(_DeviceMovesState(Reactant.to_rarray(xstart)), tdev, warm, Random.Xoshiro(33))
    @test Array(sd2.position) ≈ sh2.position rtol = 1.0e-9
    @test msd.compiled[] === k0
    @test_throws DimensionMismatch Comrade.update_context!(msd, :c, [1.0, 2.0])
    @test_throws "has no entry nope" Comrade.update_context!(msd, :nope, [1.0])

    # Invariant moves are accepted without the likelihood, in the compiled program and on
    # the host alike: `u` and `w` do not enter the sky, so random walks of them are exact.
    post_u = VLBIPosterior(SkyModel(test_model, merge(gprior, (u = VLBIGaussian(0.0, 1.0), w = VLBIGaussian(0.0, 1.0))), g), vis; admode = nothing)
    uview = CoordinateView(post_u, sp)
    keep(v, x, x′, ctx) = v
    umove = CompensatedMove("u_dev", uview, (:sky, :u), (:sky, :f2), keep; likelihood_invariant = true, traceable = true)
    wmove = CompensatedMove("w_host", uview, (:sky, :w), (:sky, :f2), keep; likelihood_invariant = true)
    tdev_u = Comrade.transport_to(Comrade.prepare_device(post_u, ReactantEx()), sp)
    thost_u = Comrade.transport_to(post_u, sp)
    xu = Comrade.inverse(thost_u, prior_sample(Random.Xoshiro(4), post_u))
    msu_h = MoveSet(post_u, (umove, wmove); space = sp, rounds = [5, 3])
    msu_d = MoveSet(post_u, (umove, wmove); space = sp, rounds = [5, 3])
    @test all(msu_d.free)
    su_h = msu_h(_DeviceMovesState(copy(xu)), thost_u, warm, Random.Xoshiro(31))
    su_d = msu_d(_DeviceMovesState(Reactant.to_rarray(xu)), tdev_u, warm, Random.Xoshiro(31))
    @test Array(su_d.position) ≈ su_h.position rtol = 1.0e-9
    @test [m.warmup for m in move_summary(msu_d)] == [m.warmup for m in move_summary(msu_h)]
    @test [m.τ for m in move_summary(msu_d)] ≈ [m.τ for m in move_summary(msu_h)]
    @test su_h.position[end] != xu[end]

    # a full run to disk in the StdNormal space: draws are finite and the stored transport
    # records the space
    dir = mktempdir()
    out = sample(
        post, ReactantNUTS(; n_adapts = 40, max_tree_depth = 4, init_step_size = 0.01), 20;
        saveto = DiskStore(name = dir, stride = 10), transport_method = sp
    ).out
    @test out.nsamples == 20
    @test all(p -> all(isfinite, values(p.sky)), Comrade.postsamples(load_samples(out)))
    @test deserialize(joinpath(dir, "transport.jls")) isa Comrade.PT.StdNormal
    rm(dir; recursive = true)
end

# Random-walk Metropolis–Hastings along the direction `v` for the log-density `target`,
# keeping its own counters.
mutable struct LineMove{F}
    target::F
    v::Vector{Float64}
    scale::Float64
    nprop::Int
    nacc::Int
    phases::Vector{Symbol}
end
LineMove(target, v, scale) = LineMove(target, v, scale, 0, 0, Symbol[])
function (m::LineMove)(state, tpost, info, rng)
    push!(m.phases, info.phase)
    z = vec(Array(state.position))
    for _ in 1:5
        z1 = z .+ m.scale * randn(rng) .* m.v
        m.nprop += 1
        if log(rand(rng)) < m.target(z1) - m.target(z)
            z = z1
            m.nacc += 1
        end
    end
    state.position = z
    return state
end

# Host-side moves between NUTS chunks (`between_chunks`). The toy target ignores the
# posterior and is a diagonal Gaussian in the latent coordinates, so the sampled moments
# are known exactly and the host can evaluate the target for a Metropolis–Hastings step.
@testset "ReactantNUTS between_chunks moves" begin
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
    n = dimension(asflat(post))
    μ = collect(range(-1.0, 1.0; length = n))
    σ = collect(range(0.5, 2.0; length = n))
    target(x) = -sum(abs2, (x .- μ) ./ σ) / 2
    ldf(x, _) = target(x)
    x0 = Comrade.transform(asflat(post), zeros(n))
    freshrng() = Reactant.ReactantRNG(Reactant.to_rarray(UInt64[3, 11]))
    s = ReactantNUTS(; n_adapts = 60, max_tree_depth = 5)
    runchain(hook; nsamples = 60, host_rng = Random.Xoshiro(2)) = sample(
        freshrng(), post, s, nsamples;
        chunk_size = 10, initial_params = x0, ldf, host_rng, between_chunks = hook
    )
    draws(res) = reduce(hcat, (Comrade.inverse(asflat(post), p) for p in Comrade.postsamples(res.out)))

    @testset "no hook and an identity hook give the same chain" begin
        r0 = runchain(nothing)
        r1 = runchain((st, _, _, _) -> st)
        @test draws(r0) == draws(r1)
        @test Array(r0.state.position) == Array(r1.state.position)
        @test Array(r0.state.step_size) == Array(r1.state.step_size)
    end

    @testset "the hook runs at every chunk end and reports there without checkpointing" begin
        calls = Tuple{Symbol, Int}[]
        reports = Int[]
        hook = function (st, _, info, _)
            push!(calls, (info.phase, info.step))
            get(info, :report, false) && push!(reports, info.step)
            return st
        end
        r = runchain(hook)
        wsteps = [step for (phase, step) in calls if phase === :warmup]
        @test !isempty(wsteps) && last(wsteps) == 60 && allequal(diff([0; wsteps]))
        @test reports == wsteps
        @test [step for (phase, step) in calls if phase === :sampling] == collect(10:10:50)
        @test length(Comrade.postsamples(r.out)) == 60
    end

    @testset "moved states restart NUTS from a recomputed gradient" begin
        # A host-array position drops the cached gradient and potential energy and comes
        # back on the device in the kernel's shape; an unmoved state keeps both.
        st = runchain(nothing; nsamples = 20).state
        info = (; phase = :sampling, step = 0, total = 0)
        kept = ext._run_between_chunks((s, _, _, _) -> s, st, nothing, info, Random.Xoshiro(1))
        @test kept.gradient !== nothing && kept.potential_energy !== nothing
        znew = vec(Array(st.position)) .+ 0.1
        moved = ext._run_between_chunks(
            (s, _, _, _) -> (s.position = znew; s), st, nothing, info, Random.Xoshiro(1)
        )
        @test moved.gradient === nothing && moved.potential_energy === nothing
        @test moved.position isa Reactant.ConcreteRArray
        @test size(moved.position) == (1, n)
        @test vec(Array(moved.position)) == znew

        # One chunk from the moved state: the kernel compiled for an empty gradient slot
        # runs, and returns a potential energy evaluated at its own position.
        st2, _, _ = ext.sample_chunked(moved, ldf, asflat(post), s; num_samples = 4, chunk_size = 4)
        @test -only(Array(st2.potential_energy)) ≈ target(vec(Array(st2.position))) rtol = 1.0e-10

        @test_throws "must return the MCMCState" ext._run_between_chunks(
            (_, _, _, _) -> nothing, st2, nothing, info, Random.Xoshiro(1)
        )
        @test_throws DimensionMismatch ext._run_between_chunks(
            (s, _, _, _) -> (s.position = zeros(n + 1); s), st2, nothing, info, Random.Xoshiro(1)
        )
    end

    @testset "a valid MH move leaves the target moments unchanged" begin
        v = normalize!([1.0; 1.0; zeros(n - 2)])
        mv = LineMove(target, v, 1.5)
        res = runchain(mv; nsamples = 3000)
        @test :warmup in mv.phases && :sampling in mv.phases
        @test count(==(:warmup), mv.phases) == s.n_adapts ÷ 10
        @test count(==(:sampling), mv.phases) == 3000 ÷ 10 - 1
        @test 0.2 < mv.nacc / mv.nprop < 0.95
        X = draws(res)
        # the stored log densities are those of the stored draws
        @test samplerstats(res.out).log_density ≈ map(target, eachcol(X)) rtol = 1.0e-10
        # NUTS on a Gaussian is close to independent draws: 6 standard errors on the
        # mean, and the sample variance within 20%
        @test all(abs.(vec(mean(X; dims = 2)) .- μ) .< 6 .* σ ./ sqrt(size(X, 2)))
        @test all(abs.(vec(var(X; dims = 2)) ./ σ .^ 2 .- 1) .< 0.2)
    end
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
