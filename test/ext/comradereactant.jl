using Reactant
using Random
using Serialization
using Distributions
using StaticArrays
import TransformVariables as TV

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
    pcd = VLBIInstrumentModels.PartiallyConditionedDist(dist, variate_index, fixed_index, fixed_values)

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

    _, vis, _, _, _ = load_data()
    g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 32, 32)
    skym = SkyModel(test_model, test_prior(), g)
    post_cpu = VLBIPosterior(vis; sky = skym, admode = nothing)
    post = todevice(post_cpu, Reactant.Sharding.NoSharding())

    tpost = asflat(post)

    # logdensity round-trips through Reactant
    x0_r = Reactant.to_rarray(prior_sample(Random.default_rng(), tpost))
    ld = @jit logdensityof(tpost, x0_r)
    @test isfinite(convert(Float64, ld))

    # Small warmup so the test stays in CI budget. The Stan windowed schedule is
    # run internally by the ProbProg NUTS pass; only n_adapts is configurable here.
    s = ReactantNUTS(; n_adapts = 50, max_tree_depth = 6)

    # MemoryStore
    chain = sample(post, s, Int32(100); chunk_size = Int32(25))
    @test chain isa PosteriorSamples
    @test length(Comrade.postsamples(chain)) == 100
    @test hasproperty(samplerstats(chain), :numerical_error)
    # divergence flags are per-draw Bools, as on the other backends
    @test eltype(samplerstats(chain).numerical_error) == Bool
    # stored log density is the flat posterior's at the draw; step size is fixed after
    # warmup; wall time is positive
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

    # the final MCMC state is the one after the last draw
    let st = samplerinfo(chain)[:final_state], tcpu = asflat(post_cpu)
        @test vec(Array(st.position)) ≈ Comrade.inverse(tcpu, last(Comrade.postsamples(chain))) rtol = 1.0e-10
    end

    show(IOBuffer(), MIME"text/plain"(), chain)

    # DiskStore + restart
    dir = mktempdir()
    out = sample(post, s, 100; saveto = DiskStore(name = dir, stride = 25))
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
        @test haskey(meta, :final_state)
    end
    let st = Reactant.ProbProg.load_state(joinpath(dir, "state.jls")), tcpu = asflat(post_cpu)
        @test vec(Array(st.position)) ≈ Comrade.inverse(tcpu, last(Comrade.postsamples(load_samples(out)))) rtol = 1.0e-10
    end

    # restart skips (already-completed) warmup, reloads state.jls, and appends
    # sampling chunks to reach the new total.
    out2 = sample(post, s, 200; saveto = DiskStore(name = dir, stride = 25), restart = true)
    @test out2.nsamples == 200
    c = load_samples(out2)
    @test length(Comrade.postsamples(c)) == 200

    # a restart that asks for no more samples than are on disk draws nothing
    out_done = @test_logs (:warn, r"already on disk") match_mode = :any sample(
        post, s, 200; saveto = DiskStore(name = dir, stride = 25), restart = true
    )
    @test out_done isa Comrade.DiskOutput
    @test out_done.nsamples == 200

    # a fresh run into a directory that already holds a chain warns, clears it, and
    # starts over; `restart = true` is the way to continue
    out3 = @test_logs (:warn, r"already contains a sampled chain") match_mode = :any sample(
        post, s, 100; saveto = DiskStore(name = dir, stride = 25)
    )
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
    g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 12, 12)
    skym = SkyModel(test_model, gprior, g)
    post = todevice(VLBIPosterior(vis; sky = skym, admode = nothing), Reactant.Sharding.NoSharding())
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

    # Checkpoints are written every `checkpoint_every` steps regardless of chunk length:
    # interrupted at step 24, warmup resumes from step 20 and matches the full warmup.
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
    g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 12, 12)
    post_cpu = VLBIPosterior(vis; sky = SkyModel(test_model, gprior, g), admode = nothing)
    post = todevice(post_cpu, Reactant.Sharding.NoSharding())
    sp = Comrade.PT.StdNormal()

    # MoveSet between NUTS chunks on a device position proposes every round
    rwmove = CompensatedMove(
        "f1_rw", CoordinateView(post_cpu, sp), (:sky, :f1), (:sky, :f2), (v, x, x′, ctx) -> v
    )
    ms = MoveSet(post_cpu, (rwmove,); space = sp, rounds = 4)
    mdir = mktempdir()
    mout = sample(
        post, ReactantNUTS(; n_adapts = 20, max_tree_depth = 4, init_step_size = 0.01), 20;
        saveto = DiskStore(name = mdir, stride = 10), transport_method = sp, between_chunks = ms
    )
    msum = only(move_summary(ms))
    @test msum.warmup.proposed == 4 * 2 && msum.sampling.proposed == 4 * 1
    @test 0 < msum.warmup.accepted + msum.sampling.accepted
    @test all(p -> all(isfinite, values(p.sky)), Comrade.postsamples(load_samples(mout)))
    rm(mdir; recursive = true)

    # Traceable (compiled) and host moves on a device position match the host kernels
    # for the same random numbers.
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

    # `update_context!` reaches the compiled program without recompiling. The host context
    # shares the moves' arrays, so a MoveSet built from the same moves sees the update.
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

    # Likelihood-invariant moves skip the likelihood on device and host; `u` and `w` do
    # not enter the sky.
    post_u = VLBIPosterior(vis; sky = SkyModel(test_model, merge(gprior, (u = VLBIGaussian(0.0, 1.0), w = VLBIGaussian(0.0, 1.0))), g), admode = nothing)
    uview = CoordinateView(post_u, sp)
    keep(v, x, x′, ctx) = v
    umove = CompensatedMove("u_dev", uview, (:sky, :u), (:sky, :f2), keep; likelihood_invariant = true, traceable = true)
    wmove = CompensatedMove("w_host", uview, (:sky, :w), (:sky, :f2), keep; likelihood_invariant = true)
    tdev_u = Comrade.transport_to(todevice(post_u, Reactant.Sharding.NoSharding()), sp)
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

    # a full run to disk: draws are finite and the stored transport records the space
    sampler = ReactantNUTS(; n_adapts = 40, max_tree_depth = 4, init_step_size = 0.01)
    dir = mktempdir()
    out = sample(post, sampler, 20; saveto = DiskStore(name = dir, stride = 10), transport_method = sp)
    @test out.nsamples == 20
    ps = Comrade.postsamples(load_samples(out))
    @test all(p -> all(isfinite, values(p.sky)), ps)
    @test deserialize(joinpath(dir, "transport.jls")) isa Comrade.PT.StdNormal
    rm(dir; recursive = true)
end

# Random-walk Metropolis-Hastings along `v` for the log density `target`.
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

# Host-side `between_chunks` moves. The target is a diagonal Gaussian in the latent
# coordinates, so the sampled moments are known.
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
    g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 12, 12)
    post = todevice(
        VLBIPosterior(vis; sky = SkyModel(test_model, gprior, g), admode = nothing), Reactant.Sharding.NoSharding()
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
    draws(res) = reduce(hcat, (Comrade.inverse(asflat(post), p) for p in Comrade.postsamples(res)))

    @testset "no hook and an identity hook give the same chain" begin
        r0 = runchain(nothing)
        r1 = runchain((st, _, _, _) -> st)
        @test draws(r0) == draws(r1)
        @test Array(samplerinfo(r0)[:final_state].position) == Array(samplerinfo(r1)[:final_state].position)
        @test Array(samplerinfo(r0)[:final_state].step_size) == Array(samplerinfo(r1)[:final_state].step_size)
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
        @test length(Comrade.postsamples(r)) == 60
    end

    @testset "moved states restart NUTS from a recomputed gradient" begin
        # A moved state drops the cached gradient and potential energy and returns to the
        # device in the kernel's shape; an unmoved state keeps both.
        st = samplerinfo(runchain(nothing; nsamples = 20))[:final_state]
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

        # one chunk from the moved state returns the potential energy at its own position
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
        @test samplerstats(res).log_density ≈ map(target, eachcol(X)) rtol = 1.0e-10
        # NUTS draws on a Gaussian are nearly independent, so plain standard errors apply
        @test all(abs.(vec(mean(X; dims = 2)) .- μ) .< 6 .* σ ./ sqrt(size(X, 2)))
        @test all(abs.(vec(var(X; dims = 2)) ./ σ .^ 2 .- 1) .< 0.2)
    end
end

# GaussMarkov instrument priors under Reactant match the CPU value and gradient. Hyperpriors
# and IID overrides are evaluated pointwise and must accept `::Number`, hence the VLBI*/PT
# distributions. Each process and form covers a different traced code path.
@testset "GaussMarkov priors under Reactant" begin
    using Enzyme

    _, vis, _, _, _ = load_data()
    g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 12, 12)
    skym = SkyModel(test_model, test_prior(), g)

    @instrument function gmint_reactant()
        return @jones begin
            lg ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = VLBIExponential(0.1), τ = VLBIExponential(2.0))))
            gp ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), OrnsteinUhlenbeck(σ = 2.0, τ = VLBIInverseGamma(3.0, 6.0))); refant = SensitiveReference(0.0))
            return Jones(exp(complex(lg, gp)))
        end
    end

    @instrument function gmint_reactant_bm()
        return @jones begin
            lg ~ ArrayPrior(GaussMarkovSitePrior(ScanSeg(), BrownianMotion(D = VLBIExponential(0.1)); init = FixedInit(0.0)))
            return Jones(exp(lg))
        end
    end

    @instrument function gmint_reactant_wb()
        return @jones begin
            gp ~ ArrayPrior(
                GaussMarkovSitePrior(ScanSeg(), WrappedBrownian(τ = VLBIInverseGamma(3.0, 6.0)); init = UniformInit());
                refant = SensitiveReference(0.0)
            )
            return Jones(exp(1im * gp))
        end
    end

    @instrument function gmint_reactant_wbc()
        return @jones begin
            gp ~ ArrayPrior(
                GaussMarkovSitePrior(
                    ScanSeg(), WrappedBrownian(τ = VLBIInverseGamma(3.0, 6.0));
                    init = UniformInit(), param = Centered()
                );
                refant = SensitiveReference(0.0)
            )
            return Jones(exp(1im * gp))
        end
    end

    @instrument function gmint_reactant_wbe()
        return @jones begin
            gp ~ ArrayPrior(
                GaussMarkovSitePrior(
                    ScanSeg(), WrappedBrownian(τ = VLBIInverseGamma(3.0, 6.0));
                    init = UniformInit(), param = AngleEmbedded()
                );
                refant = SensitiveReference(0.0)
            )
            return Jones(exp(1im * gp))
        end
    end

    @instrument function gmint_reactant_wou()
        return @jones begin
            gp ~ ArrayPrior(
                GaussMarkovSitePrior(
                    ScanSeg(),
                    WrappedOrnsteinUhlenbeck(σ = VLBIExponential(0.3), τ = VLBIInverseGamma(3.0, 6.0))
                );
                refant = SensitiveReference(0.0)
            )
            return Jones(exp(1im * gp))
        end
    end

    fgrad(x, tp) = Enzyme.gradient(
        Enzyme.set_runtime_activity(Enzyme.Reverse), Const(Base.Fix1(logdensityof, tp)), x
    )[1]

    nwhile(hlo) = length(collect(eachmatch(r"stablehlo\.while", repr(hlo))))

    function check_matches_cpu(intm)
        post_cpu = VLBIPosterior(vis; sky = skym, instrument = intm, admode = nothing)
        tpost_cpu = asflat(post_cpu)
        x0 = prior_sample(Random.Xoshiro(31), tpost_cpu)

        post = todevice(post_cpu, Reactant.Sharding.NoSharding())
        tpost = asflat(post)
        x0_r = Reactant.to_rarray(x0)

        # No `stablehlo.while` may remain: each costs ~100-200 us of device latency per
        # leapfrog step, doubled in the gradient.
        @test nwhile(Reactant.@code_hlo optimize = true logdensityof(tpost, x0_r)) == 0
        @test nwhile(Reactant.@code_hlo optimize = true fgrad(x0_r, tpost)) == 0

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

@testset "Instrument terms under Reactant" begin
    using Enzyme
    fgrad(x, tp) = Enzyme.gradient(
        Enzyme.set_runtime_activity(Enzyme.Reverse), Const(Base.Fix1(logdensityof, tp)), x
    )[1]
    check_matches_cpu(skym, intm, data; transport = asflat) =
        check_matches_cpu(VLBIPosterior(data; sky = skym, instrument = intm, admode = nothing); transport)
    function check_matches_cpu(post_cpu; transport = asflat)
        tpost_cpu = transport(post_cpu)
        x0 = prior_sample(Random.Xoshiro(41), tpost_cpu)
        tpost = transport(todevice(post_cpu, Reactant.Sharding.NoSharding()))
        x0_r = Reactant.to_rarray(x0)
        @test convert(Float64, @jit(logdensityof(tpost, x0_r))) ≈ logdensityof(tpost_cpu, x0) rtol = 1.0e-10
        @test convert(Vector{Float64}, @jit(fgrad(x0_r, tpost))) ≈ fgrad(x0, tpost_cpu) rtol = 1.0e-8
        return nothing
    end

    _, vis, _, _, _, dcoh = load_data()
    g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 12, 12)
    g0 = IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))
    Q = Jones((x, c) -> SMatrix{2, 2}(1.0 + 0im, 0, 0, cis(x.δ)); δ = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.3))))

    @testset "plain numbers" begin
        G = Jones((x, c) -> exp(complex(x.lg, x.gp)); lg = ArrayPrior(g0), gp = ArrayPrior(g0; refant = SensitiveReference(0.0)))
        check_matches_cpu(SkyModel(test_model, test_prior(), g), InstrumentModel(G), vis)
    end

    @testset "numbers" begin
        G = Jones((x, c) -> exp(complex(x.lg, x.gp)); lg = ArrayPrior(g0), gp = ArrayPrior(g0; refant = SensitiveReference(0.0)))
        S = Jones((x, c) -> exp(x.s); s = ArrayPrior(g0))
        check_matches_cpu(SkyModel(test_model, test_prior(), g), InstrumentModel(PerStation(G; AA = G * S)), vis)
    end

    @testset "matrices" begin
        G = Jones(
            (x, c) -> gains(exp(complex(x.lgR, x.gpR)), exp(complex(x.lgL, x.gpL)));
            lgR = ArrayPrior(g0), gpR = ArrayPrior(g0; refant = SingleReference(:AA, 0.0)),
            lgL = ArrayPrior(g0), gpL = ArrayPrior(g0; refant = SingleReference(:AA, 0.0)),
        )
        psky(θ, _, _) = PolarizedModel(modify(Gaussian(), Stretch(θ.σ)), 0.1 * modify(Gaussian(), Stretch(θ.σ)), 0.05 * modify(Gaussian(), Stretch(θ.σ)), 0.02 * modify(Gaussian(), Stretch(θ.σ)))
        skym = SkyModel(psky, (; σ = VLBIUniform(μas2rad(10.0), μas2rad(30.0))), g)
        check_matches_cpu(skym, InstrumentModel(PerStation(G; AA = G * Q)), dcoh)
    end

    @testset "self-calibrated station phase" begin
        intm = InstrumentModel(Jones((x, c) -> complex(exp(x.lg)); lg = ArrayPrior(g0)); selfcal = true)
        transport = p -> transport_to(p, Comrade.PT.StdNormal())
        check_matches_cpu(SkyModel(test_model, test_prior(), g), intm, vis; transport)
    end

    @testset "fields" begin
        G = Jones((x, c) -> exp(complex(x.lg, x.gp)); lg = ArrayPrior(g0), gp = ArrayPrior(g0; refant = SensitiveReference(0.0)))
        S = Jones((x, c) -> exp(x.s); s = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.1))))
        skym = SkyModel(test_model, test_prior(), g)
        post = VLBIPosterior((; A = vis, B = later(vis, 1)); sky = (; A = skym, B = skym), instrument = InstrumentModel(G * S), admode = nothing)
        check_matches_cpu(post)
    end
end

@testset "Sharded posteriors" begin
    if Reactant.XLA.REACTANT_XLA_RUNTIME != "IFRT" || length(Reactant.devices()) < 4
        @warn "Sharded posterior tests skipped: they need Reactant's IFRT runtime and XLA_FLAGS=--xla_force_host_platform_device_count=4"
    else
        using Enzyme
        mesh = Reactant.Sharding.Mesh(reshape(Reactant.devices()[1:4], 4), (:d,))
        bytime = ShardLayout(mesh; Ti = :d)
        fgrad(x, tp) = Enzyme.gradient(
            Enzyme.set_runtime_activity(Enzyme.Reverse), Const(Base.Fix1(logdensityof, tp)), x
        )[1]
        function check_matches_cpu(post_cpu, layout; transport = asflat, rtol = 1.0e-10)
            tcpu = transport(post_cpu)
            x0 = prior_sample(Random.Xoshiro(41), tcpu)
            tdev = transport(todevice(post_cpu, layout))
            xr = Reactant.to_rarray(x0; sharding = Reactant.Sharding.Replicated(mesh))
            @test convert(Float64, @jit(logdensityof(tdev, xr))) ≈ logdensityof(tcpu, x0) rtol = rtol
            @test convert(Vector{Float64}, @jit(fgrad(xr, tdev))) ≈ fgrad(x0, tcpu) rtol = 1.0e-8
            return nothing
        end

        # An image cube with `n` time planes spanning the table `t`.
        function timegrid(t, n)
            τ = Comrade.times(t)
            e = range(first(τ) - 1.0e-3, last(τ) + 1.0e-3; length = n + 1)
            return spatialgrid(μas2rad(150.0), μas2rad(150.0), 12, 12) ⊗ frames(Ti, collect(e[1:(end - 1)]), collect(e[2:end]))
        end
        pulse = BSplinePulse{3}()
        planesky(θ, grid, meta) = ContinuousImage(IntensityMap(exp.(meta.base .+ reshape(collect(values(θ)), 1, 1, :)), grid), pulse)
        function planemodel(t, n; width = 25.0)
            g = timegrid(t, n)
            base = [-(x^2 + y^2) / (2 * μas2rad(width)^2) for x in g.X, y in g.Y, _ in 1:n]
            prior = NamedTuple{ntuple(i -> Symbol(:a, i), n)}(ntuple(_ -> VLBIGaussian(0.0, 0.3), n))
            return SkyModel(planesky, prior, g; metadata = (; base))
        end

        _, vis, _, lcamp, cphase, dcoh = load_data()
        g0 = IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))
        G = Jones((x, c) -> exp(complex(x.lg, x.gp)); lg = ArrayPrior(g0), gp = ArrayPrior(g0; refant = SensitiveReference(0.0)))
        post_cpu = VLBIPosterior(vis; sky = planemodel(vis, 4), instrument = InstrumentModel(G), admode = nothing)

        # The device arrays that `x` holds.
        function devicearrays(x, found = Any[], seen = IdDict{Any, Nothing}())
            if x isa Reactant.AbstractConcreteArray
                push!(found, x)
            elseif x isa Array
                isbitstype(eltype(x)) || foreach(y -> devicearrays(y, found, seen), x)
            elseif !(isbits(x) || x isa Union{Module, Type, Symbol, AbstractString} || haskey(seen, x))
                seen[x] = nothing
                foreach(i -> isdefined(x, i) && devicearrays(getfield(x, i), found, seen), 1:fieldcount(typeof(x)))
            end
            return found
        end

        @testset "time planes" begin
            check_matches_cpu(post_cpu, bytime)
            post = todevice(post_cpu, bytime)
            @test !isempty(devicearrays(post))
            @test all(Reactant.Sharding.is_sharded, devicearrays(post))
            obs = only(post.lklhds).obs
            @test length(obs.sharding.device_to_array_slices) == 4
            @test size(obs) == size(post.instrumentmodel.layout.cell1)
            @test size(obs, 1) > size(vis, 1)
        end

        @testset "polarized" begin
            Gp = Jones(
                (x, c) -> gains(exp(complex(x.lgR, x.gpR)), exp(complex(x.lgL, x.gpL)));
                lgR = ArrayPrior(g0), gpR = ArrayPrior(g0; refant = SingleReference(:AA, 0.0)),
                lgL = ArrayPrior(g0), gpL = ArrayPrior(g0; refant = SingleReference(:AA, 0.0)),
            )
            Q = Jones((x, c) -> SMatrix{2, 2}(1.0 + 0im, 0, 0, cis(x.δ)); δ = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.3))))
            psky(θ, _, _) = PolarizedModel(modify(Gaussian(), Stretch(θ.σ)), 0.1 * modify(Gaussian(), Stretch(θ.σ)), 0.05 * modify(Gaussian(), Stretch(θ.σ)), 0.02 * modify(Gaussian(), Stretch(θ.σ)))
            skym = SkyModel(psky, (; σ = VLBIUniform(μas2rad(10.0), μas2rad(30.0))), timegrid(dcoh, 4))
            check_matches_cpu(VLBIPosterior(dcoh; sky = skym, instrument = InstrumentModel(PerStation(Gp; AA = Gp * Q)), admode = nothing), bytime)
        end

        @testset "fields" begin
            post = VLBIPosterior(
                (; A = vis, B = later(vis, 1));
                sky = (; A = planemodel(vis, 4), B = planemodel(later(vis, 1), 8)),
                instrument = InstrumentModel(G),
                admode = nothing
            )
            check_matches_cpu(post, bytime)
        end

        static = SkyModel(test_model, test_prior(), spatialgrid(μas2rad(150.0), μas2rad(150.0), 12, 12))
        @testset "static image split by record" begin
            check_matches_cpu(VLBIPosterior(vis; sky = static, instrument = InstrumentModel(G), admode = nothing), ShardLayout(mesh; Node = :d))
        end

        @testset "fixed component" begin
            fixed = FixedSkyModel(test_model(rand(NamedDist(test_prior())), nothing, nothing), timegrid(vis, 4))
            check_matches_cpu(VLBIPosterior(vis; sky = MultiSkyModel((; a = planemodel(vis, 4), b = fixed)), instrument = InstrumentModel(G), admode = nothing), bytime)
        end

        @testset "cube and static components" begin
            fixed = FixedSkyModel(test_model(rand(Random.Xoshiro(3), NamedDist(test_prior())), nothing, nothing), static.grid)
            check_matches_cpu(VLBIPosterior(vis; sky = MultiSkyModel((; a = planemodel(vis, 4), b = static)), instrument = InstrumentModel(G), admode = nothing), bytime)
            check_matches_cpu(VLBIPosterior(vis; sky = MultiSkyModel((; b = static, a = planemodel(vis, 4))), instrument = InstrumentModel(G), admode = nothing), bytime)
            check_matches_cpu(VLBIPosterior(vis; sky = MultiSkyModel((; a = planemodel(vis, 4), b = fixed)), instrument = InstrumentModel(G), admode = nothing), bytime)
        end

        # Host and device NUFFTs differ by ~1e-8 relative in the model visibilities, which closures
        # and very large χ² carry into the log density.
        @testset "closures" begin
            for c in (cphase, lcamp)
                check_matches_cpu(VLBIPosterior(c; sky = planemodel(vis, 4; width = 8.0), admode = nothing), bytime; rtol = 1.0e-8)
            end
            post = VLBIPosterior(cphase, lcamp; sky = planemodel(vis, 4; width = 8.0), admode = nothing)
            dev = todevice(post, Reactant.Sharding.NoSharding())
            θ = Comrade.transform(asflat(post), prior_sample(Random.Xoshiro(41), asflat(post)))
            _, mv = @jit Comrade.forward_model(dev, Reactant.to_rarray(θ))
            for (l, ldev) in zip(post.lklhds, dev.lklhds)
                @test convert(Float64, @jit(logdensityof(ldev, mv))) ≈ logdensityof(l, Array(baseimage(mv))) rtol = 1.0e-12
            end
        end

        @testset "self-calibration" begin
            selfcal = InstrumentModel(Jones((x, c) -> complex(exp(x.lg)); lg = ArrayPrior(g0)); selfcal = true)
            stdnormal = p -> transport_to(p, Comrade.PT.StdNormal())
            check_matches_cpu(VLBIPosterior(vis; sky = static, instrument = selfcal, admode = nothing), ShardLayout(mesh; Node = :d); transport = stdnormal)
            check_matches_cpu(VLBIPosterior(vis; sky = planemodel(vis, 4), instrument = selfcal, admode = nothing), bytime; transport = stdnormal, rtol = 1.0e-9)
        end

        @testset "errors" begin
            @test_throws "already on a device" todevice(todevice(post_cpu, Reactant.Sharding.NoSharding()), bytime)
            @test_throws "`Ti` is not a dimension" todevice(VLBIPosterior(vis; sky = static, admode = nothing), bytime)
        end
    end
end
