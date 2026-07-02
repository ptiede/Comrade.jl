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

    # flat path (the gradient path NUTS uses): TV.transform_with
    let t = asflat(pcd)
        x = rand(TV.dimension(t))
        y, _, _ = TV.transform_with(TV.LogJac(), t, x, firstindex(x))
        @test y[fixed_index] == fixed_values
        f(xx) = sum(first(TV.transform_with(TV.LogJac(), t, xx, 1)))
        @test convert(Float64, @jit f(Reactant.to_rarray(x))) ≈ sum(y)
    end

    # cube path: HypercubeTransform._step_transform must still place the fixed
    # values correctly on the CPU. (The inner cube transform itself is not yet
    # Reactant-traceable, independent of this fix, so only the flat path is jit'd.)
    let t = ascube(pcd)
        u = rand(TV.dimension(t))
        y, _ = HypercubeTransform._step_transform(t, u, firstindex(u))
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

    rm(dir, recursive = true)
end

@testset "Sharded VLBIPosterior" begin
    Sharding = Reactant.Sharding
    ndev = length(Reactant.devices())
    isifrt = Reactant.XLA.REACTANT_XLA_RUNTIME == "IFRT"

    # Build a balanced multi-frequency dataset by cycling `nfr` frequencies onto the
    # visibility points. `nfr ≥ 2` always exercises a 3D (X,Y,Fr) analytic grid, whose
    # `AnalyticAlg` must survive `prepare_device` (forcing it to a NUFFT alg would plan
    # over the now device-resident vis coords and error). When ≥2 devices are present we
    # additionally shard the Fr blocks across them.
    _, vis, _, _, _ = load_data()
    g2 = imagepixels(μas2rad(150.0), μas2rad(150.0), 32, 32)
    nfr = max(ndev, 2)
    frs = collect(range(220.0e9, 240.0e9; length = nfr))
    vismf = deepcopy(vis)
    n = length(vismf)
    arrayconfig(vismf)[:Fr] .= [frs[(i % nfr) + 1] for i in 0:(n - 1)]

    # Plain (Serial) CPU grid — sharding is requested later via `prepare_device`, never
    # baked into the grid (that would force the CPU posterior down the Reactant path).
    g3 = RectiGrid((; X = g2.X, Y = g2.Y, Fr = frs))
    skym3 = SkyModel(test_model, test_prior(), g3)

    # Step 1 regroups the data into contiguous Fr blocks at construction.
    post_cpu = VLBIPosterior(skym3, vismf; admode = nothing)
    @test issorted(arrayconfig(dataproducts(post_cpu)[begin])[:Fr])

    tpost_cpu = asflat(post_cpu)
    x0 = prior_sample(Random.default_rng(), tpost_cpu)
    ld_serial = logdensityof(tpost_cpu, x0)

    # Unsharded regression: an unsharded `ReactantEx()` over the 3D analytic grid must
    # still match the serial logdensity (guards the `AnalyticAlg`-preserve fix).
    post_r = Comrade.prepare_device(post_cpu, ReactantEx())
    tpost_r = asflat(post_r)
    ldr = @jit logdensityof(tpost_r, Reactant.to_rarray(x0))
    @test convert(Float64, ldr) ≈ ld_serial

    # Sharding requires ≥2 IFRT devices; skip the device-specific asserts otherwise so
    # the suite stays green on a single-device CI runner.
    if isifrt && ndev >= 2
        mesh = Comrade.ComradeBase.shardmesh(ndev; names = (:dev,))
        ex = ReactantEx(mesh; Fr = :dev)
        post_sh = Comrade.prepare_device(post_cpu, ex)

        # The measurement (`obs`) and noise (`kernel.S`) vectors land sharded on the mesh.
        @test Sharding.is_sharded(post_sh.lklhds[1].obs)
        @test Sharding.is_sharded(post_sh.lklhds[1].kernel.S)

        # The stored data product is sharded coherently with its likelihood (same layout as `obs`),
        # not left as an unsharded `to_rarray` copy.
        @test Sharding.is_sharded(measurement(dataproducts(post_sh)[begin]))
        @test Sharding.is_sharded(noise(dataproducts(post_sh)[begin]))

        # Sharding is layout-only: the logdensity must be unchanged.
        tpost_sh = asflat(post_sh)
        ldsh = @jit logdensityof(tpost_sh, Reactant.to_rarray(x0))
        @test convert(Float64, ldsh) ≈ ld_serial
    else
        @info "Sharded VLBIPosterior: skipping multi-device asserts" ndev isifrt
    end
end

@testset "Sharded VLBIPosterior (coherency)" begin
    Sharding = Reactant.Sharding
    ndev = length(Reactant.devices())
    isifrt = Reactant.XLA.REACTANT_XLA_RUNTIME == "IFRT"

    # Coherency data exercises the polarized sharding path: `kernel.S` is a vector of 2x2 noise
    # matrices and `obs` a vector of coherency matrices, sharded on the same flat dim-1 blocks as the
    # scalar-visibility case. Cycle `nfr` frequencies onto the coherency points and lay a polarized
    # analytic sky over an (X,Y,Fr) grid (analytic => `AnalyticAlg` must survive `prepare_device`).
    _, _, _, _, _, dcoh = load_data()
    g2 = imagepixels(μas2rad(150.0), μas2rad(150.0), 32, 32)
    nfr = max(ndev, 2)
    frs = collect(range(220.0e9, 240.0e9; length = nfr))
    dcohmf = deepcopy(dcoh)
    n = length(dcohmf)
    arrayconfig(dcohmf)[:Fr] .= [frs[(i % nfr) + 1] for i in 0:(n - 1)]

    g3 = RectiGrid((; X = g2.X, Y = g2.Y, Fr = frs))
    skym3 = SkyModel(test_skymodel_polarized, test_prior(), g3; metadata = (; lp = 0.1))
    # coherency data needs an instrument model (feed basis) to turn Stokes into coherency matrices
    intm = InstrumentModel(JonesR())

    post_cpu = VLBIPosterior(skym3, intm, dcohmf; admode = nothing)
    @test issorted(arrayconfig(dataproducts(post_cpu)[begin])[:Fr])

    tpost_cpu = asflat(post_cpu)
    x0 = prior_sample(Random.default_rng(), tpost_cpu)
    ld_serial = logdensityof(tpost_cpu, x0)

    # Unsharded regression over the coherency likelihood.
    post_r = Comrade.prepare_device(post_cpu, ReactantEx())
    tpost_r = asflat(post_r)
    ldr = @jit logdensityof(tpost_r, Reactant.to_rarray(x0))
    @test convert(Float64, ldr) ≈ ld_serial

    if isifrt && ndev >= 2
        mesh = Comrade.ComradeBase.shardmesh(ndev; names = (:dev,))
        ex = ReactantEx(mesh; Fr = :dev)
        post_sh = Comrade.prepare_device(post_cpu, ex)

        # The coherency `obs`/`kernel.S` (vectors of static matrices) shard onto the mesh, and the
        # stored data product shards coherently with them.
        @test Sharding.is_sharded(post_sh.lklhds[1].obs)
        @test Sharding.is_sharded(post_sh.lklhds[1].kernel.S)
        @test Sharding.is_sharded(measurement(dataproducts(post_sh)[begin]))

        tpost_sh = asflat(post_sh)
        ldsh = @jit logdensityof(tpost_sh, Reactant.to_rarray(x0))
        @test convert(Float64, ldsh) ≈ ld_serial
    else
        @info "Sharded VLBIPosterior (coherency): skipping multi-device asserts" ndev isifrt
    end
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
