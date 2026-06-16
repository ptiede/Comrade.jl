using Reactant
using Random
using Serialization

const ReactantEx = Comrade.ComradeBase.ReactantEx

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
