using Reactant
using Random
using Serialization

const ReactantEx = Comrade.ComradeBase.ReactantEx
const CRX = Base.get_extension(Comrade, :ComradeReactantExt)

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

    # Small warmup so the test stays in CI budget.
    s = ReactantNUTS(;
        n_adapts = 50, init_buffer = 10, term_buffer = 10, base_window = 15,
        max_tree_depth = 6,
    )

    # MemoryStore
    chain = sample(post, s, 100; chunk_size = 25)
    @test length(Comrade.postsamples(chain)) == 100
    @test hasproperty(samplerstats(chain), :numerical_error)
    @test haskey(samplerinfo(chain), :warmup_history)
    @test haskey(samplerinfo(chain), :sample_history)

    show(IOBuffer(), MIME"text/plain"(), chain)

    # DiskStore + restart
    dir = mktempdir()
    out = sample(post, s, 100; saveto = DiskStore(name = dir, stride = 25))
    @test out isa Comrade.DiskOutput
    @test out.nsamples == 100
    @test isfile(joinpath(dir, "state.jls"))
    @test isfile(joinpath(dir, "warmup_status.jls"))
    @test isfile(joinpath(dir, "metadata.jls"))

    # Warmup-status sidecar reports completion after a clean finish.
    let st = open(deserialize, joinpath(dir, "warmup_status.jls"))
        @test st.complete
        @test st.round_done == length(st.schedule)
    end

    # Metadata.jls captured the sample + warmup history.
    let meta = open(deserialize, joinpath(dir, "metadata.jls"))
        @test meta[:sampler] == :ReactantNUTS
        @test haskey(meta, :warmup_history)
        @test haskey(meta, :sample_history)
    end

    out2 = sample(post, s, 200; saveto = DiskStore(name = dir, stride = 25), restart = true)
    @test out2.nsamples == 200
    c = load_samples(out2)
    @test length(Comrade.postsamples(c)) == 200

    # Restart with a forged "incomplete warmup" status: should resume warmup
    # for the remaining windows and still produce a valid chain.
    dir2 = mktempdir()
    _ = sample(post, s, 50; saveto = DiskStore(name = dir2, stride = 25))
    let p = joinpath(dir2, "warmup_status.jls")
        st = open(deserialize, p)
        # Pretend only the first round finished.
        forged = (; st.schedule, round_done = 1, complete = false,
                  history = st.history[1:min(1, length(st.history))])
        open(io -> serialize(io, forged), p, "w")
    end
    out3 = sample(post, s, 100; saveto = DiskStore(name = dir2, stride = 25), restart = true)
    @test out3.nsamples == 100
    let st = open(deserialize, joinpath(dir2, "warmup_status.jls"))
        @test st.complete
    end

    # Schedule sanity
    let (ib, ws, tb) = CRX.warmup_schedule(s)
        @test ib + sum(ws) + tb == CRX.n_adapts(s)
    end

    rm(dir, recursive = true)
end
