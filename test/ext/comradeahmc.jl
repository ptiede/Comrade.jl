using AdvancedHMC
using Enzyme


@testset "ComradeAdvancedHMCExt" begin

    _, _, _, lcamp, cphase = load_data()
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256)
    skym = SkyModel(test_model, test_prior(), g)
    post = VLBIPosterior(skym, lcamp, cphase; admode = set_runtime_activity(Enzyme.Reverse))

    x0 = (
        sky = (
            f1 = 1.0916271439905998,
            σ1 = 8.230088139590025e-11,
            τ1 = 0.49840994275315254,
            ξ1 = -1.0489388962890198,
            f2 = 0.553944311447593,
            σ2 = 4.1283218512580634e-11,
            τ2 = 0.5076731020106428,
            ξ2 = -0.5376269092893298,
            x = 1.451956089157719e-10,
            y = 1.455983181049137e-10,
        ),
    )
    s1 = NUTS(0.65)
    hchain = sample(post, s1, 1_000; n_adapts = 500, progress = false)
    hchain = sample(post, s1, 1_000; n_adapts = 500, progress = false, initial_params = x0)
    out = sample(post, s1, 1_000; n_adapts = 500, saveto = DiskStore(name = joinpath(@__DIR__, "Test")), initial_params = x0)
    out = sample(post, s1, 1_200; n_adapts = 500, saveto = DiskStore(name = joinpath(@__DIR__, "Test")), initial_params = x0, restart = true)

    # @test logdensityof(post, hchain[1]) ≈ logdensityof(post, x0)
    # @test logdensityof(post, out[1]) ≈ logdensityof(post, x0)

    show(IOBuffer(), MIME"text/plain"(), hchain)
    hchain.sky

    c1 = load_samples(out)
    @test c1[201:451] == load_samples(out, 201:451)

    c1 = load_samples(out; table = :stats)
    @test c1[1:451] == load_samples(out, 1:451; table = :stats)

    c1 = load_samples(joinpath(@__DIR__, "Test"))

    sample(post, s1, 1_000; n_adapts = 500, saveto = DiskStore(name = joinpath(@__DIR__, "Test")), initial_params = x0, restart = true)

    @test c1[201:451] == load_samples(joinpath(@__DIR__, "Test"), 201:451)


    rm(joinpath(@__DIR__, "Test"), recursive = true)

    # Regression: sampling to a DiskStore *without* `initial_params` used to pre-sample the
    # prior into a flat vector and then double-transform it through `inverse`, which errored.
    # It should now run straight through (prior start chosen inside `initialize_params`).
    noinit_dir = joinpath(@__DIR__, "TestNoInit")
    out_noinit = sample(post, s1, 200; n_adapts = 100, saveto = DiskStore(name = noinit_dir, stride = 100))
    @test length(load_samples(out_noinit)) == 200
    rm(noinit_dir, recursive = true)

    # The DiskStore callback fires once per batch and receives the unified `info` NamedTuple
    # (shared with the ReactantNUTS disk path). Check the guaranteed fields are populated.
    seen = NamedTuple[]
    cb = info -> begin
        push!(
            seen, (;
                info.round, info.nrounds, info.num_samples,
                istime = info.time isa Real, isstep = info.step_size isa Real,
                ndiv = count(info.numerical_error), nerr = length(info.numerical_error),
                nstats = length(info.extras.stats), hassky = haskey(info.params, :sky),
            )
        )
        return nothing
    end
    cb_dir = joinpath(@__DIR__, "TestCB")
    sample(post, s1, 200; n_adapts = 100, saveto = DiskStore(name = cb_dir, stride = 100, callback = cb), initial_params = x0)
    @test length(seen) == 2
    @test [s.round for s in seen] == [1, 2]
    @test all(s -> s.nrounds == 2, seen)
    @test all(s -> s.num_samples == 100, seen)
    @test all(s -> s.nerr == 100, seen)
    @test all(s -> s.nstats == 100, seen)
    @test all(s -> s.istime && s.isstep && s.hassky, seen)
    rm(cb_dir, recursive = true)
end
