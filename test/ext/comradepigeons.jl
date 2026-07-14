using Pigeons
using Enzyme

@testset "ComradePigeonsExt" begin
    _, _, _, lcamp, cphase = load_data()
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 32, 32)
    skym = SkyModel(test_model, test_prior(), g)
    post = VLBIPosterior(skym, lcamp, cphase; admode = set_runtime_activity(Enzyme.Reverse))

    PigeonsExt = Base.get_extension(Comrade, :ComradePigeonsExt)

    # cube (`ascube`/StdUniform) path: bounded reference, so the default explorer is a slice
    # sampler and the reference is the (transport-free) StdUniform density — 0 inside the unit
    # hypercube, -Inf outside.
    @testset "cube path" begin
        cpost = ascube(post)
        @test Pigeons.default_explorer(cpost) isa SliceSampler
        ref = Pigeons.default_reference(cpost)
        @test ref isa PigeonsExt.PriorRef
        @test ref(fill(0.5, dimension(cpost))) == 0
        @test ref(fill(1.5, dimension(cpost))) == -Inf
        pt = pigeons(target = cpost, record = [traces], n_chains = 2, n_rounds = 2)
        chain = sample_array(cpost, pt)
        @test length(chain) > 0
    end

    # flat (`asflat`/TVFlat) path: unbounded, so the default explorer is gradient-based AutoMALA
    # (Enzyme). Running it exercises the general (non-cube) `sample_iid!`, the pulled-back prior
    # `PriorRef`, and the Enzyme `ADgradient`/`logdensity_and_gradient` glue.
    @testset "flat path" begin
        fpost = asflat(post)
        @test Pigeons.default_explorer(fpost) isa Pigeons.AutoMALA
        @test isfinite(Pigeons.default_reference(fpost)(prior_sample(fpost)))
        pt = pigeons(target = fpost, record = [traces], n_chains = 2, n_rounds = 2)
        @test length(sample_array(fpost, pt)) > 0
    end

    # standard-normal (`StdNormal`) path: shares the general non-cube explorer/`sample_iid!`
    # methods. Regression guard — a StdNormal target used to hit a `MethodError` in `sample_iid!`.
    @testset "StdNormal path" begin
        npost = transport_to(post, StdNormal())
        @test Pigeons.default_explorer(npost) isa Pigeons.AutoMALA
        @test Pigeons.default_reference(npost) isa PigeonsExt.PriorRef
        pt = pigeons(target = npost, record = [traces], n_chains = 2, n_rounds = 2)
        @test length(sample_array(npost, pt)) > 0
    end
end
