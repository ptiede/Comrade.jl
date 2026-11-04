@testset "VLBIPosterior data products" begin
    _, vis, amp, lcamp, cphase = load_data()
    g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 32, 32)
    skym = SkyModel(test_model, test_prior(), g)
    G = Jones(
        (x, c) -> exp(complex(x.lg, x.gp));
        lg = ArrayPrior(IIDStationPrior(ScanSeg(), VLBIGaussian(0.0, 1.0))),
        gp = ArrayPrior(IIDStationPrior(ScanSeg(), VLBIGaussian(0.0, 1.0)), refant = SensitiveReference(0.0)),
    )
    intm = InstrumentModel(G)
    fd = ImgNormalData(Comrade.flux, 1.0, 0.1)

    @test Comrade.datadomain(vis) === Comrade.VisDomain()
    @test Comrade.datadomain(cphase) === Comrade.VisDomain()
    @test Comrade.datadomain(fd) === Comrade.ImageDomain()
    @test (@inferred Comrade._splitdomains((vis, fd, lcamp))) === ((vis, lcamp), (fd,))

    @testset "single field" begin
        post = VLBIPosterior(lcamp, fd, cphase; sky = skym)
        @test dataproducts(post) == (lcamp, cphase)
        @test post.dataimg == (fd,)
        @test Comrade.datatype(typeof(post)) isa Comrade.DualData
        @test Comrade.datatype(typeof(VLBIPosterior(lcamp, cphase; sky = skym))) isa Comrade.VisData
        x = prior_sample(post)
        img = intensitymap(skymodel(post, x), g)
        lvis = logdensityof(VLBIPosterior(lcamp, cphase; sky = skym), x)
        @test logdensityof(post, x) ≈ lvis + logdensityof(Comrade.makelikelihood(fd), img)
        @inferred logdensityof(post, x)
        @test contains(sprint(show, MIME"text/plain"(), post), "Image Data Products: ImgNormalData")
    end

    @testset "several fields" begin
        sky = (; A = skym, B = skym)
        post = VLBIPosterior((; B = amp, A = (vis, fd)); sky, instrument = intm)
        @test keys(dataproducts(post)) == (:A, :B)
        @test dataproducts(post).B == (amp,)
        @test post.dataimg == (; A = (fd,), B = ())
        x = prior_sample(post)
        img = intensitymap(skymodel(post, x).A, g)
        lvis = logdensityof(VLBIPosterior((; A = vis, B = amp); sky, instrument = intm), x)
        @test logdensityof(post, x) ≈ lvis + logdensityof(Comrade.makelikelihood(fd), img)
        @inferred logdensityof(post, x)
        @test isnothing(VLBIPosterior((; A = vis, B = amp); sky).dataimg)
    end

    @testset "errors" begin
        @test_throws "the sky models and the data need the same field names" VLBIPosterior((; A = vis); sky = (; A = skym, B = skym))
        @test_throws "Float64 is not a data product" VLBIPosterior(vis, 1.0; sky = skym)
        @test_throws "a posterior needs a visibility-domain data product" VLBIPosterior(fd; sky = skym)
        @test_throws "the field B needs a visibility-domain data product" VLBIPosterior((; A = vis, B = fd); sky = (; A = skym, B = skym))
        @test_throws "takes one NamedTuple of data products" VLBIPosterior(vis; sky = (; A = skym))
        @test_throws "takes a NamedTuple of sky models" VLBIPosterior((; A = vis); sky = skym)
        @test_throws "the sky model of the field B must be an AbstractSkyModel" VLBIPosterior((; A = vis, B = amp); sky = (; A = skym, B = 1))
        @test_throws "the instrument must be an AbstractInstrumentModel" VLBIPosterior(vis; sky = skym, instrument = G)
        @test_throws "the sky must be an AbstractSkyModel" VLBIPosterior(vis; sky = test_model)
        @test_throws UndefKeywordError VLBIPosterior(skym, intm, vis)
        @test !isdefined(Comrade, :Field)
        @test !isdefined(Comrade, :ParamLogDensity)
    end
end
