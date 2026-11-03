using ChainRulesTestUtils
using ChainRulesCore
using FiniteDifferences
using Enzyme
using PythonCall
using FFTW
using StructArrays
using StaticArrays
using Distributions
using VLBIImagePriors
using StaticArrays
using Tables
import TransformVariables as TV
using VLBIImagePriors
using JET

@testset "SkyModel" begin

    f = test_model
    g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 256, 256)
    skym = SkyModel(f, test_prior(), g)
    show(IOBuffer(), MIME"text/plain"(), skym)

    x = rand(Comrade.NamedDist(test_prior()))
    m = Comrade.skymodel(skym, x)
    skyf = FixedSkyModel(m, g)

    @testset "ObservedSkyModel" begin
        _, vis, amp, lcamp, cphase = load_data()

        oskym, = Comrade.observe(skym, vis)
        oskyf, = Comrade.observe(skyf, vis)

        @test Comrade.skymodel(oskym, x) == m
        dmimg, dmvis = Comrade.idealmaps(Comrade.VisData(), oskym, (; sky = x))
        dfimg, dfvis = Comrade.idealmaps(Comrade.VisData(), oskyf, (; sky = x))
        @test dmimg ≈ dfimg
        @test dmvis ≈ dfvis
    end

    @testset "MultiSkyModel" begin
        _, vis, amp, lcamp, cphase = load_data()
        skytot = MultiSkyModel((dynamic = skym, static = skyf))
        oskytot, ptot = Comrade.observe(skytot, vis)

        show(IOBuffer(), MIME"text/plain"(), skytot)


        x = rand(ptot)

        oskym, = Comrade.observe(skym, vis)
        oskyf, = Comrade.observe(skyf, vis)

        tt = Comrade.skymodel(oskytot, x)
        @test tt.dynamic == Comrade.skymodel(oskym, x.dynamic)
        @test tt.static == Comrade.skymodel(oskyf, x.static)

        vtot = last(Comrade.idealmaps(Comrade.VisData(), oskytot, (; sky = x)))
        vdyn = last(Comrade.idealmaps(Comrade.VisData(), oskym, (; sky = x.dynamic)))
        vstat = last(Comrade.idealmaps(Comrade.VisData(), oskyf, (; sky = x.static)))

        @test_opt Comrade.idealmaps(Comrade.VisData(), oskytot, (; sky = x))

        @test vtot ≈ vdyn + vstat

        post = VLBIPosterior(skytot, vis)
        @test loglikelihood(post, (; sky = x)) ≈ logdensityof(Comrade.makelikelihood(vis), vtot)
    end
end

@testset "GMRF" begin
    δ = rand(GMRF(64.0, (128, 128)))
    m = modify(Gaussian(), Stretch(μas2rad(40.0)))
    g = spatialgrid(μas2rad(150.0), μas2rad(150.0), 128, 128)
    img = apply_fluctuations(exp, m, g, δ)
    mimg = intensitymap(m, g)
    img2 = apply_fluctuations(exp, mimg, δ)
    @test img isa IntensityMap
    @test img2 ≈ img

    img3 = apply_fluctuations(mimg, δ)

    img4 = apply_fluctuations(CenteredLR(), m, g, δ)
    img5nf = apply_fluctuations(CenteredLR(), mimg, δ)
    img5 = apply_fluctuations(CenteredLR(), mimg ./ Comrade.flux(mimg), δ)
    @test img4 ≈ img5
    @test img5nf ≈ img5

    img6 = apply_fluctuations(UnitFluxMap(exp), mimg ./ Comrade.flux(mimg), δ)
    @test img6 ≈ img5
    img7 = apply_fluctuations(UnitFluxMap(exp), mimg, δ)
    @test img7 ≈ img5

    for f in (exp, UnitFluxMap(exp), CenteredLR())
        out = similar(mimg)
        @test apply_fluctuations!(f, out, mimg, δ) === out
        @test out ≈ apply_fluctuations(f, mimg, δ)
    end
end

@testset "InstrumentModel" begin
    @testset "Hypercube Instrument" begin
        _, dvis, amp, lcamp, cphase, dcoh = load_data()
        gmap = (x, c) -> exp(x.lg + 1im .* x.gp)

        @testset "standard" begin
            intprior = (
                lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
                gp = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, inv(π^2))); refant = SensitiveReference(0.0)),
            )

            intm = InstrumentModel(Jones(gmap; intprior...))
            skym = SkyModel(test_model, test_prior(), spatialgrid(μas2rad(150.0), μas2rad(150.0), 256, 256))
            post = VLBIPosterior(skym, intm, dvis)

            cpost = ascube(post)
            x = prior_sample(cpost)
            p = Comrade.transform(cpost, x)
            @test all(x -> 0 < x < 1, x)

            l1 = logdensityof(cpost, x)
            l2 = Comrade.loglikelihood(post, p)

            @test l1 ≈ l2
        end


    end
end
