using ChainRulesTestUtils
using ChainRulesCore


function testmodel(m::Comrade.AbstractModel, npix=1024, atol=1e-4)
    plot(m)
    img = intensitymap(m, 2*Comrade.radialextent(m), 2*Comrade.radialextent(m), npix, npix)
    plot(img)
    img2 = similar(img)
    intensitymap!(img2, m)
    @test eltype(img) === Float64
    @test isapprox(flux(m), flux(img), atol=atol)
    @test isapprox(mean(img .- img2), 0, atol=1e-8)
    cache = Comrade.create_cache(Comrade.FFTAlg(padfac=4), img./flux(img)*flux(m))
    u = fftshift(fftfreq(size(img,1), 1/img.psizex))./30
    Plots.closeall()
    @test isapprox(maximum(abs, (visibility.(Ref(m), u', u) .- cache.sitp.(u', u))), 0.0, atol=atol*10)
    img = nothing
    img2 =nothing
    cache = nothing
    u = nothing
    GC.gc()
end

function testft(m, npix=256, atol=1e-4)
    mn = Comrade.NonAnalyticTest(m)
    uu = 0.25*randn(1000)
    vv = 0.25*randn(1000)
    img = intensitymap(m, 2*Comrade.radialextent(m), 2*Comrade.radialextent(m), npix, npix; pulse=DeltaPulse())
    mimg_ff = modelimage(mn, img, FFTAlg(padfac=4))
    mimg_nf = modelimage(mn, img, NFFTAlg())
    mimg_df = modelimage(mn, img, DFTAlg())

    va = visibilities(m, uu, vv)
    vff = visibilities(mimg_ff, uu, vv)
    vnf = visibilities(mimg_nf, uu, vv)
    vdf = visibilities(mimg_df, uu, vv)

    @test isapprox(maximum(abs, va-vff), 0, atol=atol*5)
    @test isapprox(maximum(abs, va-vnf), 0, atol=atol)
    @test isapprox(maximum(abs, va-vdf), 0, atol=atol)
    img = nothing
    mimg_ff = nothing
    mimg_nf = nothing
    mimg_df = nothing
    GC.gc()
end

@testset "Moments" begin
    img = IntensityMap(zeros(512, 512), 30.0, 30.0)
    m1 = Gaussian()
    intensitymap!(img, m1)
    @test isapprox(centroid(img)[1], 0.0, atol=1e-5)
    @test isapprox(centroid(img)[2], 0.0, atol=1e-5)

    I = inertia(img)
    I2 = inertia(img; center=false)
    @test isapprox(I, [1.0 0.0; 0.0 1.0], atol=1e-5)
    @test I ≈ I2

    m2 = shifted(m1, 1.0, 1.0)
    intensitymap!(img, m2)
    @test isapprox(centroid(img)[1], 1.0, atol=1e-5)
    @test isapprox(centroid(img)[2], 1.0, atol=1e-5)
    @test isapprox(inertia(img), I, atol=1e-5)

    m3 = stretched(m1, 2.0, 1.0)
    intensitymap!(img, m3)
    @test isapprox(centroid(img)[1], 0.0, atol=1e-5)
    @test isapprox(centroid(img)[2], 0.0, atol=1e-5)
    I3 = inertia(img)
    @test isapprox(I3, [4.0 0.0; 0.0 1.0], atol=1e-5)

end


@testset "FFTTest" begin
    @testset "Base" begin
        m = Gaussian()
        testft(m)
    end

    @testset "Mod" begin
        m = rotated(stretched(Gaussian(), 0.5, 1.0), π/3)
        testft(m)
        ms = shifted(m, 1.0,1.0)
        testft(ms)
    end

    @testset "Add" begin
        m1 = rotated(stretched(Gaussian(), 0.5, 1.0), π/3) + shifted(Gaussian(), 1.0, 1.0)
        testft(m1)
    end
end

@testset "Primitive models" begin

    @testset "Gaussian" begin
        m = Gaussian()
        testmodel(m, 1024, 1e-5)
    end

    @testset "Disk" begin
        m = smoothed(Disk(), 0.25)
        ComradeBase.intensity_point(Disk(), 0.0, 0.0)
        testmodel(m)
    end

    @testset "Ring" begin
        m = smoothed(Ring(), 0.25)
        ComradeBase.intensity_point(Ring(), 0.0, 0.0)
        testmodel(m, 2048)
    end

    @testset "ParabolicSegment" begin
        m = ParabolicSegment()
        m2 = ParabolicSegment(2.0, 2.0)
        @test stretched(m, 2.0, 2.0) == m2
        @test ComradeBase.intensity_point(m, 0.0, 1.0) != 0.0
        testmodel(m, 2424, 1e-3)
    end


    @testset "MRing1" begin
        α = [0.25,]
        β = [0.1,]
        test_rrule(Comrade.visibility_point, MRing(α, β), 0.5, 0.25)
        # We convolve it to remove some pixel effects
        m = convolved(MRing(α, β), stretched(Gaussian(), 0.1, 0.1))
        m2 = convolved(MRing(α[1], β[1]), stretched(Gaussian(), 0.1, 0.1))
        @test visibility(m, 0.1, 0.1) == visibility(m2, 0.1, 0.1)
        testmodel(m, 2048, 1e-3)
    end

    @testset "MRing2" begin
        α = (0.25, -0.1)
        β = (0.1, 0.2)
        test_rrule(Comrade.visibility_point, MRing(α, β), 0.5, 0.25)

        # We convolve it to remove some pixel effects
        m = convolved(MRing(α, β), stretched(Gaussian(), 0.1, 0.1))
        testmodel(m, 2048, 1e-3)
    end


    @testset "ConcordanceCrescent" begin
        m = ConcordanceCrescent(20.0, 10.0, 5.0, 0.5)
        testmodel(m)
    end


    @testset "Crescent" begin
        m = smoothed(Crescent(5.0, 2.0, 1.0, 0.5), 1.0)
        testmodel(m,1024,1e-3)
    end

    @testset "ExtendedRing" begin
        mr = ExtendedRing(8.0)
        rad = 2.5*Comrade.radialextent(mr)
        m = modelimage(mr, IntensityMap(zeros(1024,1024), rad, rad), Comrade.FFTAlg(padfac=4))
        testmodel(m)
    end
end

@testset "ModelImage" begin
    m1 = Gaussian()
    m2 = ExtendedRing(10.0)
    mimg1 = modelimage(m1)
    mimg2 = modelimage(m2)

    show(mimg1)

    img = similar(mimg2.image)
    intensitymap!(img, m2)
    @test m1 == mimg1
    @test isapprox(mean(img .- mimg2.image), 0.0, atol=1e-8)
end



@testset "Modifiers" begin
    ma = Gaussian()
    mb = ExtendedRing(8.0)
    @testset "Shifted" begin
        mas = shifted(ma, 0.5, 0.5)
        mbs = shifted(mb, 0.5, 0.5)
        testmodel(mas)
        testmodel(modelimage(mbs, IntensityMap(zeros(1024, 1024),
                                               2*Comrade.radialextent(mbs),
                                               2*Comrade.radialextent(mbs))))
    end

    @testset "Renormed" begin
        m1 = 3.0*ma
        m2 = ma*3.0
        m2inv = ma/(1/3)
        @test visibility(m1, 4.0, 0.0) == visibility(m2, 4.0, 0.0)
        @test visibility(m2, 4.0, 0.0) == visibility(m2inv, 4.0, 0.0)
        mbs = 3.0*mb
        testmodel(m1)
        testmodel(modelimage(mbs, IntensityMap(zeros(1024, 1024),
                                               2.5*Comrade.radialextent(mbs),
                                               2.5*Comrade.radialextent(mbs))))
    end

    @testset "Stretched" begin
        mas = stretched(ma, 5.0, 4.0)
        mbs = stretched(mb, 5.0, 4.0)
        testmodel(mas)
        testmodel(modelimage(mbs, IntensityMap(zeros(2024, 2024),
                                               2*Comrade.radialextent(mbs),
                                               2*Comrade.radialextent(mbs))), 1024, 1e-3)
    end

    @testset "Rotated" begin
        mas = rotated(ma, π/3)
        mbs = rotated(mb, π/3)
        testmodel(mas)
        testmodel(modelimage(mbs, IntensityMap(zeros(1024, 1024),
                                               2*Comrade.radialextent(mbs),
                                               2*Comrade.radialextent(mbs))))
    end

    @testset "AllMods" begin
        mas = rotated(stretched(shifted(ma, 0.5, 0.5), 5.0, 4.0), π/3)
        mbs = rotated(stretched(shifted(mb, 0.5, 0.5), 5.0, 4.0), π/3)
        testmodel(mas)
        testmodel(modelimage(mbs, IntensityMap(zeros(2024, 2024),
                                               2*Comrade.radialextent(mbs),
                                               2*Comrade.radialextent(mbs))), 1024, 1e-3)
    end
end

@testset "CompositeModels" begin
    m1 = Gaussian()
    m2 = ExtendedRing(8.0)

    @testset "Add models" begin
        img = IntensityMap(zeros(1024, 1024),
                                        20.0,
                                        20.0)
        mt1 = m1 + m2
        mt2 = shifted(m1, 1.0, 1.0) + m2
        mt3 = shifted(m1, 1.0, 1.0) + 0.5*stretched(m2, 0.9, 0.8)
        mc = Comrade.components(mt1)
        @test mc[1] === m1
        @test mc[2] === m2
        @test flux(mt1) ≈ flux(m1) + flux(m2)

        testmodel(modelimage(mt1, img))
        testmodel(modelimage(mt2, img))
        testmodel(modelimage(mt3, img))
    end

    @testset "Convolved models" begin
        img = IntensityMap(zeros(1024, 1024),
                                        20.0,
                                        20.0)
        mt1 = convolved(m1, m2)
        mt2 = convolved(shifted(m1, 1.0, 1.0), m2)
        mt3 = convolved(shifted(m1, 1.0, 1.0), 0.5*stretched(m2, 0.9, 0.8))
        mc = Comrade.components(mt1)
        @test mc[1] === m1
        @test mc[2] === m2

        testmodel(modelimage(mt1, img))
        testmodel(modelimage(mt2, img))
        testmodel(modelimage(mt3, img))
    end

    @testset "All composite" begin
        img = IntensityMap(zeros(1024, 1024),
                                            20.0,
                                            20.0)

        mt = m1 + convolved(m1, m2)
        mc = Comrade.components(mt)
        @test mc[1] === m1
        @test mc[2] === m1
        @test mc[3] === m2

        testmodel(modelimage(mt, img))

    end
end

@testset "PolarizedModel" begin
    mI = stretched(MRing((0.2,), (0.1,)), 20.0, 20.0)
    mQ = 0.2*stretched(MRing((0.0,), (0.6,)), 20.0, 20.0)
    mU = 0.2*stretched(MRing((0.1,), (-0.6,)), 20.0, 20.0)
    mV = 0.0*stretched(MRing((0.0,), (-0.6,)), 20.0, 20.0)
    m = PolarizedModel(mI, mQ, mU, mV)

    v = coherencymatrix(m, 0.005, 0.01)
    @test evpa(v) == evpa(m, 0.005, 0.01)
    @test m̆(v) == m̆(m, 0.005, 0.01)

    I = IntensityMap(zeros(1024,1024), 100.0, 100.0)
    Q = similar(I)
    U = similar(I)
    V = similar(I)
    pimg1 = IntensityMap(I,Q,U,V)
    intensitymap!(pimg1, m)
    pimg2 = intensitymap(m, 100.0, 100.0, 1024, 1024)
    @test isapprox(sum(abs, (stokes(pimg1, :I) .- stokes(pimg2, :I))), 0.0, atol=1e-12)
    @test isapprox(sum(abs, (stokes(pimg1, :Q) .- stokes(pimg2, :Q))), 0.0, atol=1e-12)
    @test isapprox(sum(abs, (stokes(pimg1, :U) .- stokes(pimg2, :U))), 0.0, atol=1e-12)
    @test isapprox(sum(abs, (stokes(pimg1, :V) .- stokes(pimg2, :V))), 0.0, atol=1e-12)

end

@testset "Image SqExp" begin
   c = intensitymap(rotated(stretched(Gaussian(), 2.0, 1.0), π/8), 12.0, 12.0, 12, 12; pulse=SqExpPulse(3.0))
   #mI = DImage(c.im, SqExpPulse(5.0))
   testmodel(modelimage(c, FFTAlg(padfac=3)), 1024, 1e-3)
end
#@testset "DImage Bspline0" begin
#   mI = DImage(rand(8,8), BSplinePulse{0}())
#   testmodel(mI, 1e-2)
#end
@testset "DImage BSpline1" begin
    c = intensitymap(rotated(stretched(Gaussian(), 2.0, 1.0), π/8), 12.0, 12.0, 12, 12; pulse=BSplinePulse{1}())
    testmodel(modelimage(c, FFTAlg(padfac=3)), 1024, 1e-3)
end

@testset "DImage BSpline3" begin
    c = intensitymap(rotated(stretched(Gaussian(), 2.0, 1.0), π/8), 12.0, 12.0, 12, 12; pulse=BSplinePulse{3}())
    testmodel(modelimage(c, FFTAlg(padfac=3)), 1024, 1e-3)
end

@testset "modelimage cache" begin
    img = intensitymap(rotated(stretched(Gaussian(), μas2rad(2.0), μas2rad(1.0)), π/8),
                       μas2rad(12.0), μas2rad(12.0), 24, 12; pulse=BSplinePulse{3}())
    _,_, amp, lcamp, cphase = load_data()

    cache_nf = create_cache(NFFTAlg(amp), img)
    cache_df = create_cache(DFTAlg(amp), img)
    ac_amp = arrayconfig(amp)
    ac_lcamp = arrayconfig(lcamp)
    ac_cphase = arrayconfig(cphase)

    mimg_nf = modelimage(img, cache_nf)
    mimg_df = modelimage(img, cache_df)

    vnf = visibilities(mimg_nf, ac_amp)
    vdf = visibilities(mimg_df, ac_amp)

    atol = 1e-5

    @test isapprox(maximum(abs, vnf-vdf), 0, atol=atol)


    cpnf = closure_phases(mimg_nf, ac_cphase)
    cpdf = closure_phases(mimg_df, ac_cphase)

    @test isapprox(maximum(abs, cis.(cpnf-cpdf) .- 1.0 ), 0, atol=atol)

    lcnf = logclosure_amplitudes(mimg_nf, ac_lcamp)
    lcdf = logclosure_amplitudes(mimg_df, ac_lcamp)

    @test isapprox(maximum(abs, lcnf-lcdf), 0, atol=atol)



    @testset "nuft pullback" begin
        test_rrule(Comrade.nuft, cache_nf.plan ⊢ NoTangent(), img.im')
    end
end
