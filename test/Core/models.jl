using ChainRulesTestUtils
using ChainRulesCore
using FiniteDifferences
using Zygote
using PythonCall
using FFTW


function testmodel(m::Comrade.AbstractModel, npix=1024, atol=1e-4)
    plot(m)
    g = imagepixels(4*Comrade.radialextent(m), 4*Comrade.radialextent(m), npix, npix)
    img = intensitymap(m, g)
    imgt = intensitymap(m, g, true)
    imgt2 = intensitymap(m, g, false)
    @test isapprox(maximum(img .- imgt), 0.0, atol=1e-8)
    @test isapprox(maximum(img .- imgt2), 0.0, atol=1e-8)
    plot(img)
    img2 = similar(img)
    intensitymap!(img2, m)
    @test eltype(img) === Float64
    @test isapprox(flux(m), flux(img), atol=atol)
    @test isapprox(maximum(parent(img) .- parent(img2)), 0, atol=1e-8)
    cache = Comrade.create_cache(Comrade.FFTAlg(padfac=3), img/flux(img)*flux(m))
    dx, dy = pixelsizes(img)
    u = fftshift(fftfreq(size(img,1), 1/dx))./30
    Plots.closeall()
    @test isapprox(maximum(abs, (visibility.(Ref(m), NamedTuple{(:U, :V)}.(u', u)) .- cache.sitp.(u', u))), 0.0, atol=atol*10)
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
    img = intensitymap(m, 2*Comrade.radialextent(m), 2*Comrade.radialextent(m), npix, npix)
    mimg_ff = modelimage(mn, zero(img), FFTAlg(padfac=4))
    mimg_nf = modelimage(mn, zero(img), NFFTAlg())
    mimg_df = modelimage(mn, zero(img), DFTAlg())
    cache = create_cache(FFTAlg(padfac=4), zero(img))
    cache_nf = create_cache(NFFTAlg(), zero(img))
    mimg_ff2 = modelimage(mn, cache)

    p = (U=uu, V=vv)
    va = visibilities(m, p)
    vff = visibilities(mimg_ff, p)
    vff2 = visibilities(mimg_ff2, p)
    vnf = visibilities(mimg_nf, p)
    vdf = visibilities(mimg_df, p)
    visibilities(modelimage(mn, cache_nf), p)

    @test isapprox(maximum(abs, vff2-vff), 0, atol=atol)
    @test isapprox(maximum(abs, va-vff), 0, atol=atol*5)
    @test isapprox(maximum(abs, va-vnf), 0, atol=atol)
    @test isapprox(maximum(abs, va-vdf), 0, atol=atol)
    img = nothing
    mimg_ff = nothing
    mimg_nf = nothing
    mimg_df = nothing
    GC.gc()
end


function testft_cimg(m, atol=1e-4)
    dx, dy = pixelsizes(m.img)
    u = fftshift(fftfreq(500, 1/dx))
    v = fftshift(fftfreq(500, 1/dy))
    mimg_ff = modelimage(m, FFTAlg(padfac=8))
    mimg_nf = modelimage(m, NFFTAlg(u, v))
    mimg_df = modelimage(m, DFTAlg(u, v))

    p = (U=u, V=v)
    vff = visibilities(mimg_ff, p)
    vnf = visibilities(mimg_nf, p)
    vdf = visibilities(mimg_df, p)

    @test isapprox(maximum(abs, vdf .- vnf), 0, atol=atol)
    @test isapprox(maximum(abs, vff .- vdf), 0, atol=atol)
    img = nothing
    mimg_ff = nothing
    mimg_nf = nothing
    mimg_df = nothing
    GC.gc()
end

@testset "CLEAN" begin
    _,vis, amp, lcamp, cphase = load_data()

    dirty_image(μas2rad(100.0), 128, vis)
    dirty_beam(μas2rad(100.0), 128, vis)
end

@testset "PolarizedModel Tests" begin
    u = randn(100)*0.5
    v = randn(100)*0.5
    t = sort(rand(100)*0.5)
    f = fill(230e9, 100)


    mI = stretched(MRing((0.2,), (0.1,)), 20.0, 20.0)
    mQ = 0.2*stretched(MRing((0.0,), (0.6,)), 20.0, 20.0)
    mU = 0.2*stretched(MRing((0.1,), (-0.6,)), 20.0, 20.0)
    mV = 0.0*stretched(MRing((0.0,), (-0.6,)), 20.0, 20.0)
    m = PolarizedModel(mI, mQ, mU, mV)
    @inferred visibility(m, (U=0.0, V=0.0))
    @inferred Comrade.intensity_point(m, (X=0.0, Y=0.0))


    function foo(x)
        m = PolarizedModel(
            stretched(Gaussian(), x[1], x[2]),
            stretched(Gaussian(), x[3], x[4]),
            shifted(Gaussian(), x[5], x[6]),
            x[7]*Gaussian()
        )
        vis = Comrade._coherency(Comrade.visibilities_analytic(m, u, v, t, f), CirBasis)
        Σ = map(x->real.(x .+ 1), zero.(vis))
        l = Comrade.CoherencyLikelihood(vis, Σ, 0.0)
        return logdensityof(l, zero.(vis))
    end

    x = rand(7)
    foo(x)
    testgrad(foo, x)
end


@testset "methods " begin
    _,_, amp, lcamp, cphase = load_data()

    m = rotated(stretched(Gaussian(), μas2rad(2.0), μas2rad(1.0)), π/8)
    u1, v1, u2, v2, u3, v3 = uvpositions(cphase[1])
    @test closure_phase(m, (U=u1, V=v1), (U=u2, V=v2), (U=u3, V=v3)) ≈ 0.0

    u1, v1, u2, v2, u3, v3, u4, v4 = uvpositions(lcamp[1])
    logclosure_amplitude(m, (U=u1, V=v1), (U=u2, V=v2), (U=u3, V=v3), (U=u4, V=v4))

end

@testset "Closure Cache" begin
    img = intensitymap(rotated(stretched(Gaussian(), μas2rad(2.0), μas2rad(1.0)), π/8),
                       μas2rad(12.0), μas2rad(12.0), 24, 12)
    _,_, amp, lcamp, cphase = load_data()

    cimg = ContinuousImage(img, DeltaPulse())
    cache_nf = create_cache(NFFTAlg(amp), img, DeltaPulse())
    cimg2 = ContinuousImage(img, cache_nf)
    cache_df = create_cache(DFTAlg(amp), img, DeltaPulse())
    cimg3 = ContinuousImage(img, cache_df)

    ac_amp = arrayconfig(amp)
    ac_lcamp = arrayconfig(lcamp)
    ac_cphase = arrayconfig(cphase)

    mimg_nf = modelimage(cimg, cache_nf)
    mimg_df = modelimage(cimg, cache_df)

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
end
