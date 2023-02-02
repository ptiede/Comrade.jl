function testpol(m)
    g = GriddedKeys(imagepixels(5.0, 5.0, 128, 128))
    img = intensitymap(m, g)
    img2 = zero(img)
    intensitymap!(img2, m)

    @test all(==(1), img .≈ img2)

    u = fftshift(fftfreq(length(g.X), 1/step(g.X)))
    uv = (U=u, V=u)
    v  = visibilities(m, uv)

    plot(m)
    plot(img)
    Plots.closeall()
    return v
end

@testset "Polarized Analytic" begin
    m = PolarizedModel(Gaussian(), 0.1*Gaussian(), 0.1*Gaussian(), 0.1*Gaussian())
    testpol(m)
end

@testset "Polarized Semi Analytic" begin
    m = PolarizedModel(ExtendedRing(8.0), 0.1*Gaussian(), 0.1*Gaussian(), 0.1*Gaussian())
    g = GriddedKeys(imagepixels(10.0, 10.0, 512, 512))
    s = map(length, dims(g))
    vff = testpol(modelimage(m, IntensityMap(zeros(StokesParams{Float64}, s), g), FFTAlg()))
    vnf = testpol(modelimage(m, IntensityMap(zeros(StokesParams{Float64}, s), g), NFFTAlg()))
    vdf = testpol(modelimage(m, IntensityMap(zeros(StokesParams{Float64}, s), g), DFTAlg()))

    @test isapprox(vff, vnf, atol=1e-6)
    @test isapprox(vff, vdf, atol=1e-6)
    @test isapprox(vnf, vdf, atol=1e-6)
end

@testset "Polarized Modified" begin
    g = GriddedKeys(imagepixels(5.0, 5.0, 128, 128))
    s = map(length, dims(g))
    m0 = PolarizedModel(ExtendedRing(2.0), 0.1*Gaussian(), 0.1*Gaussian(), 0.1*Gaussian())
    m = shifted(m0, 0.1 ,0.1)
    testpol(modelimage(m, IntensityMap(zeros(StokesParams{Float64}, s), g)))

    m = rotated(m0, 0.1)
    testpol(modelimage(m, IntensityMap(zeros(StokesParams{Float64}, s), g)))

    m = renormed(m0, 0.1)
    testpol(modelimage(m, IntensityMap(zeros(StokesParams{Float64}, s), g)))

    m = stretched(m0, 0.1, 0.4)
    testpol(modelimage(m, IntensityMap(zeros(StokesParams{Float64}, s), g)))

end

@testset "Polarized Combinators" begin
    m1 = PolarizedModel(Gaussian(), 0.1*Gaussian(), 0.1*Gaussian(), 0.1*Gaussian())
    m2 = PolarizedModel(Disk(), shifted(Disk(), 0.1, 1.0), ZeroModel(), ZeroModel())
    testpol(convolved(m1,m2))
end

@testset "Polarized All Mod" begin
    m1 = PolarizedModel(Gaussian(), 0.1*Gaussian(), 0.1*Gaussian(), 0.1*Gaussian())
    m2 = PolarizedModel(ExtendedRing(8.0), shifted(Disk(), 0.1, 1.0), ZeroModel(), ZeroModel())
    m = convolved(m1,m2)+m1
    g = GriddedKeys(imagepixels(5.0, 5.0, 128, 128))
    s = map(length, dims(g))
    testpol(modelimage(m, IntensityMap(zeros(StokesParams{Float64}, s), g)))
end

@testset "Rotation" begin
    m = PolarizedModel(Gaussian(), Gaussian(), ZeroModel(), 0.1*Gaussian())
    g = GriddedKeys(imagepixels(5.0, 5.0, 128, 128))
    img1 = intensitymap(m, g)
    @test all(==(1), stokes(img1, :Q) .≈ stokes(img1, :I))
    @test all(==(1), stokes(img1, :U) .≈ 0.0*stokes(img1, :I))
    @test all(==(1), stokes(img1, :V) .≈ 0.1*stokes(img1, :I))

    # Now we rotate into perfect U
    img2 = intensitymap(rotated(m,  π/4), g)
    @test all(==(1), stokes(img2, :U) .≈ stokes(img1, :I))
    @test all(==(1), isapprox.(stokes(img2, :Q), 0.0, atol=1e-16))
    @test all(==(1), stokes(img1, :V) .≈ stokes(img2, :V))

    # Now we rotate into perfectly -Q
    img3 = intensitymap(rotated(m,  π/2), g)
    @test all(==(1), stokes(img3, :Q) .≈ -stokes(img3, :I))
    @test all(==(1), isapprox.(stokes(img3, :U), 0.0, atol=1e-16))
    @test all(==(1), stokes(img1, :V) .≈ stokes(img3, :V))

    # Now we rotate into perfectly -U
    img4 = intensitymap(rotated(m,  3π/4), g)
    @test all(==(1), stokes(img4, :U) .≈ -stokes(img4, :I))
    @test all(==(1), isapprox.(stokes(img4, :Q), 0.0, atol=1e-16))
    @test all(==(1), stokes(img1, :V) .≈ stokes(img4, :V))

    # Now make sure it is π periodic
    img5 = intensitymap(rotated(m,  -π/4), g)
    @test all(==(1), img4 .≈ img5)
end
