@testset "image modifiers" begin

    m = Gaussian()
    img = intensitymap(m, 20.0, 20.0, 128, 128)

    @testset "Rotate invariant" begin
        img2 = rotated(img, pi/4)
        @test isapprox(img2, img, rtol=1e-4)
    end

    @testset "Stretched" begin
        m2 = stretched(m, 2.0, 1.0)
        imgs  = intensitymap(m2, 20.0, 20.0, 128, 128)
        imgs2 = stretched(img, 2.0, 1.0)
        @test isapprox(imgs2, imgs, rtol=1e-4)
    end

    @testset "Stretch and rotate" begin
        m2 = modify(m, Stretch(2.0, 1.0), Rotate(π/4))
        imgs = intensitymap(m2, axiskeys(img))

        imgs2 = modify(img, Stretch(2.0, 1.0), Rotate(π/4))

        @test isapprox(imgs2, imgs, rtol=1e-4)
    end

    @testset "convolve" begin
        cimg = Comrade.convolve(img, Gaussian())
        img2 = modify(img, Stretch(√(2.0)))
        @test isapprox(cimg, img2, rtol=1e-4)
    end
end

@testset "CLEAN image" begin
    _,vis, _, _, _ = load_data()

    dirty_image(μas2rad(100.0), 128, vis)
    dirty_beam(μas2rad(100.0), 128, vis)

end
