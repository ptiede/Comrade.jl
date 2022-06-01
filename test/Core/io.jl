@testset "io.jl" begin
    imc, head = load(joinpath(@__DIR__, "../example_image.fits"), IntensityMap)
    ime = ehtim.image.load_image(joinpath(@__DIR__, "../example_image.fits"))

    @test size(ime.imarr("I")) == size(imc)
    @test flux(imc) ≈ ime.total_flux()
    @test imc.fovx ≈ ime.fovx()
    save("test.fits", imc)
    rm("test.fits")
end
