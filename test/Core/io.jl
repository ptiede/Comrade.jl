using FileIO
@testset "io.jl" begin
    imc, head = load(joinpath(@__DIR__, "../example_image.fits"), IntensityMap)
    ime = ehtim.image.load_image(joinpath(@__DIR__, "../example_image.fits"))
    data = load_data()
    @test size(ime.imarr("I")) == size(imc)
    @test flux(imc) ≈ ime.total_flux()
    @test imc.fov[1] ≈ ime.fovx()
    @test imc.fov[2] ≈ ime.fovx()
    save("test.fits", imc)
    rm("test.fits")
    save("test.fits", imc, last(data))
    rm("test.fits")
end
