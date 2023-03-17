using FileIO
@testset "io.jl" begin
    imc = Comrade.load(joinpath(@__DIR__, "../example_image.fits"), IntensityMap)
    ime = ehtim.image.load_image(joinpath(@__DIR__, "../example_image.fits"))
    data = load_data()
    @test size(ime.imarr("I")) == size(imc)
    @test flux(imc) ≈ ime.total_flux()
    fov = fieldofview(imc)
    @test fov.X ≈ pyconvert(Float64, ime.fovx())
    @test fov.Y ≈ pyconvert(Float64, ime.fovx())
    Comrade.save("test.fits", imc)
    rm("test.fits")
end
