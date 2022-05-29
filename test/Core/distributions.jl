using MeasureBase
using Distributions

@testset "distributions" begin
    # d = Rice(ν=2.0, σ=0.5)
    # logdensityof(d, 1.0)

    # samples = [rand(d) for _ in 1:100_000]
    # ms = mean(samples)
    # vs = var(samples)
    # @test isapprox(ms, mean(d); atol=1e-2)
    # @test isapprox(vs, var(d); atol=1e-2)

    # dhsnr = Rice(ν=5.0, σ=1e-3)
    # μ = sqrt(dhsnr.ν^2 - dhsnr.σ^2)
    # @test isapprox(logdensityof(dhsnr, 4.9), logpdf(Normal(μ, dhsnr.σ), 4.9); atol=1e-1)

    dvm1 = CPVonMises(μ=0.0, σ=0.1)
    dvm2 = CPVonMises(μ=0.0, κ=1/0.1^2)

    # Take off the normalization piece
    ds = VonMises(0.0, 1/0.1^2)
    norm = -log(ds.I0κx) - log(2π)
    @test logpdf(ds, 0.5) - norm ≈ logdensityof(dvm2, 0.5)


    @test logdensityof(dvm1, 0.5) ≈ logdensityof(dvm2, 0.5)
end
