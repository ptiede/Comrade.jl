using Comrade, Optimization
using Pyehtim, OptimizationOptimJL, Distributions, VLBIImagePriors
using Zygote
using Test

include(joinpath(@__DIR__, "../test/test_util.jl"))

@testset "ComradeOptimizationExt.jl" begin
    _, _, _, lcamp, cphase = load_data()
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256)
    skym = SkyModel(test_model, test_prior(), g)
    post = VLBIPosterior(skym, lcamp, cphase)

    x0 = transform(tpost, [
            0.0,
            -0.4,
            0.0,
            2.0,
            0.0,
            -1.0,
            0.0,
            0.5,
            2.0,
            2.0,
        ])

    xopt, sol = comrade_opt(post, LBFGS(), AutoZygote(); initial_params=x0, maxiters=10_000)

    @test isapprox(xopt.f1/xopt.f2, 2.0, atol=1e-3)
    @test isapprox(xopt.σ1*2*sqrt(2*log(2)), μas2rad(40.0), rtol=1e-3)
    @test isapprox(xopt.σ1*xopt.τ1*2*sqrt(2*log(2)), μas2rad(20.0), rtol=1e-3)
    @test isapprox(xopt.ξ1, π/3, atol=1e-3)
    @test isapprox(xopt.σ2*2*sqrt(2*log(2)), μas2rad(20.0), atol=1e-3)
    @test isapprox(xopt.σ2*xopt.τ2*2*sqrt(2*log(2)), μas2rad(10.0), rtol=1e-3)
    @test isapprox(xopt.ξ2, π/6, atol=1e-3)
    @test isapprox(xopt.x, μas2rad(30.0), rtol=1e-3)
    @test isapprox(xopt.y, μas2rad(30.0), rtol=1e-3)
    @test chi2(post, xopt) < 0.1

end
