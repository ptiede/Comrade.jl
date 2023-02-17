using Comrade, ComradeDynesty, Distributions
using Test

load_ehtim()
include(joinpath(@__DIR__, "../../../test/test_util.jl"))

@testset "ComradeDynesty.jl" begin
    m, vis, amp, lcamp, cphase = load_data()
    prior = test_prior()
    lklhd = RadioLikelihood(test_model, lcamp, cphase)
    post = Posterior(lklhd, prior)
    a1 = NestedSampler(dimension(post))
    a2 = DynamicNestedSampler(dimension(post))

    samples, dy = sample(post, a1; dlogz=0.1, print_progress=false)
    samples, dy = sample(post, a2; print_progress=false)

    esamples = ComradeDynesty.equalresample(samples, 1000)

    cpost = ascube(post)
    xopt = samples.posterior[draw=end, chain=1]
    @test isapprox(xopt.f1/xopt.f2, 2.0, atol=1e-2)
    @test isapprox(xopt.σ1*2*sqrt(2*log(2)), μas2rad(40.0), rtol=1e-2)
    @test isapprox(xopt.σ1*xopt.τ1*2*sqrt(2*log(2)), μas2rad(20.0), rtol=1e-2)
    @test isapprox(xopt.ξ1, π/3, atol=1e-2)
    @test isapprox(xopt.σ2*2*sqrt(2*log(2)), μas2rad(20.0), atol=1e-2)
    @test isapprox(xopt.σ2*xopt.τ2*2*sqrt(2*log(2)), μas2rad(10.0), rtol=1e-2)
    @test isapprox(xopt.ξ2, π/6, atol=1e-2)
    @test isapprox(xopt.x, μas2rad(30.0), rtol=1e-2)
    @test isapprox(xopt.y, μas2rad(30.0), rtol=1e-2)


end
