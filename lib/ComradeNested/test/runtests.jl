using Comrade, ComradeNested, Distributions
using Test

load_ehtim()
include(joinpath(@__DIR__, "../../../test/test_util.jl"))

@testset "ComradeNested.jl" begin
    m, vis, amp, lcamp, cphase = load_data()
    prior = test_prior()
    lklhd = RadioLikelihood(test_model, lcamp, cphase)
    post = Posterior(lklhd, prior)
    a1 = Nested(dimension(post), 1000)

    chain, stats = sample(post, a1; dlogz=0.01, progress=false)


    #cpost = ascube(post)
    #xopt = chain[end]
    #@test isapprox(xopt.f1/xopt.f2, 2.0, atol=1e-2)
    #@test isapprox(xopt.σ1*2*sqrt(2*log(2)), μas2rad(40.0), rtol=1e-3)
    #@test isapprox(xopt.σ1*xopt.τ1*2*sqrt(2*log(2)), μas2rad(20.0), rtol=1e-3)
    #@test isapprox(-xopt.ξ1, π/3, atol=1e-2)
    #@test isapprox(xopt.σ2*2*sqrt(2*log(2)), μas2rad(20.0), atol=1e-2)
    #@test isapprox(xopt.σ2*xopt.τ2*2*sqrt(2*log(2)), μas2rad(10.0), rtol=1e-2)
    #@test isapprox(-xopt.ξ2, π/6, atol=1e-2)
    #@test isapprox(xopt.x, μas2rad(30.0), rtol=1e-2)
    #@test isapprox(xopt.y, μas2rad(30.0), rtol=1e-2)


end
