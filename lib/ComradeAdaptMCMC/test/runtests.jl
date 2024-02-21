using Pyehtim, Comrade, ComradeAdaptMCMC, Distributions, VLBIImagePriors
using Test

include(joinpath(@__DIR__, "../../../test/test_util.jl"))

@testset "ComradeAdaptMCMC.jl" begin

    m, vis, amp, lcamp, cphase = load_data()
    prior = test_prior()
    lklhd = RadioLikelihood(test_model, lcamp, cphase)
    post = Posterior(lklhd, prior)
    a1 = AdaptMCMC(ntemp=5)
    a2 = AdaptMCMC(ntemp=5, all_levels=true)

    x0 = prior_sample(post)

    chain = sample(post, a1, 100_000; thin=5)
    chain = sample(post, a2, 100_000; thin=5)
    chain = sample(post, a1, 100_000; thin=5, initial_params=x0)

    cpost = ascube(post)
    l0 = logdensityof(cpost, Comrade.HypercubeTransform.inverse(cpost, x0))
    @test l0 < logdensityof(cpost, Comrade.HypercubeTransform.inverse(cpost, chain[end]))
end
