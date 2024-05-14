using Pyehtim, Comrade, ComradeAdaptMCMC, Distributions, VLBIImagePriors
using Test

include(joinpath(@__DIR__, "../../../test/test_util.jl"))

@testset "ComradeAdaptMCMC.jl" begin

    m, vis, amp, lcamp, cphase = load_data()
    prior = test_prior()
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 256, 256)
    skym = SkyModel(test_model, test_prior(), g)
    post = VLBIPosterior(skym, lcamp, cphase)
    a1 = AdaptMCMC(ntemp=5)
    a2 = AdaptMCMC(ntemp=5, all_levels=true)

    x0 = prior_sample(post)

    chain = sample(post, a1, 10_000; thin=5)
    chain = sample(post, a2, 10_000; thin=5)
    chain = sample(post, a1, 10_000; thin=5, initial_params=x0)

    cpost = ascube(post)
    l0 = logdensityof(cpost, Comrade.HypercubeTransform.inverse(cpost, x0))
    @test l0 < logdensityof(cpost, Comrade.HypercubeTransform.inverse(cpost, chain[end]))
end
