using Comrade, ComradeAHMC, StatsBase
using Test

include(joinpath(@__DIR__, "../../../test/test_util.jl"))

@testset "ComradeAHMC.jl" begin

    _, _, lcamp, cphase = load_data()
    lklhd = RadioLikelihood(lcamp, cphase)

    prior = test_prior()
    post = Posterior(lklhd, prior, test_model)

    hchain, hstats, res = sample(post, AHMC(metric=DiagEuclideanMetric(dimension(post))),
                                       10_000; nadapts=2000, init_params=nchain[end], progress=true)
    mh = Comrade.rmap(mean, hchain)[1]
    sh = Comrade.rmap(std, hchain)[1]

    @test isapprox(collect(values(mh)), collect(values(mn)), rtol=1e-3)
    @test isapprox(collect(values(sh)), collect(values(sn)), rtol=1e-1)



end
