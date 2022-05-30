using Comrade, ComradeAHMC
using Test

load_ehtim()
include(joinpath(@__DIR__, "../../../test/test_util.jl"))

@testset "ComradeAHMC.jl" begin

    _, _, _, lcamp, cphase = load_data()
    lklhd = RadioLikelihood(lcamp, cphase)

    prior = test_prior()
    post = Posterior(lklhd, prior, test_model)
    ndim = dimension(post)
    x0 = (f1 = 1.0916271439905998,
          σ1 = 8.230088139590025e-11,
          τ1 = 0.49840994275315254,
          ξ1 = -1.0489388962890198,
          f2 = 0.553944311447593,
          σ2 = 4.1283218512580634e-11,
          τ2 = 0.5076731020106428,
          ξ2 = -0.5376269092893298,
          x = 1.451956089157719e-10,
          y = 1.455983181049137e-10)
    s1 = AHMC(metric=DiagEuclideanMetric(ndim))
    s2 = AHMC(metric=DenseEuclideanMetric(ndim))
    s3 = AHMC(metric=DenseEuclideanMetric(ndim), autodiff=Comrade.AD.ZygoteBackend())
    hchain, hstats = sample(post, s1, 3_000; nadapts=2_000, progress=false)
    hchain, hstats = sample(post, s1, 3_000; nadapts=2_000, progress=false, init_params=x0)
    hchain, hstats = sample(post, s2, 3_000; nadapts=2_000, progress=false, init_params=x0)

    hchain, hstats = sample(post, s2, Comrade.AbstractMCMC.MCMCThreads(), 3_000, 2; nadapts=2_000, progress=false)
    hchain, hstats = sample(post, s2, Comrade.AbstractMCMC.MCMCThreads(), 3_000, 2; nadapts=2_000, progress=false,init_params=[x0,x0])

    cpost = asflat(post)
    l0 = logdensityof(cpost, Comrade.HypercubeTransform.inverse(cpost, x0))
    @test 10*l0  < mean(hstats[1].log_density[2000:end])

end
