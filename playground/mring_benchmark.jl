using Pkg; Pkg.activate(@__DIR__)
using Comrade
using VLBIImagePriors
using BenchmarkTools
using Distributions, DistributionsAD


# To download the data visit https://doi.org/10.25739/g85n-f134
obs = load_ehtim_uvfits(joinpath(@__DIR__, "../examples/", "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obs = scan_average(obs)
# extract log closure amplitudes and closure phases
dvis = extract_vis(obs)

function model(θ)
    (;radius, width, α, β, f, σG, τG, ξG, xG, yG) = θ
    ring = f*smoothed(stretched(MRing((α,), (β,)), radius, radius), width)
    g = (1-f)*shifted(rotated(stretched(Gaussian(), σG, σG*(1+τG)), ξG), xG, yG)
    return ring + g
end

# To construct our likelihood `p(V|M)` where `V` is our data and `M` is our model, we use the `RadioLikelihood` function.
# The first argument of `RadioLikelihood` is always a function that constructs our Comrade
# model from the set of parameters `θ`.
lklhd = RadioLikelihood(model, dvis)

# We now need to specify the priors for our model. The easiest way to do this is to
# specify a NamedTuple of distributions:

prior = (
          radius = Uniform(μas2rad(10.0), μas2rad(30.0)),
          width = Uniform(μas2rad(1.0), μas2rad(10.0)),
          α = Uniform(-0.5, 0.5),
          β = Uniform(-0.5, 0.5),
          f = Uniform(0.0, 1.0),
          σG = Uniform(μas2rad(1.0), μas2rad(40.0)),
          τG = Uniform(0.0, 0.75),
          ξG = Uniform(0.0, 1π),
          xG = Uniform(-μas2rad(80.0), μas2rad(80.0)),
          yG = Uniform(-μas2rad(80.0), μas2rad(80.0))
        )

# To form the posterior we now call

post = Posterior(lklhd, prior)

# This constructs a posterior density that can be evaluated by calling `logdensityof`.
# For example,

x0 = prior_sample(post)

@benchmark logdensityof($post, $x0)

using ForwardDiff
using AbstractDifferentiation
tpost = asflat(post)
ℓ = logdensityof(tpost)
x0 = prior_sample(tpost)
ab = AD.ForwardDiffBackend{10}()
@benchmark AD.gradient($ab, $ℓ, $x0)


using Zygote
@benchmark Zygote.gradient($ℓ, $x0)
