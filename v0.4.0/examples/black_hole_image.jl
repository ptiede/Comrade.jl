# # Making an Image of a Black Hole

# `Comrade` has been designed to work with the EHT and ngEHT.
# In this tutorial we will show how to reproduce some of the results
# from [EHTC VI 2019](https://iopscience.iop.org/article/10.3847/2041-8213/ab1141).

# In EHTC VI, they considered fitting simple geometric models to the data
# to estimate the image size, shape, brightness profile etc of the black hole.
# In this page we will construct a similar model and fit it to the data in under
# 50 lines of code (sans comments). To start we load some packages we will need

using Comrade
using Plots
using Distributions
using Pathfinder
using AdvancedHMC

# Pathfinder and AdvancedHMC will be required for sampling the posterior.
# Distributions is a general purpose probability distribution package in Julia.
# The next step is to load the data. For this we will use the publically
# available M 87 data which can be downloaded
# from [cyverse](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/EHTC_FirstM87Results_Apr2019).
# For an introduction to data loading see [Loading Data into Comrade](@ref).

load_ehtim()
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "../assets/SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
obs.add_scans()
# kill 0-baselines since we don't care about large scale flux and make scan-average data
obs = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true)
# grab data products we want to fit: log closure amplitudes and closure phases
dlcamp = extract_lcamp(obs; count="min")
dcphase = extract_cphase(obs, count="min")

# For this demo only consider closure products since these are invariant to station specific
# gain systematics. Given these data products we can then form our visibility likelihood:

lklhd = RadioLikelihood(dlcamp, dcphase)

# The `lklhd` constructs a measure using `MeasureTheory.jl`. To evaluate the likelihood
# we need to pass it a type that implements the `Comrade` model interface which is described
# in [Model Interface](@ref).

# To finish the construction of our posterior we need to specify an image model and a prior.
# For the image model we will be using a modified `MRing`, which is a
# infinitely thin delta ring with an azimuthal structure given by a Fourier expansion.
# To give the ring some width we will convolve the ring with a gaussian, and add an
# additional gaussian to the image to model any non-ring flux. The resulting model function
# is given by

function model(θ)
  (;radius, width, α, β, f, σG, τG, ξG, xG, yG) = θ
  ring = f*smoothed(stretched(MRing((α,), (β,)), radius, radius), width)
  g = (1-f)*shifted(rotated(stretched(Gaussian(), σG, σG*(1+τG)), ξG), xG, yG)
  return ring + g
end

# We now need to specify the priors for our model. The easiest way to do this is to
# specify a NamedTuple of distributions.

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

post = Posterior(lklhd, prior, model)

# This constructs a posterior density that can be evaluated by calling `logdensity`.
# For example

logdensityof(post, (radius = μas2rad(20.0),
                  width = μas2rad(10.0),
                  α = 0.3,
                  β = 0.3,
                  f = 0.6,
                  σG = μas2rad(20.0),
                  τG = 0.1,
                  ξG = 0.5,
                  xG = 0.0,
                  yG = 0.0))

# We can now try to sample from our posterior `post` so that we can make probabilstic
# inferences about our data.

# Our strategy here will be to use Hamiltonian Monte Carlo. However, to lower burn-in time
# we will first use the excellent Pathfinder algorithm[^Zhang2021], which uses variational
# inference to find a starting point that is roughly drawn from the posterior. The benefit
# this is that we will greatly minimize burn-in time, and we get a variational approximation
# for free.

q, ϕ, _ = multipathfinder(post, 100)

# The q's are a uniformly weighted mixture of multivariate normal distributions, while
# the `ϕ`'s are our approximate posterior draws. We can then start an HMC run by calling

ndim = dimension(post)
chain, stats = sample(post, HMC(metric=DiagEuclideanMetric(ndim)), 2000; nadapts=1000, init_params=ϕ[1])

# That's it! We now have a converged chain. To finish it up we can then plot some diagnostics.

# First to plot the image we call

plot(model(chain[end]))

# That looks similar to the EHTC VI, and it took us no time at all!. To see how well the
# model is fitting the data we can plot the model and data products

plot(model(chain[end]), dlcamp)

# or even better the residuals

residual(model(chain[end]), dlcamp)
















# [^Zhang2021]: Lu Zhang, Bob Carpenter, Andrew Gelman, Aki Vehtari (2021).
# Pathfinder: Parallel quasi-Newton variational inference.
# arXiv: [2108.03782](https://arxiv.org/abs/2108.03782) [stat.ML].
# [Code](https://github.com/LuZhangstat/Pathfinder))
