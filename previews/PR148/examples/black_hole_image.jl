# # Making an Image of a Black Hole

# `Comrade` has been designed to work with the EHT and ngEHT.
# In this tutorial we will show how to reproduce some of the results
# from [EHTC VI 2019](https://iopscience.iop.org/article/10.3847/2041-8213/ab1141).

# In EHTC VI, they considered fitting simple geometric models to the data
# to estimate the image size, shape, brightness profile etc of the black hole.
# In this page we will construct a similar model and fit it to the data in under
# 50 lines of code (sans comments). To start we load some packages we will need

using Comrade

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
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)

# For this demo only consider closure products since these are invariant to station specific
# gain systematics. Given these data products we can then form our radio likelihood:

lklhd = RadioLikelihood(dlcamp, dcphase)

# The `lklhd` constructs a measure using `MeasureTheory.jl`. To evaluate the likelihood
# we need to pass it a type that implements the `Comrade` model interface which is described
# in [Model Interface](@ref).

# To finish the construction of our posterior we need to specify an image model and a prior.
# For the image model we will be using a modified `MRing`, which is a
# infinitely thin delta ring with an azimuthal structure given by a Fourier expansion.
# To give the ring some width we will convolve the ring with a gaussian, and add an
# additional gaussian to the image to model any non-ring flux. For the model a user
# must give a function that accepts a named tuple and return the constructed model:

function model(θ)
  (;radius, width, α, β, f, σG, τG, ξG, xG, yG) = θ
  ring = f*smoothed(stretched(MRing((α,), (β,)), radius, radius), width)
  g = (1-f)*shifted(rotated(stretched(Gaussian(), σG, σG*(1+τG)), ξG), xG, yG)
  return ring + g
end

# We now need to specify the priors for our model. The easiest way to do this is to
# specify a NamedTuple of distributions:

using Distributions
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

# This constructs a posterior density that can be evaluated by calling `logdensityof`.
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

# Now this model is in **parameter** space. Often optimization and sampling algorithms
# want it in some modified space. For example, nested sampling algorithms want the
# parameters in the unit hypercube. To transform the posterior to the unit hypercube we
# can use the `ascube` function

cpost = ascube(post)

# If we want to flatten the parameter space and move to (-∞, ∞) support we can use the
# `asflat` function

fpost = asflat(post)

# These transformed posterior expect a vector of parameters. That is we can evaluate the
# transformed log density by calling

logdensityof(cpost, rand(dimension(cpost)))
logdensityof(fpost, randn(dimension(fpost)))

# Note that this automatically takes care of the jacobian in the parameter transformation.

# ## Sampling the posterior

# Our strategy here will be to use Hamiltonian Monte Carlo. However, to lower burn-in time
# we will first use an optimizer to find a reasonable starting location. Since this is a lower
# dimensional problem we will use BlackboxOptim or the GalacticBBO package

using ComradeGalactic
using GalacticBBO

ndim = dimension(fpost)
f = OptimizationFunction(fpost)
prob = OptimizationProblem(f, randn(ndim), nothing, lb=fill(-5.0, ndim), ub=fill(5.0, ndim))
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=50_000)

# The sol vector is in the transformed space, so first we need to transform back to parameter space

xopt = transform(fpost, sol)

# And we can also plot the map

using Plots
plot(model(xopt), title="MAP image", xlims=(-60.0,50.0), ylims=(-60.0,50.0))

# The main goal of `Comrade` is to explore the posterior of the model parameters.
# Currently the go to tool is [AdvancedHMC.jl](https://github.com/TuringLang/AdvancedHMC.jl).
# To sample from the posterior you can use the following:

using ComradeAHMC
chain, stats = sample(post, AHMC(metric=DiagEuclideanMetric(ndim)), 2000; nadapts=1000, init_params=xopt)

# That's it! To finish it up we can then plot some simple visual fit diagnostics.

# First to plot the image we call

plot(model(chain[end]), title="Random image", xlims=(-60.0,50.0), ylims=(-60.0,50.0))

# What about the mean image? Well let's grab 100 images from the chain
meanimg = mean(intensitymap.(model.(sample(chain, 100)), μas2rad(120.0), μas2rad(120.0), 128, 128))
plot(sqrt.(meanimg), title="Mean Image") #plot on a sqrt color scale to see the Gaussian

# That looks similar to the EHTC VI, and it took us no time at all!. To see how well the
# model is fitting the data we can plot the model and data products

plot(model(xopt), dlcamp)

# or even better the residuals

residual(model(xopt), dlcamp)

# These residuals are a little ratty which suggests our model is missing some emission.
# In fact, this model is slightly too simple to explain the data.
# Check out [EHTC VI 2019](https://iopscience.iop.org/article/10.3847/2041-8213/ab1141)
# for some ideas what features need to be added to the model to get a better fit!
