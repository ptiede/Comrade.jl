# # Making an Image of a Black Hole

# `Comrade` has been designed to work with the EHT and ngEHT.
# In this tutorial, we will show how to reproduce some of the results
# from [EHTC VI 2019](https://iopscience.iop.org/article/10.3847/2041-8213/ab1141).

# In EHTC VI, they considered fitting simple geometric models to the data
# to estimate the black hole's image size, shape, brightness profile, etc.
# In this tutorial, we will construct a similar model and fit it to the data in under
# 50 lines of code (sans comments). To start, we load Comrade and some other packages we need.

using Comrade

using Pkg #hide
Pkg.activate(joinpath(dirname(pathof(Comrade)), "..", "examples")) #hide

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(42)
#-

# The next step is to load the data. We will use the publically
# available M 87 data which can be downloaded
# from [cyverse](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/EHTC_FirstM87Results_Apr2019).
# For an introduction to data loading, see [Loading Data into Comrade](@ref).

obs = load_ehtim_uvfits(joinpath(dirname(pathof(Comrade)), "..", "examples", "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
# Now we will kill 0-baselines since we don't care about large-scale flux and
# since we know that the gains in this dataset are coherent across a scan, we make scan-average data
obs = scan_average(obs.flag_uvdist(uv_min=0.1e9))
# Now we extract the data products we want to fit:
#   1. log closure amplitudes
#   2. closure phases
dvis = extract_vis(obs)

# For the image model, we will use a modified `MRing`, a
# infinitely thin delta ring with an azimuthal structure given by a Fourier expansion.
# To give the MRing some width, we will convolve the ring with a Gaussian and add an
# additional gaussian to the image to model any non-ring flux.
# Comrade expects that any model function must accept a named tuple and returns  must always return an object that implements the Comrade [Model Interface](@ref)
#-
function model(θ)
    (;σG, τG, ξG, xG, yG) = θ
    g = shifted(rotated(stretched(Gaussian(), σG, σG*(1+τG)), ξG), xG, yG)
    return g
end

# To construct our likelihood `p(V|M)` where `V` is our data and `M` is our model, we use the `RadioLikelihood` function.
# The first argument of `RadioLikelihood` is always a function that constructs our Comrade
# model from the set of parameters `θ`.
lklhd = RadioLikelihood(model, dvis)

# We now need to specify the priors for our model. The easiest way to do this is to
# specify a NamedTuple of distributions:

using Distributions
prior = (
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

# If we want to flatten the parameter space and move from constrained parameters to (-∞, ∞)
# support we can use the `asflat` function

fpost = asflat(post)

struct TransformPrior{T,P}
    t::T
    p::P
end

import TransformVariables as TV
function (θ::TransformPrior)(x)
    y, ly = TV.transform_and_logjac(θ.t, x)
    return logdensityof(θ.p, y) + ly
end

# These transformed posterior expect a vector of parameters. That is we can evaluate the
# transformed log density by calling

# logdensityof(cpost, rand(rng, dimension(cpost)))
logdensityof(fpost, randn(rng, dimension(fpost)))

# note that `cpost` logdensity vector expects that each element lives in `[0,1]`.


# ### Finding the Optimal Image

# Typically, most VLBI modeling codes only care about finding the optimal or best guess
# image of our posterior `post` To do this, we will use [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/) and
# specifically the [`BlackBoxOptim.jl`](https://github.com/robertfeldt/BlackBoxOptim.jl) package. For Comrade, this workflow is
# very similar to the usual `Optimization.jl` workflow. The only thing to keep in
# mind is that `Optimization.jl` expects that the function we are evaluating expects the
# parameters to be represented as a flat `Vector` of float. Therefore, we must use
# one of our transformed posteriors, `cpost` or `fpost`. For this example
#, we will use `cpost` since it restricts the domain to live within the compact unit hypercube
#, which is easier to explore for non-gradient-based optimizers like `BBO`.

using ComradeOptimization
using OptimizationBBO

ndim = dimension(fpost)
f = OptimizationFunction(fpost)
prob = Optimization.OptimizationProblem(f, randn(rng, ndim), nothing, lb=fill(-5.0, ndim), ub=fill(5.0, ndim))

# Now we solve for our optimial image.

sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=50_000)

# The sol vector is in the transformed space, so first we need to transform back to parameter space
# to that we can interpret the solution.

xopt = transform(fpost, sol)

# Given this we can now plot the optimal image or the *maximum a posteriori* (MAP) image.

using Plots
plot(model(xopt), title="MAP image", xlims=(-60.0,50.0), ylims=(-60.0,50.0))

# ### Quantifying the Uncertainty of the Reconstruction

# While finding the optimal image is often helpful, in science, the most important thing is to
# quantify the certainty of our inferences. This is the goal of Comrade. In the language
# of Bayesian statistics, we want to find a representation of the posterior of possible image
# reconstructions given our choice of model and the data.
#
# Comrade provides several sampling and other posterior approximation tools. To see the
# list, please see the Libraries section of the docs. For this example, we will be using
# [AdvancedHMC.jl](https://github.com/TuringLang/AdvancedHMC.jl), which uses
# an adaptive Hamiltonian Monte Carlo sampler called NUTS to approximate the posterior.
# Most of Comrade's external libraries follow a similar interface. To use AdvancedHMC
# do the following:

using ComradeAHMC
chain, stats = sample(rng, post, AHMC(metric=DiagEuclideanMetric(ndim)), 2000; nadapts=1000, init_params=xopt)

# That's it! To finish it up we can then plot some simple visual fit diagnostics.

# First to plot the image we call

plot(model(chain[end]), title="Random image", xlims=(-60.0,60.0), ylims=(-60.0,60.0))

# What about the mean image? Well let's grab 100 images from the chain, where we first remove the
# adaptation steps since they don't sample from the correct posterior distribution
meanimg = mean(intensitymap.(model.(sample(chain[1000:end], 100)), μas2rad(120.0), μas2rad(120.0), 128, 128))
plot(sqrt.(max.(meanimg, 0.0)), title="Mean Image") #plot on a sqrt color scale to see the Gaussian

# That looks similar to the EHTC VI, and it took us no time at all!. To see how well the
# model is fitting the data we can plot the model and data products

plot(model(xopt), dlcamp, label="MAP")

# We can also plot what many draws from the posterior look like
p = plot(dlcamp);
uva = [sqrt.(uvarea(dlcamp[i])) for i in 1:length(dlcamp)]
for i in 1:10
    m = logclosure_amplitudes(model(chain[rand(rng, 1000:2000)]), arrayconfig(dlcamp))
    scatter!(uva, m, color=:grey, label=:none, alpha=0.1)
end
p

# Finally, we can also put everything onto a common scale and plot the normalized residuals.
# The normalied residuals are the difference between the data
# and the model, divided by the data's error:

residual(model(xopt), dlcamp)

# All diagnostic plots suggest that the model is missing some emission sources.
# In fact, this model is too simple to explain the data.
# Check out [EHTC VI 2019](https://iopscience.iop.org/article/10.3847/2041-8213/ab1141)
# for some ideas about what features need to be added to the model to get a better fit!
