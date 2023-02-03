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

obs = load_ehtim_uvfits(joinpath(@__DIR__, "../assets/SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
# Now we will kill 0-baselines since we don't care about large scale flux and
# since we know that the gains in this dataset are coherent across a scan we make scan-average data
obs = scan_average(obs.flag_uvdist(uv_min=0.1e9))
# Now we extract the data products we want to fit:
#   1. log closure amplitudes
#   2. closure phases
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)

# Now the next thing we need to specify is our image model.
# !!! note
#     If we were fitting other data products like complex visibilities we would also need
#     to include a intrument or RIME model to forward model the gains.

# For the image model we will be using a modified `MRing`, which is a
# infinitely thin delta ring with an azimuthal structure given by a Fourier expansion.
# To give the ring some width we will convolve the ring with a gaussian, and add an
# additional gaussian to the image to model any non-ring flux. For the model a user
# must give a function that accepts a named tuple and return the constructed model:
# !!! note
#    The function model must always return an object that implements the Comrade [`Model Interface`](@ref)

function model(θ)
    (;radius, width, α, β, f, σG, τG, ξG, xG, yG) = θ
    ring = f*smoothed(stretched(MRing((α,), (β,)), radius, radius), width)
    g = (1-f)*shifted(rotated(stretched(Gaussian(), σG, σG*(1+τG)), ξG), xG, yG)
    return ring + g
end

# Now we can construct our likelihood `p(V|M)` where `V` is our data and `M` is our model.
# The first argument of `RadioLikelihood` is always a function that contructs our Comrade
# model from the set of parameters `θ`
lklhd = RadioLikelihood(model, dlcamp, dcphase)

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

post = Posterior(lklhd, prior)

# This constructs a posterior density that can be evaluated by calling `logdensityof`.
# For example,

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

# ## Reconstruction
# Now that we have fully specified our model we now will try to find the optimal reconstruction
# of our model given our observed data.

# Currently `post` is in **parameter** space. Often optimization and sampling algorithms
# want it in some modified space. For example, nested sampling algorithms want the
# parameters in the unit hypercube. To transform the posterior to the unit hypercube we
# can use the `ascube` function

cpost = ascube(post)

# If we want to flatten the parameter space and move from constrained parameters to (-∞, ∞)
# support we can use the `asflat` function

fpost = asflat(post)

# These transformed posterior expect a vector of parameters. That is we can evaluate the
# transformed log density by calling

logdensityof(cpost, rand(dimension(cpost)))
logdensityof(fpost, randn(dimension(fpost)))

# note that `cpost` logdensity vector expects that each element lives in `[0,1]`.


# ### Finding the Optimal Image

# Typically most VLBI modeling codes only care about finding the optimal or best guess
# image of our posterior `post` To do this we will use [`Optimization.jl`](@ref) and
# specifically the [`BlackBoxOptim.jl`](@ref) package. For Comrade this workflow is
# very similar to the usual `Optimization.jl` workflow. The only thing to keep in
# mind is that `Optimization.jl` expects that the function we are evaluating expects the
# parameters to be represented as a flat `Vector` of float. Therefore, we must use
# one of our transformed posteriors `cpost` or `fpost`. For the purposes of this example
# we will use `cpost` since it restricts the domain to live within the compact Unit hypercube
# which is easier to explore for non-gradient based optimizers like `BBO`.

using ComradeOptimization
using OptimizationBBO

ndim = dimension(fpost)
f = OptimizationFunction(fpost)
prob = OptimizationProblem(f, randn(ndim), nothing, lb=fill(-5.0, ndim), ub=fill(5.0, ndim))

# Now we solve for our optimial image.

sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=50_000)

# The sol vector is in the transformed space, so first we need to transform back to parameter space
# to that we can give meaning back to our solution.

xopt = transform(fpost, sol)

# Given this we can now plot the optimal image or the *maximum a posteriori* (MAP) image.

using Plots
plot(model(xopt), title="MAP image", xlims=(-60.0,50.0), ylims=(-60.0,50.0))

# ### Quantifying the Uncertainty of the Reconstruction

# While finding the optimal image is useful often in science the most important thing is to
# quantify the certainty of our inferences. This is really the goal of Comrade. In the language
# of Bayesian statistics we want to find representation of the posterior of possible image
# reconstructions given our choice of model and the data.
#
# Comrade provides a number of sampling and other posterior approximation tools. To see the
# list please see [`Libraries`](@ref). For this example we will be using
# [AdvancedHMC.jl](https://github.com/TuringLang/AdvancedHMC.jl) which uses
# an adaptive Hamiltonian Monte Carlo sampler called NUTS to approximate the posterior.
# Most of `Comrade`'s external libraries follow a very similar interface. To use AdvancedHMC
# do the following:

using ComradeAHMC
chain, stats = sample(post, AHMC(metric=DiagEuclideanMetric(ndim)), 2000; nadapts=1000, init_params=xopt)

# That's it! To finish it up we can then plot some simple visual fit diagnostics.

# First to plot the image we call

plot(model(chain[end]), title="Random image", xlims=(-60.0,50.0), ylims=(-60.0,50.0))

# What about the mean image? Well let's grab 100 images from the chain, where we first remove the
# adaptation steps since they don't sample from the correct posterior distribution
meanimg = mean(intensitymap.(model.(sample(chain[1000:end], 100)), μas2rad(120.0), μas2rad(120.0), 128, 128))
plot(sqrt.(max.(meanimg, 0.0)), title="Mean Image") #plot on a sqrt color scale to see the Gaussian

# That looks similar to the EHTC VI, and it took us no time at all!. To see how well the
# model is fitting the data we can plot the model and data products

plot(model(xopt), dlcamp, label="MAP")

# We can also plot what many draws from the posterior look like
p = plot(dlcamp)
uva = [sqrt.(uvarea(dlcamp[i])) for i in 1:length(dlcamp)]
for i in 1:10
    m = logclosure_amplitudes(model(chain[rand(1000:2000)]), arrayconfig(dlcamp))
    scatter!(uva, m, color=:grey, label=:none, alpha=0.1)
end
p

# Finally, we can also put everything onto a common scale and plot the normalized residuals.
# The normalied residuals are the difference between the data
# and the model, divided by the data's error:

residual(model(xopt), dlcamp)

# All of the diagnostic plots suggest that the model is missing some emission.
# In fact, this model is slightly too simple to explain the data.
# Check out [EHTC VI 2019](https://iopscience.iop.org/article/10.3847/2041-8213/ab1141)
# for some ideas what features need to be added to the model to get a better fit!


# For a real run we should also check that the MCMC chain has converged. For
# this we can use MCMCDiagnosticTools
using MCMCDiagnosticTools, Tables
# First lets look at the effective sample size (ESS) and R̂. This is important since
# the Monte Carlo standard error for MCMC estimates is proportional to 1/√ESS (for some problems)
# and R̂ is a measure of chain convergence. To find both we can use:
essrhat = map(ess_rhat∘(x->reshape(x, :, 1, 1)), Tables.columns(chain))
# Here the first value is the ESS and the second is the R̂. Note that we typically want R̂ < 1.01
# for all parameters, but you should also be running the problem at least 4 times from 4 different
# starting locations.

# In our example here we see that we have an ESS > 100 for all parameters and the R̂ < 1.01
# meaning that our MCMC chain mixing well. For more diagnostics see
# [`MCMCDiagnosticTools.jl`](https://turinglang.github.io/MCMCDiagnosticTools.jl/stable/).


# Computing information
# ```
# Julia Version 1.7.3
# Commit 742b9abb4d (2022-05-06 12:58 UTC)
# Platform Info:
#   OS: Linux (x86_64-pc-linux-gnu)
#   CPU: 11th Gen Intel(R) Core(TM) i7-1185G7 @ 3.00GHz
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-12.0.1 (ORCJIT, tigerlake)
# ```
