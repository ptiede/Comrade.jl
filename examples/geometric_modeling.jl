# # Geometric Modeling of EHT Data

# `Comrade` has been designed to work with the EHT and ngEHT.
# In this tutorial, we will show how to reproduce some of the results
# from [EHTC VI 2019](https://iopscience.iop.org/article/10.3847/2041-8213/ab1141).

# In EHTC VI, they considered fitting simple geometric models to the data
# to estimate the black hole's image size, shape, brightness profile, etc.
# In this tutorial, we will construct a similar model and fit it to the data in under
# 50 lines of code (sans comments). To start, we load Comrade and some other packages we need.



using Comrade

# ## Load the Data
using Pkg #hide
Pkg.activate(joinpath(dirname(pathof(Comrade)), "..", "examples")) #hide

using Pyehtim

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(42)
#-
# The next step is to load the data. We will use the publically
# available M 87 data which can be downloaded
# from [cyverse](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/EHTC_FirstM87Results_Apr2019).
# For an introduction to data loading, see [Loading Data into Comrade](@ref).

obs = load_uvfits_and_array(joinpath(dirname(pathof(Comrade)), "..", "examples", "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
# Now we will kill 0-baselines since we don't care about large-scale flux and
# since we know that the gains in this dataset are coherent across a scan, we make scan-average data
obs = Pyehtim.scan_average(obs.flag_uvdist(uv_min=0.1e9))

# Now we extract the data products we want to fit
dlcamp, dcphase = extract_table(obs, LogClosureAmplitudes(;snrcut=3.0), ClosurePhases(;snrcut=3.0))

# !!!warn
#    We remove the low-snr closures since they are very non-gaussian. This can create rather
#    large biases in the model fitting since the likelihood has much heavier tails that the
#    usual Gaussian approximation.
#-

# For the image model, we will use a modified `MRing`, a
# infinitely thin delta ring with an azimuthal structure given by a Fourier expansion.
# To give the MRing some width, we will convolve the ring with a Gaussian and add an
# additional gaussian to the image to model any non-ring flux.
# Comrade expects that any model function must accept a named tuple and returns  must always
# return an object that implements the [VLBISkyModels Interface](https://ehtjulia.github.io/VLBISkyModels.jl/stable/interface/)
#-
function model(θ)
    (;radius, width, ma, mp, τ, ξτ, f, σG, τG, ξG, xG, yG) = θ
    α = ma.*cos.(mp .- ξτ)
    β = ma.*sin.(mp .- ξτ)
    ring = f*smoothed(modify(MRing(α, β), Stretch(radius, radius*(1+τ)), Rotate(ξτ)), width)
    g = (1-f)*shifted(rotated(stretched(Gaussian(), σG, σG*(1+τG)), ξG), xG, yG)
    return ring + g
end

# To construct our likelihood `p(V|M)` where `V` is our data and `M` is our model, we use the `RadioLikelihood` function.
# The first argument of `RadioLikelihood` is always a function that constructs our Comrade
# model from the set of parameters `θ`.
lklhd = RadioLikelihood(model, dlcamp, dcphase)

# We now need to specify the priors for our model. The easiest way to do this is to
# specify a NamedTuple of distributions:

using Distributions, VLBIImagePriors
prior = NamedDist(
          radius = Uniform(μas2rad(10.0), μas2rad(30.0)),
          width = Uniform(μas2rad(1.0), μas2rad(10.0)),
          ma = (Uniform(0.0, 0.5), Uniform(0.0, 0.5)),
          mp = (Uniform(0, 2π), Uniform(0, 2π)),
          τ = Uniform(0.0, 1.0),
          ξτ= Uniform(0.0, π),
          f = Uniform(0.0, 1.0),
          σG = Uniform(μas2rad(1.0), μas2rad(100.0)),
          τG = Uniform(0.0, 1.0),
          ξG = Uniform(0.0, 1π),
          xG = Uniform(-μas2rad(80.0), μas2rad(80.0)),
          yG = Uniform(-μas2rad(80.0), μas2rad(80.0))
        )


# Note that for `α` and `β` we use a product distribution to signify that we want to use a
# multivariate uniform for the mring components `α` and `β`. In general the structure of the
# variables is specified by the prior. Note that this structure must be compatible with the
# model definition `model(θ)`.

# To form the posterior we now call

post = Posterior(lklhd, prior)

# !!!warn
#    As of Comrade 0.9 we have switched to the proper covariant closure likelihood.
#    This is slower than the naieve diagonal liklelihood, but takes into account the
#    correlations between closures that share the same baselines.

# This constructs a posterior density that can be evaluated by calling `logdensityof`.
# For example,

logdensityof(post, (radius = μas2rad(20.0),
                  width = μas2rad(10.0),
                  ma = (0.3, 0.3),
                  mp = (π/2, π),
                  τ = 0.1,
                  ξτ= π/2,
                  f = 0.6,
                  σG = μas2rad(50.0),
                  τG = 0.1,
                  ξG = 0.5,
                  xG = 0.0,
                  yG = 0.0))

# ## Reconstruction

# Now that we have fully specified our model, we now will try to find the optimal reconstruction
# of our model given our observed data.

# Currently, `post` is in **parameter** space. Often optimization and sampling algorithms
# want it in some modified space. For example, nested sampling algorithms want the
# parameters in the unit hypercube. To transform the posterior to the unit hypercube, we
# can use the `ascube` function

cpost = ascube(post)

# If we want to flatten the parameter space and move from constrained parameters to (-∞, ∞)
# support we can use the `asflat` function

fpost = asflat(post)

# These transformed posterior expect a vector of parameters. That is we can evaluate the
# transformed log density by calling

logdensityof(cpost, rand(rng, dimension(cpost)))
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

sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=50_000);

# The sol vector is in the transformed space, so first we need to transform back to parameter space
# to that we can interpret the solution.

xopt = transform(fpost, sol)

# Given this we can now plot the optimal image or the *maximum a posteriori* (MAP) image.

import CairoMakie as CM
g = imagepixels(μas2rad(200.0), μas2rad(200.0), 256, 256)
fig, ax, plt = CM.image(g, model(xopt); axis=(xreversed=true, aspect=1, xlabel="RA (μas)", ylabel="Dec (μas)"), figure=(;resolution=(650,500),) ,colormap=:afmhot)

# ### Quantifying the Uncertainty of the Reconstruction

# While finding the optimal image is often helpful, in science, the most important thing is to
# quantify the certainty of our inferences. This is the goal of Comrade. In the language
# of Bayesian statistics, we want to find a representation of the posterior of possible image
# reconstructions given our choice of model and the data.
#
# Comrade provides several sampling and other posterior approximation tools. To see the
# list, please see the Libraries section of the docs. For this example, we will be using
# [Pigeons.jl](https://github.com/Julia-Tempering/Pigeons.jl) which is a state-of-the-art
# parallel tempering sampler that enables global exploration of the posterior. For smaller dimension
# problems (< 100) we recommend using this sampler especially if you have access to > 1 thread/core.
using Pigeons
pt = pigeons(target=cpost, explorer=SliceSampler(), record=[traces, round_trip, log_sum_ratio], n_chains=18, n_rounds=9)
chain = sample_array(cpost, pt)


# That's it! To finish it up we can then plot some simple visual fit diagnostics.

# First to plot the image we call
imgs = intensitymap.(skymodel.(Ref(post), sample(chain, 100)), μas2rad(200.0), μas2rad(200.0), 128, 128)
imageviz(imgs[end], colormap=:afmhot)

# What about the mean image? Well let's grab 100 images from the chain, where we first remove the
# adaptation steps since they don't sample from the correct posterior distribution
meanimg = mean(imgs)
imageviz(meanimg, colormap=:afmhot)

# That looks similar to the EHTC VI, and it took us no time at all!. To see how well the
# model is fitting the data we can plot the model and data products
using Plots
plot(model(xopt), dlcamp, label="MAP")

# We can also plot random draws from the posterior predictive distribution.
# The posterior predictive distribution create a number of synthetic observations that
# are marginalized over the posterior.
p = plot(dlcamp);
uva = [sqrt.(uvarea(dlcamp[i])) for i in 1:length(dlcamp)]
for i in 1:10
    m = simulate_observation(post, sample(chain, 1)[1])[1]
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


# For a real run we should also check that the MCMC chain has converged. For
# this we can use MCMCDiagnosticTools
using MCMCDiagnosticTools, Tables
# First, lets look at the effective sample size (ESS) and R̂. This is important since
# the Monte Carlo standard error for MCMC estimates is proportional to 1/√ESS (for some problems)
# and R̂ is a measure of chain convergence. To find both, we can use:
compute_ess(x::NamedTuple) = map(compute_ess, x)
compute_ess(x::AbstractVector{<:Number}) = ess_rhat(x)
compute_ess(x::AbstractVector{<:Tuple}) = map(ess_rhat, Tables.columns(x))
compute_ess(x::Tuple) = map(compute_ess, x)
essrhat = compute_ess(Tables.columns(chain))
# Here, the first value is the ESS, and the second is the R̂. Note that we typically want R̂ < 1.01
# for all parameters, but you should also be running the problem at least four times from four different
# starting locations. In the future we will write an extension that works with Arviz.jl.

# In our example here, we see that we have an ESS > 100 for all parameters and the R̂ < 1.01
# meaning that our MCMC chain is a reasonable approximation of the posterior. For more diagnostics, see
# [`MCMCDiagnosticTools.jl`](https://turinglang.github.io/MCMCDiagnosticTools.jl/stable/).
