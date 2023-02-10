# # Hybrid Imaging of a Black Hole

# In this tutorial we will use **hybrid imaging** to analyze the 2017 EHT data.
# By hybrid imaging we mean decomposing the model into simple geometric models, e.g., rings
# and such, plus a rasterized image model to soak up additional structure in the image.
# This was first developed in [`BB20`](https://iopscience.iop.org/article/10.3847/1538-4357/ab9c1f)
# and applied to EHT 2017 data. In this work we will use a modified model to analyze the data.

# ## Introduction to Hybrid modeling and imaging
# The benfits of using a hybrid based modeling approach is the effective compression of
# information/parameters when fitting the data. Hybrid modeling requires the user to
# incorporate specific knowledge of how you expect the source to look like. For instance
# for M87 we expect the image to be dominated by a ring-like structure. Therefore, rather
# than using a high-dimensional raster to recover the ring we can use a ring model plus
# a very low-dimensional or large pixel size raster to soak up the rest of the parameters.
# This is the approach we will take in this tutorial to analyze the April 11 2017 EHT data
# of M87.


# ## Loading the Data

# To get started we will load Comrade
# ## Load the Data
using Pkg; Pkg.activate(@__DIR__)

using Comrade


# To download the data visit https://doi.org/10.25739/g85n-f134
# To load the eht-imaging obsdata object we do:
obs = load_ehtim_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))

# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      coherent.
obs = scan_average(obs).add_fractional_noise(0.02)

# For this tutorial we will stick to fitting closure only data, although we can
# get better results by also modeling gains, since closure only modeling is equivalent
# to assuming infinite gain priors.
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)

# ## Building the Model/Posterior

# Now we must build our intensity/visibility model. That is, the model that takes in a
# named tuple of parameters and perhaps some metadata required to construct the model.
# For our model we will be using a raster or `ContinuousImage` model plus a `m-ring` model
# and a large asymmetric gaussian component to model the unresolved short-baseline flux.

function model(θ, metadata)
    (;c, f, r, σ, ma, mp, fg, σg, τg, ξg) = θ
    (; grid, cache) = metadata
    img = IntensityMap(f*c, grid)
    mimg = ContinuousImage(img, cache)
    s,c = sincos(mp)
    α = ma*c
    β = ma*s
    ring = (1-f)*smoothed(stretched(MRing(α, β), r, r),σ)
    gauss = fg*rotated(stretched(Gaussian(), σg, σg*(1+τg)), ξg)
    return (1-fg)*(ring + mimg) + gauss
end

# Before we move on let's go into the `model` function a bit. This function takes two arguments
# `θ` and `metadata`. The `θ` argument is a named tuple of parameters that will be fit directly
# to the data. The `metadata` argument is all the ancillary information we need to construct the model.
# For our hybrid model we will need 2 variables for the metadata which will be defined below.
# Note that unlike geometric model fitting here we also use the [`modelimage`](@ref) function.
# This wraps our `ContinuousImage` model with additional information contained in the cache
# that allows us to move from our grid of image fluxes to the model visibilities.
# To combine the models we use `Comrade`'s overloaded `+` operators which will combine the
# images such that there intensities and visibilities are added pointwise.

# Now let's define our metadata. First we will define the cache for the image. This is
# required to compute the numerical Fourier transform.
fovxy  = μas2rad(90.0)
npix   = 6
grid   = imagepixels(fovxy, fovxy, npix, npix)
buffer = IntensityMap(zeros(npix,npix), grid)
# For our image we will be using the
# discrete Fourier transform (`DFTAlg`) since we are using such a small dimensional image. For
# larger rasters (8x8 and above) we recommend to use `NFFTAlg` instead of `DFFTAlg`
# due to the improved scaling. The last argument to the `create_cache` call is the image
# *kernel* or *pulse* and is defines the continuous function that we convolve our image with
# to produce a continuous on sky-image.
cache  = create_cache(DFTAlg(dlcamp), buffer, BSplinePulse{3}())

# Now we form the metadata
metadata = (;grid, cache)

# This is everything we need to form our likelihood, note the first two arguments must be
# the model and then the metadata for the likelihood. The rest of the arguments are required
# to be [`EHTObservations`](@ref)
lklhd = RadioLikelihood(model, metadata, dlcamp, dcphase)

# This forms our model. The next step is defining our image priors.
# For our raster `c` we will use a *Dirichlet* prior which is a multivariate prior
# that exists on the simples. That is, the sum of all the numbers from a `Dirichlet`
# distribution always equals unity. The first paramter is the concentration parameter `α`.
# As `α→0` the images tend to become very sparse, while for `α >> 1` the images tend to
# have uniform brightness. The `α=1` distribution is then just the uniform distribution
# on the simplex. For our work here we use the uniform distribution.

# !!! Warning
#    As α gets small sampling get very difficult and quite multimodals due to the nature
#    of the sparsity prior, be careful when checking convergence when using such a prior.

using VLBIImagePriors
using Distributions
prior = (
          c  = ImageDirichlet(1.0, npix, npix),
          f  = Uniform(0.0, 1.0),
          r  = Uniform(μas2rad(10.0), μas2rad(30.0)),
          σ  = Uniform(μas2rad(0.5), μas2rad(20.0)),
          ma = Uniform(0.0, 0.5),
          mp = Uniform(0.0, 2π),
          fg = Uniform(0.2, 1.0),
          σg = Uniform(μas2rad(50.0), μas2rad(500.0)),
          τg = Uniform(0.0, 1.0),
          ξg = Uniform(0, π)
        )

# This is everything we need to specify our posterior distribution, which our is the main
# object of interest in image reconstructions when using Bayesian inference.
post = Posterior(lklhd, prior)

# To sample from our prior we can do
xrand = prior_sample(post)
# and then plot the results
using Plots
img = intensitymap(model(xrand, metadata), 1.5*fovxy, 1.5*fovxy, 128, 128)
plot(img, title="Random sample")

# ## Reconstructing the Image

# To sample from this posterior it is convienent to first move from our constrained paramter space
# to a unconstrained one (i.e., the support of the transformed posterior is (-∞, ∞)). This is
# done using the `asflat` function.
tpost = asflat(post)

# We can now also find the dimension of our posterior, or the number of parameters we are going to sample.
# !!! Warning
#    This can often be different from what you would expect. This is especially true when using
#    angular variables where to make sampling easier we often artifically increase the dimension
#    of the parameter space.
ndim = dimension(tpost)

# Now we optimize. First we will use BlackBoxOptim which is a genetic algorithm to get us
# in the region of the best fit model.
using ComradeOptimization
using OptimizationBBO
f = OptimizationFunction(tpost, Optimization.AutoForwardDiff())
prob = OptimizationProblem(f, prior_sample(tpost), nothing, lb=fill(-5.0, ndim), ub=fill(5.0,ndim))
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=100_000)

# Alright now we can zoom to the peak!
using OptimizationOptimJL
prob = OptimizationProblem(f, sol.u, nothing)
ℓ = logdensityof(tpost)
sol = solve(prob, LBFGS(), maxiters=1_000, callback=((x,p)->(@info ℓ(x);false)), g_tol=1e-1)

# Before we analyze our solution we first need to transform back to parameter space.
xopt = transform(tpost, sol)

# First we will evaluate our fit by plotting the residuals
residual(model(xopt, metadata), dlcamp)
residual(model(xopt, metadata), dcphase)

# These look pretty reasonable, although maybe they are a bit high. This could probably be
# improved in a few ways, but that is beyond the goal of this quick tutorial.
# Plotting the image we see that we have a ring and an image that looks like a sharper version
# of the original M87 image. This is because we used a more physically motivated model, namely assuming that
# the image should have a ring component in it.
img = intensitymap(model(xopt, metadata), 1.5*fovxy, 1.5*fovxy, 128, 128)
plot(img, title="MAP Image")

# now we sample using hmc
using ComradeAHMC
metric = DiagEuclideanMetric(ndim)
chain, stats = sample(post, AHMC(;metric, autodiff=Comrade.AD.ForwardDiffBackend()), 2500; nadapts=1500, init_params=xopt)

# this took 25 minutes on my desktop which has a Ryzen 7950x cpu. Let's analyze the resulting
# approximate posterior samples

# First we plot the mean image and standard deviation images.
# To do this we first clip the first 1,500 MCMC steps since that is during tuning and
# so the posterior is not sampling from the correct stationary distribution.
using StatsBase
msamples = model.(chain[1501:5:end], Ref(metadata))

# The mean image is then given by
imgs = intensitymap.(msamples, 1.5*fovxy, 1.5*fovxy, 128, 128)
plot(mean(imgs), title="Mean Image")
plot(std(imgs), title="Std Dev.")

# We can also split up the model into its components and analyze each separately
comp = Comrade.components.(Comrade.basemodel.(first.(Comrade.components.(msamples))))
ring_samples = first.(comp)
rast_samples = last.(comp)
ring_imgs = intensitymap.(ring_samples, fovxy, fovxy, 128, 128)
rast_imgs = intensitymap.(rast_samples, fovxy, fovxy, 128, 128)

ring_mean, ring_std = mean_and_std(ring_imgs)
rast_mean, rast_std = mean_and_std(rast_imgs)

p1 = plot(ring_mean, title="Ring Mean", clims=(0.0, maximum(ring_mean)), colorbar=:none)
p2 = plot(ring_std, title="Ring Std. Dev.", clims=(0.0, maximum(ring_mean)), colorbar=:none)
p3 = plot(rast_mean, title="Raster Mean", clims=(0.0, maximum(ring_mean)), colorbar=:none)
p4 = plot(rast_std,  title="Raster Std. Dev.", clims=(0.0, maximum(ring_mean)), colorbar=:none)

plot(p1,p2,p3,p4, layout=(2,2), size=(650, 650))

# Finally, let's take a look at some of the ring parameters
using StatsPlots
p1 = density(rad2μas(chain.r)*2, xlabel="Ring Diameter (μas)")
p2 = density(rad2μas(chain.σ)*2*sqrt(2*log(2)), xlabel="Ring FWHM (μas)")
p3 = density(-rad2deg.(chain.mp) .+ 360.0, xlabel = "Ring PA (deg) E of N")
p4 = density(2*chain.ma, xlabel="Brightness asymmetry")
p5 = density(1 .- chain.f, xlabel="Ring flux fraction")
plot(p1, p2, p3, p4, p5, size=(900, 600), legend=nothing)

# This is very consistent with the original M87 results and it only took 20 minutes compared to the week it used
# to take using old imaging tools.

# ## Computing information
# ```
# Julia Version 1.8.5
# Commit 17cfb8e65ea (2023-01-08 06:45 UTC)
# Platform Info:
#   OS: Linux (x86_64-linux-gnu)
#   CPU: 32 × AMD Ryzen 9 7950X 16-Core Processor
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-13.0.1 (ORCJIT, znver3)
#   Threads: 1 on 32 virtual cores
# Environment:
#   JULIA_EDITOR = code
#   JULIA_NUM_THREADS = 1
# ```
