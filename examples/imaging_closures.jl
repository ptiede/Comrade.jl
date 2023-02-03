# # Imaging a Black Hole using only Closure Quantities

# In this tutorial we will create a preliminary reconstruction of the 2017 M87 data on April 6
# using closure only imaging. This serves as a general introduction to imaging in Comrade,
# but it ignores the need to also model the instrument. For an introduction to simulaneous
# image and instrument modeling see [`Stokes I Simultaneuous Image and Instrument Modeling`](@ref)


# ## Introduction to Closure Imaging
# The EHT is the highest resolution telescope ever created. Its resolution is equivalent
# to roughly tracking a hockey puck on the moon when viewing it from earth. However,
# the EHT is also a unique interferometer. For one the data is produces is incredible sparse.
# The telescope it only from 8 geographic locations around the planet, each with its unique
# telescope. Additionally, the EHT observes at a much higher frequency than typical interferometers.
# As a result, the EHT data is often poorly calibrated. Meaning there can be large instrumental effects
# often called *gains* that can corrupt our signal. One way to deal with this is to fit quantities
# that are independent of gains, these are often called **closure quantities**. The types of
# closure quantities are briefly described in [Introduction to the VLBI Imaging Problem](@ref).
#
# I this tutorial we will do closure only modeling of M87 to produce preliminary images of M87.
# Note that for publication ready images we really should also consider fitting complex visibilities.

using Pkg; Pkg.activate(@__DIR__)

# To get started we will load Comrade
using Comrade

# ## Load the Data
# To download the data visit https://doi.org/10.25739/g85n-f134
# To load the eht-imaging obsdata object we do:
obs = load_ehtim_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))

# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      coherent.
#   - Add 1% systematic noise to deal with calibration issues that cause 1% non-closing errors.
obs = scan_average(obs).add_fractional_noise(0.015)

# Now we extract our closure quantities from the EHT data set.
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)

# ## Build the Model/Posterior
# For our model we will be using a image model that consists of a raster of point sources,
# convolved with some pulse or kernel to make a `ContinuousImage` object with it `Comrade's`
# generic image model. Note that `ContinuousImage(img, cache)` actually creates a `ModelImage`
# object that allows `Comrade` to numerically compute the Fourier transform of the image.
function model(θ, metadata)
    (;c) = θ
    (; grid, cache) = metadata
    # Construct the image model
    img = IntensityMap(c, grid)
    return  ContinuousImage(img, cache)
end


# Now 'et's set up our image model. The EHT's nominal resolution is 20-25 μas. Additionally,
# the EHT is not very sensitive to larger field of views, typically 60-80 μas is enough to
# describe the compact flux of M87. Given this we only need to use a small number of pixels
# to describe our image.
npix = 8
fovxy = μas2rad(75.0)

# Now we can feed in the array information to form the cache. We will be using a DFT since
# it is efficient for so few pixels
# We will use a Dirichlet prior to enforce that the flux sums to unity since closures are
# degenerate to total flux.
grid = imagepixels(fovxy, fovxy, npix, npix)
buffer = IntensityMap(zeros(npix,npix), grid)
cache = create_cache(DFTAlg(dlcamp), buffer, BSplinePulse{3}())
metadata = (;grid, cache)

# Now we need to specify our image prior. For this work we use a very simple Dirichlet prior
using VLBIImagePriors
(;X, Y) = grid
prior = (c = ImageDirichlet(1.0, npix, npix), )

lklhd = RadioLikelihood(model, metadata, dlcamp, dcphase)
post = Posterior(lklhd, prior)

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

# Alright now let's zoom to the peak
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
img = intensitymap(model(xopt, metadata), fovxy, fovxy, 128, 128)
plot(img, title="MAP Image")

# To sample from the posterior we will use HMC and more specifically the NUTS algorithm. For information about NUTS
# see Michael Betancourt's [notes](https://arxiv.org/abs/1701.02434).
# !!! note
#    For our `metric` we use a diagonal matrix due to easier tuning.
using ComradeAHMC
metric = DiagEuclideanMetric(ndim)
chain, stats = sample(post, AHMC(;metric, autodiff=AD.ZygoteBackend()), 3000; nadapts=2000, init_params=xopt)

# !!! note
#    This takes roughly 6 minutes on my machine but depending on your computer your run time may be a bit longer.
#
# Now that we have our posterior we can start to assess which parts of the image is strongly inferred by the
# data. This is rather unique to `Comrade` where more traditional imaging algorithms like CLEAN and RML are inherently
# unable to assess uncertainty in their reconstructions.
#
# To explore our posterior let's first create images from a bunch of draws from the posterior
msamples = model.(chain[2000:10:end], Ref(metadata))

# The mean image is then given by
using StatsBase
imgs = intensitymap.(msamples, fovxy, fovxy, 128, 128)
mimg, simg = mean_and_std(imgs)
p1 = plot(mimg, title="Mean Image")
p2 = plot(simg./mimg, title="1/SNR")
p3 = plot(imgs[1], title="Draw 1")
p4 = plot(imgs[end], title="Draw 2")
plot(p1, p2, p3, p4, size=(800,800), colorbar=:none)

# And viola you have a quick and prelminary image of M87 fitting only closure products.
# For a publication level version we would recommend
#    1. Running the chain longer and multiple time to properly assess things like ESS and R̂ (see [`Making an Image of a Black Hole`](@ref))
#    2. Also fitting gains, typically gain amplitudes are good to 10-20% for the EHT not the infinite uncertainty closures implicitly assume
#    3. Making sure the posterior is unimodel (hint for this example it isn't!). The EHT image posteriors can be quite complicated so typically
#       you want to use a sampler that can deal with multi-modal posteriors. Check out the package [`Pigeons.jl`](https://github.com/Julia-Tempering/Pigeons.jl)
#       for an **in-development** package that should easily enable this type of sampling.



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
