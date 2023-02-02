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

# To get started we will load Comrade
using Comrade

# and eht-imaging so we can load the data
load_ehtim()

# To download the data visit https://doi.org/10.25739/g85n-f134
obs = load_ehtim_uvfits(joinpath(@__DIR__, "SR1_M87_2017_101_lo_hops_netcal_StokesI.uvfits"))

# Now we do some minor preprocessing:
#   1. We kill 0-baselines since we don't care about large scale flux
#   2. Scan average the data since the data have been preprocessed so that the gain phases
#      coherent.
obs = scan_average(obs.flag_uvdist(uv_min=0.1e9)).add_fractional_noise(0.01)

# For this tutorial we will stick to fitting closure only data.
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)

# Now we must build our intensity/visibility model. That is, the model that takes in a
# named tuple of parameters and perhaps some metadata required to construct the model.
# For our model we will be using a raster or `ContinuousImage` model plus a `m-ring` model.
# See XYZ for an introduction to the image model, and XYZ for an introduction to the m-ring
# model.

function model(θ, metadata)
    (;c, f, r, σ, τ, ξ, ma, mp) = θ
    (; grid, pulse, cache) = metadata
    img = IntensityMap(f*c, grid)
    cimg = ContinuousImage(img, pulse)
    mimg = modelimage(cimg, cache)
    s,c = sincos(mp)
    α = ma*c
    β = ma*s
    ring = (1-f)*rotated(smoothed(stretched(MRing(α, β), r, r*τ),σ), ξ)
    return ring + mimg
end

# Before we move on let's go into the `model` function a bit. This function takes two arguments
# `θ` and `metadata`. The `θ` argument is a named tuple of parameters that will be fit directly
# to the data. The `metadata` argument is all the ancillary information we need to construct the model.
# For our hybrid model we will need 3 variables for the metadata which will be defined below.
# Note that unlike geometric model fitting here we also use the [`modelimage`](@ref) function.
# This wraps our `ContinuousImage` model with additional information contained in the cache
# that allows us to move from our grid of image fluxes to the model visibilities.
# To combine the models we use `Comrade`'s overloaded `+` operators which will combine the
# images such that there intensities and visibilities are added pointwise.


# To finish specifying our model we will now specify our metadata
# First we need to define the `grid` variable
# Now form model we are going to fit
fovxy = μas2rad(120.0)
npix = 6
mms = Model(dlcamp, fovxy, npix)


prior = (
          c = ImageDirichlet(0.5, npix, npix),
          f = Uniform(0.0, 1.0),
          r = Uniform(μas2rad(10.0), μas2rad(30.0)),
          σ = Uniform(μas2rad(0.5), μas2rad(20.0)),
          τ = Uniform(0.5, 1.0),
          ξ = Uniform(-π/2, π/2),
          ma = Uniform(0.0, 0.5),
          mp = Uniform(0.0, 2π),
        )

# Now we can build the posterior
post = Posterior(lklhd, prior, mms)

# We will use HMC to explore the posterior. However, we can really help sampling if we first
# get a good starting location. As such we will first optimize. First lets transform the
# posterior parameters to the unconstrained space.
tpost = asflat(post)


ndim = dimension(tpost)

# Now we optimize. First we will use BlackBoxOptim which is a genetic algorithm to get us
# in the region of the best fit model.
f = OptimizationFunction(tpost, Optimization.AutoForwardDiff())
prob = OptimizationProblem(f, randn(ndim), nothing, lb=fill(-5.0, ndim), ub=fill(5.0, ndim))
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=100_000)

# Alright now we can zoom to the peak! But we can even do better, we will also compute
# the laplace approximation to the posterior as a byproduct.
prob = OptimizationProblem(f, sol.u, nothing)
ℓ = logdensityof(tpost)
sol = solve(prob, LBFGS(), maxiters=1_000, callback=((x,p)->(@info ℓ(x);false)), g_tol=1e-1)

# tranform back to parameter space
xopt = transform(tpost, sol)

# plot the residuals and the intensity map
residual(mms(xopt), dlcamp)
residual(mms(xopt), dcphase)
img = intensitymap(mms(xopt), μas2rad(160.0), μas2rad(160.0), 512, 512)
plot(abs.(img) , xlims=(-60.0,60.0), ylims=(-60.0,60.0))

# now we sample using hmc
metric = DenseEuclideanMetric(ndim)
chain, stats = sample(post, AHMC(;metric, autodiff=Comrade.AD.ForwardDiffBackend()), 4000; nadapts=3000, init_params=xopt)

# this took 25 minutes on my laptop which has a 11 gen core i7

# Plot the mean image and standard deviation
# Adaptation ruins detailed balance/reversibility in the chain!
# We will also split the uncertainty and mean maps into ring and raster
using StatsBase
msamples = Comrade.components.(mms.(sample(chain, 500)))
ring_samples = first.(msamples)
rast_samples = last.(msamples)
ring_imgs = intensitymap.(ring_samples, fovxy, fovxy, 256, 256)
rast_imgs = intensitymap.(rast_samples, fovxy, fovxy, 256, 256)

ring_mean, ring_std = mean_and_std(ring_imgs)
rast_mean, rast_std = mean_and_std(rast_imgs)
both_mean, both_std = mean_and_std(ring_imgs .+ rast_imgs)

p1 = plot(ring_mean, title="Ring Mean", clims=(0.0, maximum(ring_mean)))
p2 = plot(ring_std,  title="Ring Std. Dev.", clims=(0.0, maximum(ring_mean)))
p3 = plot(rast_mean, title="Raster Mean", clims=(0.0, maximum(rast_mean)))
p4 = plot(rast_std,  title="Raster Std. Dev.", clims=(0.0, maximum(rast_mean)))
p5 = plot(both_mean, title="Both Mean", clims=(0.0, maximum(both_mean)))
p6 = plot(both_std,  title="Both Std. Dev.", clims=(0.0, maximum(both_mean)))

plot(p1,p2,p3,p4,p5,p6, layout=(3,2), size=(950,1000))

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
