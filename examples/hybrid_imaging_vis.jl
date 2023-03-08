# # Hybrid Imaging of a Black Hole

# In this tutorial, we will use **hybrid imaging** to analyze the 2017 EHT data.
# By hybrid imaging, we mean decomposing the model into simple geometric models, e.g., rings
# and such, plus a rasterized image model to soak up the additional structure.
# This approach was first developed in [`BB20`](https://iopscience.iop.org/article/10.3847/1538-4357/ab9c1f)
# and applied to EHT 2017 data. We will use a similar model in this tutorial.

# ## Introduction to Hybrid modeling and imaging
# The benefit of using a hybrid-based modeling approach is the effective compression of
# information/parameters when fitting the data. Hybrid modeling requires the user to
# incorporate specific knowledge of how you expect the source to look like. For instance
# for M87, we expect the image to be dominated by a ring-like structure. Therefore, instead
# of using a high-dimensional raster to recover the ring, we can use a ring model plus
#, a very low-dimensional or large pixel size raster to soak up the rest of the emission.
# This is the approach we will take in this tutorial to analyze the April 6 2017 EHT data
# of M87.


# ## Loading the Data

# To get started we will load Comrade
using Comrade

# ## Load the Data
using Pkg #hide
Pkg.activate(joinpath(dirname(pathof(Comrade)), "..", "examples")) #hide

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(1234)


# To download the data visit https://doi.org/10.25739/g85n-f134
# To load the eht-imaging obsdata object we do:
obs = load_ehtim_uvfits(joinpath(dirname(pathof(Comrade)), "..", "examples", "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))

# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      coherent.
obs = scan_average(obs).add_fractional_noise(0.01)

# For this tutorial we will stick to fitting closure only data, although we can
# get better results by also modeling gains, since closure only modeling is equivalent
# to assuming infinite gain priors.
dvis = extract_vis(obs)

# ## Building the Model/Posterior

# Now we build our intensity/visibility model. That is, the model that takes in a
# named tuple of parameters and perhaps some metadata required to construct the model.
# For our model, we will use a raster or `ContinuousImage` model, an `m-ring` model,
# and a large asymmetric Gaussian component to model the unresolved short-baseline flux.

function model(θ, metadata)
    (;c, f, r, σ, ma, mp, fg, σg, τg, ξg, lgamp, gphase) = θ
    (; grid, cache, gcache) = metadata
    ## Form the image model
    ## We multiple the flux by 1.1 since that is the measured total flux of M87
    img = IntensityMap(1.1*f*(1-fg)*c, grid)
    mimg = ContinuousImage(img, cache)
    ## Form the ring model
    s,c = sincos(mp)
    α = ma*c
    β = ma*s
    ring = (1.1*(1-f)*(1-fg))*smoothed(stretched(MRing(α, β), r, r),σ)
    gauss = (1.1*fg)*rotated(stretched(Gaussian(), σg, σg*(1+τg)), ξg)
    m = mimg + (ring + gauss)
    ## Construct the gain model
    gvis = exp.(lgamp .+ 1im.*gphase)
    jg = jonesStokes(gvis, gcache)
    return JonesModel(jg, m)
end

# Before we move on, let's go into the `model` function a bit. This function takes two arguments
# `θ` and `metadata`. The `θ` argument is a named tuple of parameters that are fit
# to the data. The `metadata` argument is all the ancillary information we need to construct the model.
# For our hybrid model, we will need two variables for the metadata, a `grid` that specifies
# the locations of the image pixels and a `cache` that defines the algorithm used to calculate
# the visibilities given the image model. This is required since `ContinuousImage` is most easily
# computed using number Fourier transforms like the [`NFFT`](https://github.com/JuliaMath/NFFT.jl)
# or [FFT](https://github.com/JuliaMath/FFTW.jl).
# To combine the models, we use `Comrade`'s overloaded `+` operators, which will combine the
# images such that their intensities and visibilities are added pointwise.

# Now let's define our metadata. First we will define the cache for the image. This is
# required to compute the numerical Fourier transform.
fovxy  = μas2rad(90.0)
npix   = 6
grid   = imagepixels(fovxy, fovxy, npix, npix)
buffer = IntensityMap(zeros(npix,npix), grid)
# For our image, we will use the
# discrete Fourier transform (`DFTAlg`) since we use such a small dimensional image. For
# larger rasters (8x8 and above) we recommend using `NFFTAlg` instead of `DFFTAlg`
# due to the improved scaling. The last argument to the `create_cache` call is the image
# *kernel* or *pulse* defines the continuous function we convolve our image with
# to produce a continuous on-sky image.
cache  = create_cache(DFTAlg(dvis), buffer, BSplinePulse{3}())

# Now we construct our gain segmenting cache
gcache = jonescache(dvis, ScanSeg())

# Now we form the metadata
metadata = (;grid, cache, gcache)

# This is everything we need to form our likelihood. Note the first two arguments must be
# the model and then the metadata for the likelihood. The rest of the arguments are required
# to be [`Comrade.EHTObservation`](@ref)
lklhd = RadioLikelihood(model, metadata, dvis)

# Moving onto our prior, we first focus on the instrument model priors.
# Each station requires its own prior on both the amplitudes and phases.
# For the amplitudes
# we assume that the gains are apriori well calibrated around unit gains (or 0 log gain amplitudes)
# which corresponds to no instrument corruption. The gain dispersion is then set to 10% for
# all stations except LMT, representing that we expect 10% deviations from scan-to-scan. For LMT
# we let the prior expand to 100% due to the known pointing issues LMT had in 2017.
using Distributions
using DistributionsAD
distamp = (AA = Normal(0.0, 0.1),
           AP = Normal(0.0, 0.1),
           LM = Normal(0.0, 1.0),
           AZ = Normal(0.0, 0.1),
           JC = Normal(0.0, 0.1),
           PV = Normal(0.0, 0.1),
           SM = Normal(0.0, 0.1),
           )

# For the phases, we use a wrapped von Mises prior to respect the periodicity of the variable.
# !!! warning
#     We use AA (ALMA) as a reference station so we use a more restrictied gain prior.
#     We recommend that you do not set the reference prior tighter than this since it
#     makes sampling very difficult.
#-
using VLBIImagePriors
distphase = (AA = DiagonalVonMises(0.0, inv(0.01)),
             AP = DiagonalVonMises(0.0, inv(π^2)),
             LM = DiagonalVonMises(0.0, inv(π^2)),
             AZ = DiagonalVonMises(0.0, inv(π^2)),
             JC = DiagonalVonMises(0.0, inv(π^2)),
             PV = DiagonalVonMises(0.0, inv(π^2)),
             SM = DiagonalVonMises(0.0, inv(π^2)),
           )


# The next step is defining our image priors.
# For our raster `c`, we will use a *Dirichlet* prior, a multivariate prior
# that exists on the simplex. That is, the sum of all the numbers from a `Dirichlet`
# distribution always equals unity. The first parameter is the concentration parameter `α`.
# As `α→0`, the images tend to become very sparse, while for `α >> 1`, the images tend to
# have uniform brightness. The `α=1` distribution is the uniform distribution
# on the simplex. For our work here, we use the uniform simplex distribution.


# !!! Warning
#    As α gets small sampling, it gets very difficult and quite multimodal due to the nature
#    of the sparsity prior, be careful when checking convergence when using such a prior.

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
          ξg = Uniform(0, π),
          lgamp = CalPrior(distamp, gcache),
          gphase = CalPrior(distphase, gcache)
        )

# This is everything we need to specify our posterior distribution, which our is the main
# object of interest in image reconstructions when using Bayesian inference.
post = Posterior(lklhd, prior)

# To sample from our prior we can do
xrand = prior_sample(rng, post)
# and then plot the results
using Plots
img = intensitymap(model(xrand, metadata), μas2rad(120.0), μas2rad(120.0), 128, 128)
plot(img, title="Random sample")

# ## Reconstructing the Image

# To sample from this posterior, it is convenient to first move from our constrained parameter space
# to an unconstrained one (i.e., the support of the transformed posterior is (-∞, ∞)). This is
# done using the `asflat` function.
tpost = asflat(post)

# We can now also find the dimension of our posterior or the number of parameters we will sample.
# !!! Warning
#    This can often be different from what you would expect. This is especially true when using
#    angular variables, where we often artificially increase the dimension
#    of the parameter space to make sampling easier.
ndim = dimension(tpost)

# Now we optimize.
using OptimizationOptimJL
f = OptimizationFunction(tpost, Optimization.AutoZygote())
prob = Optimization.OptimizationProblem(f, prior_sample(rng, tpost), nothing)
ℓ = logdensityof(tpost)
sol = solve(prob, LBFGS(), maxiters=5_000, g_tol=1e-1)

# Before we analyze our solution we first need to transform back to parameter space.
xopt = transform(tpost, sol)

# First we will evaluate our fit by plotting the residuals
residual(model(xopt, metadata), dvis)

# These look reasonable, although they are a bit high. This could be
# improved in a few ways, but that is beyond the goal of this quick tutorial.
# Plotting the image, we see that we have a ring and an image that looks like a sharper version
# of the original M87 image. This is because we used a more physically motivated model by assuming that
# the image should have a ring component.
img = intensitymap(model(xopt, metadata), 1.5*fovxy, 1.5*fovxy, 128, 128)
plot(img, title="MAP Image")

# now we sample using hmc
using ComradeAHMC
metric = DiagEuclideanMetric(ndim)
chain, stats = sample(rng, post, AHMC(;metric, autodiff=Val(:Zygote)), 1000; nadapts=700, init_params=xopt)

# !!! warning
#     This should be run for likely 2-3x more steps to properly estimate expectations of the posterior
#-

# Now lets plot the mean image and standard deviation images.
# To do this we first clip the first 250 MCMC steps since that is during tuning and
# so the posterior is not sampling from the correct stationary distribution.
using StatsBase
msamples = model.(chain[701:2:end], Ref(metadata))

# The mean image is then given by
imgs = intensitymap.(msamples, 1.5*fovxy, 1.5*fovxy, 128, 128)
plot(mean(imgs), title="Mean Image")
plot(std(imgs), title="Std Dev.")

# We can also split up the model into its components and analyze each separately
comp = Comrade.components.(Comrade.basemodel.(first.(Comrade.components.(msamples))))
ring_samples = getindex.(comp,2)
rast_samples = first.(comp)
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

# And viola hybrid imaging is now simple and takes under 10 minutes

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
