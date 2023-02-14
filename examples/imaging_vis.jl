# # Stokes I Simultaneous Image and Instrument Modeling

# In this tutorial, we will create a preliminary reconstruction of the 2017 M87 data on April 6
# by simultaneously creating an image and model for the instrument. By instrument model, we
# mean something akin to self-calibration in traditional VLBI imaging terminology. However,
# unlike traditional self-cal, we will at each point in our parameter space effectively explore
# the possible self-cal solutions. This will allow us to constrain and marginalize over the
# instrument effects, such as time variable gains.

# To get started we load Comrade.
using Comrade


using Pkg #hide
Pkg.activate(joinpath(dirname(pathof(Comrade)), "..", "examples")) #hide


# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(123)



# ## Load the Data


# To download the data visit https://doi.org/10.25739/g85n-f134
# First we will load our data:
obs = load_ehtim_uvfits(joinpath(dirname(pathof(Comrade)), "..", "examples", "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))

# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      coherent.
#   - Add 1% systematic noise to deal with calibration issues that cause 1% non-closing errors.
obs = scan_average(obs.flag_uvdist(0.1e9).add_fractional_noise(0.01))

# Now we extract our complex visibilities.
dvis = extract_vis(obs)

# ##Building the Model/Posterior

# Now, we must build our intensity/visibility model. That is, the model that takes in a
# named tuple of parameters and perhaps some metadata required to construct the model.
# For our model, we will use a raster or `ContinuousImage` for our image model.
# Unlike other imaging examples
# (e.g., [Imaging a Black Hole using only Closure Quantities](@ref)) we also need to include
# a model for the instrument, i.e., gains as well. The gains will be broken into two components
#   - Gain amplitudes which are typically known to 10-20%, except for LMT, which has amplitudes closer to 50-100%.
#   - Gain phases which are more difficult to constrain and can shift rapidly.
# The model is given below:

function model(θ, metadata)
    (;c, lgamp, gphase) = θ
    (; grid, cache, gcache, gcacher) = metadata
    ## Construct the image model we fix the flux to 0.6 Jy in this case
    img = IntensityMap(0.6.*c, grid)
    m = ContinuousImage(img,cache)
    ## Now form our instrument model
    jp = @fastmath jonesStokes(exp, gphase.*1im, gcacher)
    jg = @fastmath jonesStokes(exp, lgamp, gcache)
    return JonesModel(jp*jg, m)
end

# The model construction is very similar to [Imaging a Black Hole using only Closure Quantities](@ref),
# except we fix the compact flux to 0.6 Jy for simplicity in this run. For more information about the image model
# please read the closure-only example. Let's discuss the instrument model [`JonesModel`](@ref).
# Thanks to the EHT pre-calibration, the gains are stable over scans. Therefore, we can
# model the gains on a scan-by-scan basis. To form the instrument model, we need our
#   1. Our (log) gain amplitudes and phases are given below by `lgamp` and `gphase`
#   2. Our function or cache that maps the gains from a list to the stations they impact `gcache.`
#   3. The set of [`JonesPairs`](@ref) produced by [`jonesStokes`](@ref)
# These three ingredients then specify our instrument model. The instrument model can then be
# combined with our image model `cimg` to form the total `JonesModel`.




# Now, let's set up our image model. The EHT's nominal resolution is 20-25 μas. Additionally,
# the EHT is not very sensitive to a larger field of view. Typically 60-80 μas is enough to
# describe the compact flux of M87. Given this, we only need to use a small number of pixels
# to describe our image.
npix = 8
fovxy = μas2rad(67.5)

# Now let's form our cache's. First, we have our usual image cache which is needed to numerically
# compute the visibilities.
grid = imagepixels(fovxy, fovxy, npix, npix)
buffer = IntensityMap(zeros(npix, npix), grid)
cache = create_cache(NFFTAlg(dvis), buffer, BSplinePulse{3}())
# Second, we now construct our instrument model cache. This tells us how to map from the gains
# to the model visibilities. However, to construct this map, we also need to specify the observation
# segmentation over which we expect the gains to change. This is specified in the second argument
# to `jonescache`, and currently, there are two options
#   - `ScanSeg()`: which forces the corruptions to only change from scan-to-scan
#   - `TrackSeg()`: which forces the corruptions to be constant over a night's observation
# For this work, we use the scan segmentation since that is roughly the timescale we expect the
# complex gains to vary. Finally, we model the phases and amplitudes separately. For the phases
# we use a segmented cache where ϕₜ₋₁ + ϵₜ gives the gain at timestamp t, ϕₜ, and our
# model parameters are the ϵₜ. This is specified in [`jonescache`](@ref) by a third argument
# that specifies whether to use this segmented prior.
gcacher = jonescache(dvis, ScanSeg(), true) ## use segmented prior
gcache = jonescache(dvis, ScanSeg())

# Now we can form our metadata we need to fully define our model.
metadata = (;grid, cache, gcache, gcacher)

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

# For the phases, as mentioned above, we will use a segmented gain prior.
# This means that rather than the parameters
# being directly the gains, we fit the first gain for each site, and then
# the other parameters are the segmented gains compared to the previous time. To model this
#, we break the gain phase prior into two parts. The first is the prior
# for the first observing timestamp of each site, `distphase0`, and the second is the
# prior for segmented gain ϵₜ from time i to i+1, given by `distphase`. For the EHT, we are
# dealing with pre-calibrated data, so often, the gain phase jumps from scan to scan are
# minor. As such, we can put a more informative prior on `distphase`.
using VLBIImagePriors
distphase0 = (AA = DiagonalVonMises(0.0, inv(1e-6)),
             AP = DiagonalVonMises(0.0, inv(π^2)),
             LM = DiagonalVonMises(0.0, inv(π^2)),
             AZ = DiagonalVonMises(0.0, inv(π^2)),
             JC = DiagonalVonMises(0.0, inv(π^2)),
             PV = DiagonalVonMises(0.0, inv(π^2)),
             SM = DiagonalVonMises(0.0, inv(π^2)),
           )

distphase = (AA = DiagonalVonMises(0.0, inv(1e-6)),
             AP = DiagonalVonMises(0.0, inv(0.2^2)),
             LM = DiagonalVonMises(0.0, inv(0.2^2)),
             AZ = DiagonalVonMises(0.0, inv(0.2^2)),
             JC = DiagonalVonMises(0.0, inv(0.2^2)),
             PV = DiagonalVonMises(0.0, inv(0.2^2)),
             SM = DiagonalVonMises(0.0, inv(0.2^2)),
           )


# We can now form our model parameter priors. Like our other imaging examples, we use a
# Dirichlet prior for our image pixels. For the log gain amplitudes, we use the `CalPrior`
# which automatically constructs the prior for the given jones cache `gcache`.
(;X, Y) = grid
prior = (
         c = ImageDirichlet(1.0, npix, npix),
         lgamp = CalPrior(distamp, gcache),
         gphase = CalPrior(distphase0, distphase, gcacher),
        )


# Putting it all together we form our likelihood and posterior objects for optimization and
# sampling.
lklhd = RadioLikelihood(model, metadata, dvis)
post = Posterior(lklhd, prior)

# ## Reconstructing the Image and Instrument Effects

# To sample from this posterior, it is convenient to move from our constrained parameter space
# to an unconstrained one (i.e., the support of the transformed posterior is (-∞, ∞)). This is
# done using the `asflat` function.
tpost = asflat(post)

# We can now also find the dimension of our posterior or the number of parameters we are going to sample.
# !!! Warning
#     This can often be different from what you would expect. This is especially true when using
#     angular variables where we often artificially increase the dimension
#     of the parameter space to make sampling easier.
#-
ndim = dimension(tpost)

# Now we optimize. Unlike other imaging examples here we move straight to gradient optimizers
# due to the higher dimension of the space.
using ComradeOptimization
using OptimizationOptimJL
using Zygote
f = OptimizationFunction(tpost, Optimization.AutoZygote())
prob = Optimization.OptimizationProblem(f, prior_sample(rng, tpost), nothing)
ℓ = logdensityof(tpost)
sol = solve(prob, LBFGS(), maxiters=10_000, g_tol=1e-1, callback=(x,p)->(@info ℓ(x); false) )

# !!! Warning
#    Fitting gains tends to be very difficult, meaning that optimization can take a lot longer.
#    The upside is that we usually get nicer images.

# Before we analyze our solution, we first need to transform back to parameter space.
xopt = transform(tpost, sol)

# First we will evaluate our fit by plotting the residuals
using Plots
residual(model(xopt, metadata), dvis)

# These look reasonable, although there may be some minor overfitting. This could be
# improved in a few ways, but that is beyond the goal of this quick tutorial.
# Plotting the image, we see that we have a much cleaner version of the closure-only image from
# [Imaging a Black Hole using only Closure Quantities](@ref).
img = intensitymap(model(xopt, metadata), fovxy, fovxy, 128, 128)
plot(img, title="MAP Image")


# Because we also fit the instrument model, we can inspect their parameters.
# To do this, `Comrade` provides a `caltable` function that converts the flattened gain parameters
# to a tabular format based on the time and its segmentation.
gt = Comrade.caltable(gcacher, xopt.gphase)
plot(gt, layout=(3,3), size=(600,500))
# The gain phases are pretty random, although much of this is due to us picking a random
# reference station for each scan.

# Moving onto the gain amplitudes, we see that most of the gain variation is within 10% as expected
# except LMT, which has massive variations.
gt = Comrade.caltable(gcache, exp.(xopt.lgamp))
plot(gt, layout=(3,3), size=(600,500))


# To sample from the posterior, we will use HMC, specifically the NUTS algorithm. For information about NUTS,
# see Michael Betancourt's [notes](https://arxiv.org/abs/1701.02434).
# !!! note
#    For our `metric,` we use a diagonal matrix due to easier tuning
#-
# However, due to the need to sample a large number of gain parameters, constructing the posterior
# is rather time-consuming. Therefore, for this tutorial, we will only do a quick preliminary run, and any posterior
# inferences should be appropriately skeptical.
#-
using ComradeAHMC
metric = DiagEuclideanMetric(ndim)
chain, stats = sample(rng, post, AHMC(;metric, autodiff=AD.ZygoteBackend()), 400; nadapts=200, init_params=xopt)
#-
# !!! warning
#     This should be run for likely an order of magnitude more steps to properly estimate expectations of the posterior
#-


# Now that we have our posterior, we can put error bars on all of our plots above.
# Let's start by finding the mean and standard deviation of the gain phases
gphase  = hcat(chain.gphase...)
mgphase = mean(gphase, dims=2)
sgphase = std(gphase, dims=2)

# and now the gain amplitudes
gamp  = exp.(hcat(chain.lgamp...))
mgamp = mean(gamp, dims=2)
sgamp = std(gamp, dims=2)

# Now we can use the measurements package to automatically plot everything with error bars.
# First we create a `caltable` the same way but making sure all of our variables have errors
# attached to them.
using Measurements
gmeas_am = measurement.(mgamp, sgamp)
ctable_am = caltable(gcache, vec(gmeas_am)) # caltable expects gmeas_am to be a Vector
gmeas_ph = measurement.(mgphase, sgphase)
ctable_ph = caltable(gcache, vec(gmeas_ph))

# Now let's plot the phase curves
plot(ctable_ph, layout=(3,3), size=(600,500))
#-
# and now the amplitude curves
plot(ctable_am, layout=(3,3), size=(600,500))

# Finally let's construct some representative image reconstructions.
samples = model.(chain[100:5:end], Ref(metadata))
imgs = intensitymap.(samples, fovxy, fovxy, 128,  128);

mimg = mean(imgs)
simg = std(imgs)
p1 = plot(mimg, title="Mean", clims=(0.0, maximum(mimg)));
p2 = plot(simg,  title="Std. Dev.", clims=(0.0, maximum(mimg)));
p3 = plot(imgs[begin],  title="Draw 1", clims = (0.0, maximum(mimg)));
p4 = plot(imgs[end],  title="Draw 2", clims = (0.0, maximum(mimg)));
plot(p1,p2,p3,p4, layout=(2,2), size=(800,800))

# And viola, you have just finished making a preliminary image and instrument model reconstruction.
# In reality, you should run the `sample` step for many more MCMC steps to get a reliable estimate
# for the reconstructed image and instrument model parameters.

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
