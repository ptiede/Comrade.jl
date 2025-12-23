import Pkg #hide
__DIR = @__DIR__ #hide
pkg_io = open(joinpath(__DIR, "pkg.log"), "w") #hide
Pkg.activate(__DIR; io = pkg_io) #hide
Pkg.develop(; path = joinpath(__DIR, "..", "..", ".."), io = pkg_io) #hide
Pkg.instantiate(; io = pkg_io) #hide
Pkg.precompile(; io = pkg_io) #hide
close(pkg_io) #hide


# # Stokes I Simultaneous Image and Instrument Modeling

# In this tutorial, we will create a preliminary reconstruction of the 2017 M87 data on April 6
# by simultaneously creating an image and model for the instrument. By instrument model, we
# mean something akin to self-calibration in traditional VLBI imaging terminology. However,
# unlike traditional self-cal, we will solve for the gains each time we update the image
# self-consistently. This allows us to model the correlations between gains and the image.

# To get started we load Comrade.
using Comrade


using Pyehtim
using LinearAlgebra

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(42)


# ## Load the Data


# To download the data visit https://doi.org/10.25739/g85n-f134
# First we will load our data:
obs = ehtim.obsdata.load_uvfits(joinpath(__DIR, "..", "..", "Data", "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      coherent.
#   - Add 1% systematic noise to deal with calibration issues that cause 1% non-closing errors.
obs = scan_average(obs).add_fractional_noise(0.02)

# Now we extract our complex visibilities.
dvis = extract_table(obs, Visibilities())

# ##Building the Model/Posterior

# Now, we must build our intensity/visibility model. That is, the model that takes in a
# named tuple of parameters and perhaps some metadata required to construct the model.
# For our model, we will use a raster or `ContinuousImage` for our image model.
# The model is given below:


# The model construction is very similar to [Imaging a Black Hole using only Closure Quantities](@ref),
# except we include a large scale gaussian since we want to model the zero baselines.
# For more information about the image model please read the closure-only example.
function sky(θ, metadata)
    (; fg, c, σimg) = θ
    (; ftot, mimg) = metadata
    ## Apply the GMRF fluctuations to the image
    rast = apply_fluctuations(CenteredLR(), mimg, σimg .* c.params)
    pimg = parent(rast)
    @. pimg = (ftot * (1 - fg)) * pimg
    m = ContinuousImage(rast, BSplinePulse{3}())
    ## Add a large-scale gaussian to deal with the over-resolved mas flux
    g = modify(Gaussian(), Stretch(μas2rad(500.0), μas2rad(500.0)), Renormalize(ftot * fg))
    x, y = centroid(m)
    return shifted(m, -x, -y) + g
end


# Now, let's set up our image model. The EHT's nominal resolution is 20-25 μas. Additionally,
# the EHT is not very sensitive to a larger field of view. Typically 60-80 μas is enough to
# describe the compact flux of M87. Given this, we only need to use a small number of pixels
# to describe our image.
npix = 48
fovx = μas2rad(200.0)
fovy = μas2rad(200.0)

# Now let's form our cache's. First, we have our usual image cache which is needed to numerically
# compute the visibilities.
grid = imagepixels(fovx, fovy, npix, npix)

# Now we need to specify our image prior. For this work we will use a Gaussian Markov
# Random field prior
# Since we are using a Gaussian Markov random field prior we need to first specify our `mean`
# image. This behaves somewhat similary to a entropy regularizer in that it will
# start with an initial guess for the image structure. For this tutorial we will use a
# a symmetric Gaussian with a FWHM of 50 μas
using VLBIImagePriors
using Distributions
fwhmfac = 2 * sqrt(2 * log(2))
mpr = modify(Gaussian(), Stretch(μas2rad(60.0) ./ fwhmfac))
mimg = intensitymap(mpr, grid)


# Now we can form our metadata we need to fully define our model.
# We will also fix the total flux to be the observed value 1.1. This is because
# total flux is degenerate with a global shift in the gain amplitudes making the problem
# degenerate. To fix this we use the observed total flux as our value.
skymeta = (; ftot = 1.1, mimg = mimg ./ flux(mimg))


# To make the Gaussian Markov random field efficient we first precompute a bunch of quantities
# that allow us to scale things linearly with the number of image pixels. The returns a
# functional that accepts a single argument related to the correlation length of the field.
# The second argument defines the underlying random field of the Markov process. Here
# we are using a zero mean and unit variance Gaussian Markov random field.
# For this tutorial we will use the first order random field
cprior = corr_image_prior(grid, dvis)


# Putting everything together the total prior is then our image prior, a prior on the
# standard deviation of the MRF, and a prior on the fractional flux of the Gaussian component.
prior = (
    c = cprior,
    σimg = truncated(Normal(0.0, 0.5); lower = 0.0),
    fg = Uniform(0.0, 1.0),
)

# Now we can construct our sky model.
skym = SkyModel(sky, prior, grid; metadata = skymeta)

# Unlike other imaging examples
# (e.g., [Imaging a Black Hole using only Closure Quantities](@ref)) we also need to include
# a model for the instrument, i.e., gains. The gains will be broken into two components
#   - Gain amplitudes which are typically known to 10-20%, except for LMT, which has amplitudes closer to 50-100%.
#   - Gain phases which are more difficult to constrain and can shift rapidly.

G = SingleStokesGain() do x
    lg = x.lg
    gp = x.gp
    return exp(lg + 1im * gp)
end

intpr = (
    lg = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.2)); LM = IIDSitePrior(ScanSeg(), Normal(0.0, 1.0))),
    gp = ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))); refant = SEFDReference(0.0), phase = true),
)
intmodel = InstrumentModel(G, intpr)


# To form the posterior we just combine the skymodel, instrument model and the data. To utilize
# gradients of the posterior we also need to load `Enzyme.jl`. Under the hood, Comrade will use
# Enzyme to compute the gradients of the posterior.
using Enzyme
post = VLBIPosterior(skym, intmodel, dvis)

# ## Optimization and Sampling

# Now we need to actually compute our image. For this we will first follow standard approaches in
# VLBI and find the maximum a posteriori (MAP) estimate of the image and instrument model.
# This is done using the `comrade_opt` function which accepts the posterior, an optimization
# algorithm, and some keyword arguments. For this tutorial we will use the L-BFGS algorithm.
# For more information about the optimization algorithms available see the Optimization.jl
# [docs](https://docs.sciml.ai/Optimization/stable/).
using Optimization, OptimizationLBFGSB
xopt, sol = comrade_opt(
    post, LBFGSB(); initial_params = prior_sample(rng, post),
    maxiters = 2000, g_tol = 1.0e-1
);

# !!! warning
#     Fitting gains tends to be very difficult, meaning that optimization can take a lot longer.
#     The upside is that we usually get nicer images.
#-
# First we will evaluate our fit by plotting the residuals
using CairoMakie
using DisplayAs
res = residuals(post, xopt)
plotfields(res[1], :uvdist, :res) |> DisplayAs.PNG |> DisplayAs.Text

# These look reasonable, although there may be some minor overfitting. This could be
# improved in a few ways, but that is beyond the goal of this quick tutorial.
# Plotting the image, we see that we have a much cleaner version of the closure-only image from
# [Imaging a Black Hole using only Closure Quantities](@ref).
g = imagepixels(fovx, fovy, 128, 128)
img = intensitymap(skymodel(post, xopt), g)
imageviz(img, size = (500, 400)) |> DisplayAs.PNG |> DisplayAs.Text


# Because we also fit the instrument model, we can inspect their parameters.
# First, let's query the posterior object with the optimal parameters to get
# the instrument model.
intopt = instrumentmodel(post, xopt)
# This returns a `SiteArray` object which contains the gains as a flat vector with
# metadata about the sites, time, and frequency. To visualize the gains we can use
# the `plotcaltable` function which automatically plots the gains. Since the gains are complex
# we will first plot the phases.
plotcaltable(angle.(intopt)) |> DisplayAs.PNG |> DisplayAs.Text

# Due to the a priori calibration of the data, the gain phases are quite stable and just drift
# over time. Note that the one outlier around 2.5 UT is when ALMA leaves the array. As a result
# we shift the reference antenna to APEX on that scan causing the gain phases to jump.


# We can also plot the gain amplitudes
plotcaltable(abs.(intopt)) |> DisplayAs.PNG |> DisplayAs.Text
# Here we find relatively stable gains for most stations. The exception is LMT which
# has a large offset after 2.5 UT. This is a known issue with the LMT in 2017 and is
# due to pointing issues. However, we can see the power of simultaneous imaging and
# instrument modeling as we are able to solve for the gain amplitudes and get a reasonable
# image.


# One problem with the MAP estimate is that it does not provide uncertainties on the image.
# That is we are unable to statistically assess which components of the image are certain.
# Comrade is really a Bayesian imaging and calibration package for VLBI. Therefore, our
# goal is to sample from the posterior distribution of the image and instrument model.
# This is a very high-dimensional distribution with typically 1,000 - 100,000 parameters.
# To sample from this very high dimensional distribution, Comrade has an array of samplers
# that can be used. However, note that `Comrade` also satisfies the `LogDensityProblems.jl`
# interface. Therefore, you can use any package that supports `LogDensityProblems.jl` if you
# have your own fancy sampler.
# -
# For this example, we will use HMC, specifically the NUTS algorithm. For
# However, due to the need to sample a large number of gain parameters, constructing the posterior
# can take a few minutes. Therefore, for this tutorial, we will only do a quick preliminary
# run.
#-
using AdvancedHMC
chain = sample(rng, post, NUTS(0.8), 700; n_adapts = 500, initial_params = xopt)
#-
# !!! note
#     The above sampler will store the samples in memory, i.e. RAM. For large models this
#     can lead to out-of-memory issues. To avoid this we recommend using the
#     `saveto = DiskStore()` kwargs which periodically saves the samples to disk limiting memory
#     useage. You can load the chain using `load_samples(diskout)` where `diskout` is
#     the object returned from sample.
#-


# Now we prune the adaptation phase
chain = chain[begin+500:end]

#-
# !!! warning
#     This should be run for likely an order of magnitude more steps to properly estimate expectations of the posterior
#-


# Now that we have our posterior, we can put error bars on all of our plots above.
# Let's start by finding the mean and standard deviation of the gain phases
mchain = Comrade.rmap(mean, chain);
schain = Comrade.rmap(std, chain);
# Now we can use the measurements package to automatically plot everything with error bars.
# First we create a `caltable` the same way but making sure all of our variables have errors
# attached to them.
using Measurements
gmeas = instrumentmodel(post, (; instrument = map((x, y) -> Measurements.measurement.(x, y), mchain.instrument, schain.instrument)))
ctable_am = caltable(abs.(gmeas))
ctable_ph = caltable(angle.(gmeas))

# Now let's plot the phase curves
plotcaltable(ctable_ph) |> DisplayAs.PNG |> DisplayAs.Text
#-
# and now the amplitude curves
plotcaltable(ctable_am) |> DisplayAs.PNG |> DisplayAs.Text

# Finally let's construct some representative image reconstructions.
samples = skymodel.(Ref(post), chain[begin:5:end])
imgs = intensitymap.(samples, Ref(g))

mimg = mean(imgs)
simg = std(imgs)
fig = Figure(; resolution = (700, 700));
axs = [Axis(fig[i, j], xreversed = true, aspect = 1) for i in 1:2, j in 1:2]
image!(axs[1, 1], mimg, colormap = :afmhot); axs[1, 1].title = "Mean"
image!(axs[1, 2], simg ./ (max.(mimg, 1.0e-8)), colorrange = (0.0, 2.0), colormap = :afmhot);axs[1, 2].title = "Frac. Uncer."
image!(axs[2, 1], imgs[1], colormap = :afmhot);
image!(axs[2, 2], imgs[end], colormap = :afmhot);
hidedecorations!.(axs)
fig |> DisplayAs.PNG |> DisplayAs.Text


# And viola, you have just finished making a preliminary image and instrument model reconstruction.
# In reality, you should run the `sample` step for many more MCMC steps to get a reliable estimate
# for the reconstructed image and instrument model parameters.
