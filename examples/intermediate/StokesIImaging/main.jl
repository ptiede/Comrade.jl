# # Stokes I Simultaneous Image and Instrument Modeling

# In this tutorial, we will create a preliminary reconstruction of the 2017 M87 data on April 6
# by simultaneously creating an image and model for the instrument. By instrument model, we
# mean something akin to self-calibration in traditional VLBI imaging terminology. However,
# unlike traditional self-cal, we will solve for the gains each time we update the image
# self-consistently. This allows us to model the correlations between gains and the image.

# To get started we load Comrade.
import Pkg #hide
__DIR = @__DIR__ #hide
pkg_io = open(joinpath(__DIR, "pkg.log"), "w") #hide
Pkg.activate(__DIR; io=pkg_io) #hide
Pkg.develop(; path=joinpath(__DIR, "..", "..", ".."), io=pkg_io) #hide
Pkg.instantiate(; io=pkg_io) #hide
Pkg.precompile(; io=pkg_io) #hide
close(pkg_io) #hide


using Comrade



using Pyehtim
using LinearAlgebra

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(12)



# ## Load the Data


# To download the data visit https://doi.org/10.25739/g85n-f134
# First we will load our data:
obs = ehtim.obsdata.load_uvfits(joinpath(__DIR, "..", "..", "Data", "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
# obs = ehtim.obsdata.load_uvfits("~/Dropbox (Smithsonian External)/M872021Project/Data/2021/CASA/e21e18/V4/M87_calibrated_b3.uvf+EVPA_rotation+netcal_10savg+flag.uvfits")
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
    (;fg, c, σimg) = θ
    (;ftot, mimg) = metadata
    ## Apply the GMRF fluctuations to the image
    rast = apply_fluctuations(CenteredLR(), mimg, σimg.*c.params)
    pimg = parent(rast)
    @. pimg = (ftot*(1-fg))*pimg
    m = ContinuousImage(rast, BSplinePulse{3}())
    x0, y0 = centroid(m)
    ## Add a large-scale gaussian to deal with the over-resolved mas flux
    g = modify(Gaussian(), Stretch(μas2rad(500.0), μas2rad(500.0)), Renormalize(ftot*fg))
    return shifted(m, -x0, -y0) + g
end


# Now, let's set up our image model. The EHT's nominal resolution is 20-25 μas. Additionally,
# the EHT is not very sensitive to a larger field of view. Typically 60-80 μas is enough to
# describe the compact flux of M87. Given this, we only need to use a small number of pixels
# to describe our image.
npix = 32
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
fwhmfac = 2*sqrt(2*log(2))
mpr  = modify(Gaussian(), Stretch(μas2rad(50.0)./fwhmfac))
mimg = intensitymap(mpr, grid)


# Now we can form our metadata we need to fully define our model.
# We will also fix the total flux to be the observed value 1.1. This is because
# total flux is degenerate with a global shift in the gain amplitudes making the problem
# degenerate. To fix this we use the observed total flux as our value.
skymeta = (;ftot = 1.1, mimg = mimg./flux(mimg))


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
         σimg = truncated(Normal(0.0, 0.5); lower=0.0),
         fg = Uniform(0.0, 1.0),
        )

# Now we can construct our sky model.
skym = SkyModel(sky, prior, grid; metadata=skymeta)

# Unlike other imaging examples
# (e.g., [Imaging a Black Hole using only Closure Quantities](@ref)) we also need to include
# a model for the instrument, i.e., gains. The gains will be broken into two components
#   - Gain amplitudes which are typically known to 10-20%, except for LMT, which has amplitudes closer to 50-100%.
#   - Gain phases which are more difficult to constrain and can shift rapidly.

G = SingleStokesGain() do x
    lg = x.lgμ + x.lgσ*x.lgz
    gp = x.gp
    return exp(lg + 1im*gp)
end

intpr = (
    lgμ = ArrayPrior(IIDSitePrior(TrackSeg(), Normal(0.0, 0.2)); LM = IIDSitePrior(TrackSeg(), Normal(0.0, 1.0))),
    lgσ = ArrayPrior(IIDSitePrior(TrackSeg(), Exponential(0.1))),
    lgz = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 1.0))),
    gp= ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))); refant=SEFDReference(0.0), phase=true)
        )
intmodel = InstrumentModel(G, intpr)


# To form the posterior we just combine the skymodel, instrument model and the data. Additionally,
# since we want to use gradients we need to specify the AD mode. Essentially for all modes we recommend
# using `Enzyme.set_runtime_activity(Enzyme.Reverse)`. Eventually as Comrade and Enzyme matures we will 
# no need `set_runtime_activity`.
using Enzyme
post = VLBIPosterior(skym, intmodel, dvis; admode=set_runtime_activity(Enzyme.Reverse))
# done using the `asflat` function.
tpost = asflat(post)
ndim = dimension(tpost)

# We can now also find the dimension of our posterior or the number of parameters we are going to sample.
# !!! warning
#     This can often be different from what you would expect. This is especially true when using
#     angular variables where we often artificially increase the dimension
#     of the parameter space to make sampling easier.
#-

# To initialize our sampler we will use optimize using Adam
using Optimization
using OptimizationOptimisers
xopt, sol = comrade_opt(post, Optimisers.Adam(); initial_params=prior_sample(rng, post), maxiters=20_000, g_tol=1e-1)

# !!! warning
#     Fitting gains tends to be very difficult, meaning that optimization can take a lot longer.
#     The upside is that we usually get nicer images.
#-
# First we will evaluate our fit by plotting the residuals
using Plots
using DisplayAs
residual(post, xopt) |> DisplayAs.PNG |> DisplayAs.Text

# These look reasonable, although there may be some minor overfitting. This could be
# improved in a few ways, but that is beyond the goal of this quick tutorial.
# Plotting the image, we see that we have a much cleaner version of the closure-only image from
# [Imaging a Black Hole using only Closure Quantities](@ref).
import CairoMakie as CM
CM.activate!(type = "png", px_per_unit=3) #hide
g = imagepixels(fovx, fovy, 128, 128)
img = intensitymap(skymodel(post, xopt), g)
imageviz(img, size=(500, 400))|> DisplayAs.PNG |> DisplayAs.Text



# Because we also fit the instrument model, we can inspect their parameters.
# To do this, `Comrade` provides a `caltable` function that converts the flattened gain parameters
# to a tabular format based on the time and its segmentation.
intopt = instrumentmodel(post, xopt)
gt = Comrade.caltable(angle.(intopt))
plot(gt, layout=(3,3), size=(600,500)) |> DisplayAs.PNG |> DisplayAs.Text

# The gain phases are pretty random, although much of this is due to us picking a random
# reference sites for each scan.

# Moving onto the gain amplitudes, we see that most of the gain variation is within 10% as expected
# except LMT, which has massive variations.
gt = Comrade.caltable(abs.(intopt))
plot(gt, layout=(3,3), size=(600,500)) |> DisplayAs.PNG |> DisplayAs.Text


# To sample from the posterior, we will use HMC, specifically the NUTS algorithm. For
# information about NUTS,
# see Michael Betancourt's [notes](https://arxiv.org/abs/1701.02434).
# However, due to the need to sample a large number of gain parameters, constructing the posterior
# is rather time-consuming. Therefore, for this tutorial, we will only do a quick preliminary
# run
#-
using AdvancedHMC
chain = sample(rng, post, NUTS(0.8), 1_000; n_adapts=500, progress=false, initial_params=xopt)
#-
# !!! note
#     The above sampler will store the samples in memory, i.e. RAM. For large models this
#     can lead to out-of-memory issues. To fix that you can include the keyword argument
#     `saveto = DiskStore()` which periodically saves the samples to disk limiting memory
#     useage. You can load the chain using `load_samples(diskout)` where `diskout` is
#     the object returned from sample.
#-



# Now we prune the adaptation phase
chain = chain[501:end]

#-
# !!! warning
#     This should be run for likely an order of magnitude more steps to properly estimate expectations of the posterior
#-


# Now that we have our posterior, we can put error bars on all of our plots above.
# Let's start by finding the mean and standard deviation of the gain phases
mchain = Comrade.rmap(mean, chain)
schain = Comrade.rmap(std, chain)
# Now we can use the measurements package to automatically plot everything with error bars.
# First we create a `caltable` the same way but making sure all of our variables have errors
# attached to them.
using Measurements
gmeas = instrumentmodel(post, (;instrument= map((x,y)->Measurements.measurement.(x,y), mchain.instrument, schain.instrument)))
ctable_am = caltable(abs.(gmeas))
ctable_ph = caltable(angle.(gmeas))

# Now let's plot the phase curves
plot(ctable_ph, layout=(4,3), size=(600,500)) |> DisplayAs.PNG |> DisplayAs.Text
#-
# and now the amplitude curves
plot(ctable_am, layout=(4,3), size=(600,500)) |> DisplayAs.PNG |> DisplayAs.Text

# Finally let's construct some representative image reconstructions.
samples = skymodel.(Ref(post), chain[begin:5:end])
imgs = intensitymap.(samples, Ref(g))

mimg = mean(imgs)
simg = std(imgs)
fig = CM.Figure(;resolution=(700, 700));
axs = [CM.Axis(fig[i, j], xreversed=true, aspect=1) for i in 1:2, j in 1:2]
CM.image!(axs[1,1], mimg, colormap=:afmhot); axs[1, 1].title="Mean"
CM.image!(axs[1,2], simg./(max.(mimg, 1e-8)), colorrange=(0.0, 2.0), colormap=:afmhot);axs[1,2].title = "Std"
CM.image!(axs[2,1], imgs[1],   colormap=:afmhot);
CM.image!(axs[2,2], imgs[end], colormap=:afmhot);
CM.hidedecorations!.(axs)
fig |> DisplayAs.PNG |> DisplayAs.Text



# And viola, you have just finished making a preliminary image and instrument model reconstruction.
# In reality, you should run the `sample` step for many more MCMC steps to get a reliable estimate
# for the reconstructed image and instrument model parameters.
