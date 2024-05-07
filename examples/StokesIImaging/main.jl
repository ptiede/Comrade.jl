# # Stokes I Simultaneous Image and Instrument Modeling

# In this tutorial, we will create a preliminary reconstruction of the 2017 M87 data on April 6
# by simultaneously creating an image and model for the instrument. By instrument model, we
# mean something akin to self-calibration in traditional VLBI imaging terminology. However,
# unlike traditional self-cal, we will at each point in our parameter space effectively explore
# the possible self-cal solutions. This will allow us to constrain and marginalize over the
# instrument effects, such as time variable gains.

# To get started we load Comrade.
import Pkg #hide
__DIR = @__DIR__ #hide
pkg_io = open(joinpath(__DIR, "pkg.log"), "w") #hide
Pkg.activate(__DIR; io=pkg_io) #hide
Pkg.develop(; path=joinpath(__DIR, "..", ".."), io=pkg_io) #hide
Pkg.instantiate(; io=pkg_io) #hide
Pkg.precompile(; io=pkg_io) #hide
close(pkg_io) #hide


ENV["GKSwstype"] = "nul" #hide
using Comrade



using Pyehtim
using LinearAlgebra

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(12)



# ## Load the Data


# To download the data visit https://doi.org/10.25739/g85n-f134
# First we will load our data:
obs = ehtim.obsdata.load_uvfits(joinpath(__DIR, "../Data/SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))

# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      coherent.
#   - Add 1% systematic noise to deal with calibration issues that cause 1% non-closing errors.
obs = scan_average(obs).add_fractional_noise(0.01)

# Now we extract our complex visibilities.
dvis = extract_table(obs, ComplexVisibilities())

# ##Building the Model/Posterior

# Now, we must build our intensity/visibility model. That is, the model that takes in a
# named tuple of parameters and perhaps some metadata required to construct the model.
# For our model, we will use a raster or `ContinuousImage` for our image model.
# The model is given below:


function sky(θ, metadata)
    (;fg, c, σimg) = θ
    (;ftot, meanpr, grid) = metadata
    ## Transform to the log-ratio pixel fluxes
    cp = meanpr .+ σimg.*c.params
    ## Transform to image space
    rast = (ftot*(1-fg))*(to_simplex(CenteredLR(), cp))
    m = ContinuousImage(rast, grid, BSplinePulse{3}())
    x0, y0 = centroid(m.img)
    ## Add a large-scale gaussian to deal with the over-resolved mas flux
    g = modify(Gaussian(), Stretch(μas2rad(250.0), μas2rad(250.0)), Renormalize(ftot*fg))
    return shifted(m, -x0, -y0) + g
end

# Unlike other imaging examples
# (e.g., [Imaging a Black Hole using only Closure Quantities](@ref)) we also need to include
# a model for the instrument, i.e., gains as well. The gains will be broken into two components
#   - Gain amplitudes which are typically known to 10-20%, except for LMT, which has amplitudes closer to 50-100%.
#   - Gain phases which are more difficult to constrain and can shift rapidly.

using VLBIImagePriors
using Distributions
G = SingleStokesGain() do x
    lg = x.lg
    gp = x.gp
    return exp(lg + 1im*gp)
end

intpr = (
    # lg0= ArrayPrior(IIDSitePrior(TrackSeg(), Normal(0.0, 0.1)); LM = IIDSitePrior(TrackSeg(), Normal(0.0, 1.0))),
    # lgσ= ArrayPrior(IIDSitePrior(TrackSeg(), Normal(-1.0,1.0)); LM = IIDSitePrior(TrackSeg(), Normal(-1.0, 1.0))),
    lg= ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.2))),
    # gpσ= ArrayPrior(IIDSitePrior(TrackSeg(), Normal(0.0, 1.0)); refant=SEFDReference(0.0)),
    gp= ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))); refant=SEFDReference(0.0))
        )
intmodel = InstrumentModel(G, intpr)

# The model construction is very similar to [Imaging a Black Hole using only Closure Quantities](@ref),
# except we include a large scale gaussian since we want to model the zero baselines.
# For more information about the image model
# please read the closure-only example. Let's discuss the instrument model [`Comrade.JonesModel`](@ref).
# Thanks to the EHT pre-calibration, the gains are stable over scans. Therefore, we can
# model the gains on a scan-by-scan basis. To form the instrument model, we need our
#   1. Our (log) gain amplitudes and phases are given below by `lgamp` and `gphase`
#   2. Our function or cache that maps the gains from a list to the sites they impact `gcache.`
#   3. The set of [`Comrade.JonesPairs`](@ref) produced by [`jonesStokes`](@ref)
# These three ingredients then specify our instrument model. The instrument model can then be
# combined with our image model `cimg` to form the total `JonesModel`.




# Now, let's set up our image model. The EHT's nominal resolution is 20-25 μas. Additionally,
# the EHT is not very sensitive to a larger field of view. Typically 60-80 μas is enough to
# describe the compact flux of M87. Given this, we only need to use a small number of pixels
# to describe our image.
npix = 32
fovx = μas2rad(150.0)
fovy = μas2rad(150.0)

# Now let's form our cache's. First, we have our usual image cache which is needed to numerically
# compute the visibilities.
grid = imagepixels(fovx, fovy, npix, npix)

# Now we need to specify our image prior. For this work we will use a Gaussian Markov
# Random field prior
# Since we are using a Gaussian Markov random field prior we need to first specify our `mean`
# image. This behaves somewhat similary to a entropy regularizer in that it will
# start with an initial guess for the image structure. For this tutorial we will use a
# a symmetric Gaussian with a FWHM of 50 μas
fwhmfac = 2*sqrt(2*log(2))
mpr = modify(Gaussian(), Stretch(μas2rad(40.0)./fwhmfac))
imgpr = intensitymap(mpr, grid)

# Now since we are actually modeling our image on the simplex we need to ensure that
# our mean image has unit flux
imgpr ./= flux(imgpr)
# and since our prior is not on the simplex we need to convert it to `unconstrained or real space`.
meanpr = to_real(CenteredLR(), Comrade.baseimage(imgpr))

# Now we can form our metadata we need to fully define our model.
# We will also fix the total flux to be the observed value 1.1. This is because
# total flux is degenerate with a global shift in the gain amplitudes making the problem
# degenerate. To fix this we use the observed total flux as our value.
skymeta = (;ftot = 1.1, meanpr, grid)



# In addition we want a reasonable guess for what the resolution of our image should be.
# For radio astronomy this is given by roughly the longest baseline in the image. To put this
# into pixel space we then divide by the pixel size.
beam = beamsize(dvis)
rat = (beam/(step(grid.X)))

# To make the Gaussian Markov random field efficient we first precompute a bunch of quantities
# that allow us to scale things linearly with the number of image pixels. The returns a
# functional that accepts a single argument related to the correlation length of the field.
# The second argument defines the underlying random field of the Markov process. Here
# we are using a zero mean and unit variance Gaussian Markov random field. The keyword
# argument specifies the order of the Gaussian field. Currently, we recommend using order
#  - 1 which is identical to TSV variation and L₂ regularization
#  - 2 which is identical to a Matern 1 process in 2D and is really the convolution of two
#    order 1 processes
# For this tutorial we will use the first order random field
crcache = ConditionalMarkov(GMRF, grid; order=1)

# To demonstrate the prior let create a few random realizations



# Now we can finally form our image prior. For this we use a heirarchical prior where the
# inverse correlation length is given by a Half-Normal distribution whose peak is at zero and
# standard deviation is `0.1/rat` where recall `rat` is the beam size per pixel.
# For the variance of the random field we use another
# half normal prior with standard deviation 0.1. The reason we use the half-normal priors is
# to prefer "simple" structures. Gaussian Markov random fields are extremly flexible models,
# and to prevent overfitting it is common to use priors that penalize complexity. Therefore, we
# want to use priors that enforce similarity to our mean image. If the data wants more complexity
# then it will drive us away from the prior.
cprior = HierarchicalPrior(crcache, truncated(InverseGamma(2.0, -log(0.1)*rat); upper=2*npix))


# We can now form our model parameter priors. Like our other imaging examples, we use a
# Dirichlet prior for our image pixels. For the log gain amplitudes, we use the `CalPrior`
# which automatically constructs the prior for the given jones cache `gcache`.
prior = NamedDist(
         c = cprior,
         fg = Uniform(0.0, 1.0),
         σimg = Exponential(0.1),
        )

skym = SkyModel(sky, prior, grid; metadata=skymeta)

post = VLBIPosterior(skym, intmodel, dvis)
# done using the `asflat` function.
tpost = asflat(post)
ndim = dimension(tpost)

# We can now also find the dimension of our posterior or the number of parameters we are going to sample.
# !!! warning
#     This can often be different from what you would expect. This is especially true when using
#     angular variables where we often artificially increase the dimension
#     of the parameter space to make sampling easier.
#-

# To initialize our sampler we will use optimize using LBFGS
using ComradeOptimization
using OptimizationOptimJL
using Zygote
# Enzyme.API.runtimeActivity!(true)
f = OptimizationFunction(tpost, Optimization.AutoZygote())
prob = Optimization.OptimizationProblem(f, randn(rng, ndim), nothing)
ℓ = logdensityof(tpost)
sol = solve(prob, LBFGS(), maxiters=5000, g_tol=1e-1);

# Now transform back to parameter space
xopt = transform(tpost, sol.u)

# !!! warning
#     Fitting gains tends to be very difficult, meaning that optimization can take a lot longer.
#     The upside is that we usually get nicer images.
#-
# First we will evaluate our fit by plotting the residuals
using Plots
residual(post, xopt)

# These look reasonable, although there may be some minor overfitting. This could be
# improved in a few ways, but that is beyond the goal of this quick tutorial.
# Plotting the image, we see that we have a much cleaner version of the closure-only image from
# [Imaging a Black Hole using only Closure Quantities](@ref).
import CairoMakie as CM
CM.activate!(type = "png", px_per_unit=3) #hide
g = imagepixels(fovx, fovy, 128, 128)
img = intensitymap(skymodel(post, xopt), g)
imageviz(img, size=(500, 400))


# Because we also fit the instrument model, we can inspect their parameters.
# To do this, `Comrade` provides a `caltable` function that converts the flattened gain parameters
# to a tabular format based on the time and its segmentation.
gt = Comrade.caltable(xopt.instrument.gp)
plot(gt, layout=(3,3), size=(600,500))

# The gain phases are pretty random, although much of this is due to us picking a random
# reference sites for each scan.

# Moving onto the gain amplitudes, we see that most of the gain variation is within 10% as expected
# except LMT, which has massive variations.
gt = Comrade.caltable(exp.(xopt.instrument.lg))
plot(gt, layout=(3,3), size=(600,500))


# To sample from the posterior, we will use HMC, specifically the NUTS algorithm. For
# information about NUTS,
# see Michael Betancourt's [notes](https://arxiv.org/abs/1701.02434).
# !!! note
#     For our `metric,` we use a diagonal matrix due to easier tuning
#-
# However, due to the need to sample a large number of gain parameters, constructing the posterior
# is rather time-consuming. Therefore, for this tutorial, we will only do a quick preliminary
# run, and any posterior
# inferences should be appropriately skeptical.
#-
using ComradeAHMC
metric = DiagEuclideanMetric(ndim)
chain = sample(rng, post, AHMC(;metric), 2500; n_adapts=1500, progress=true)
#-
# !!! note
#     The above sampler will store the samples in memory, i.e. RAM. For large models this
#     can lead to out-of-memory issues. To fix that you can include the keyword argument
#     `saveto = DiskStore()` which periodically saves the samples to disk limiting memory
#     useage. You can load the chain using `load_table(diskout)` where `diskout` is
#     the object returned from sample. For more information please see [ComradeAHMC](@ref).
#-



# Now we prune the adaptation phase
chain = chain[1501:end]

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
gmeas_am = Measurements.measurement.(mchain.instrument.lg, schain.instrument.lg)
ctable_am = caltable(exp.(gmeas_am)) # caltable expects gmeas_am to be a Vector
gmeas_ph = Measurements.measurement.(mchain.instrument.gp, schain.instrument.gp)
ctable_ph = caltable(gmeas_ph)

# Now let's plot the phase curves
plot(ctable_ph, layout=(3,3), size=(600,500))
#-
# and now the amplitude curves
plot(ctable_am, layout=(3,3), size=(600,500))

# Finally let's construct some representative image reconstructions.
samples = skymodel.(Ref(post), chain[begin:20:end])
imgs = intensitymap.(samples, Ref(g))

mimg = mean(imgs)
simg = std(imgs)
fig = CM.Figure(;resolution=(400, 400))
CM.image(fig[1,1], mimg,
                   axis=(xreversed=true, aspect=1, title="Mean Image"),
                   colormap=:afmhot)
CM.image(fig[1,2], simg./mimg,
                   axis=(xreversed=true, aspect=1, title="1/SNR",),
                   colormap=:afmhot)
CM.image(fig[2,1], imgs[1],
                   axis=(xreversed=true, aspect=1,title="Draw 1"),
                   colormap=:afmhot)
CM.image(fig[2,2], imgs[end],
                   axis=(xreversed=true, aspect=1,title="Draw 2"),
                   colormap=:afmhot)
CM.hidedecorations!.(fig.content)
fig



# And viola, you have just finished making a preliminary image and instrument model reconstruction.
# In reality, you should run the `sample` step for many more MCMC steps to get a reliable estimate
# for the reconstructed image and instrument model parameters.
