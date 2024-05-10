# # Imaging a Black Hole using only Closure Quantities

# In this tutorial, we will create a preliminary reconstruction of the 2017 M87 data on April 6
# using closure-only imaging. This tutorial is a general introduction to closure-only imaging in Comrade.
# For an introduction to simultaneous
# image and instrument modeling, see [Stokes I Simultaneous Image and Instrument Modeling](@ref)


# ## Introduction to Closure Imaging
# The EHT is the highest-resolution telescope ever created. Its resolution is equivalent
# to roughly tracking a hockey puck on the moon when viewing it from the earth. However,
# the EHT is also a unique interferometer. For one, the data it produces is incredibly sparse.
# The array is formed from only eight geographic locations around the planet, each with its unique
# telescope. Additionally, the EHT observes at a much higher frequency than typical interferometers.
# As a result, it is often difficult to directly provide calibrated data since the source model can be complicated. This implies there can be large instrumental effects
# often called *gains* that can corrupt our signal. One way to deal with this is to fit quantities
# that are independent of gains. These are often called **closure quantities**. The types of
# closure quantities are briefly described in [Introduction to the VLBI Imaging Problem](@ref).
#
# In this tutorial, we will do closure-only modeling of M87 to produce preliminary images of M87.

import Pkg #hide
__DIR = @__DIR__ #hide
pkg_io = open(joinpath(__DIR, "pkg.log"), "w") #hide
Pkg.activate(__DIR; io=pkg_io) #hide
Pkg.develop(; path=joinpath(__DIR, "..", ".."), io=pkg_io) #hide
Pkg.instantiate(; io=pkg_io) #hide
Pkg.precompile(; io=pkg_io) #hide
close(pkg_io) #hide
ENV["GKSwstype"] = "nul"; #hide



# To get started, we will load Comrade
using Comrade

# Pyehtim loads eht-imaging using PythonCall this is necessary to load uvfits files
# currently.
using Pyehtim

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(123)


# ## Load the Data
# To download the data visit https://doi.org/10.25739/g85n-f134
# To load the eht-imaging obsdata object we do:
obs = ehtim.obsdata.load_uvfits(joinpath(__DIR, "../Data/SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))

# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      are coherent.
#   - Add 2% systematic noise to deal with calibration issues such as leakage.
obs = scan_average(obs).add_fractional_noise(0.02)

# Now, we extract our closure quantities from the EHT data set.
dlcamp, dcphase  = extract_table(obs, LogClosureAmplitudes(;snrcut=3), ClosurePhases(;snrcut=3))

# ## Build the Model/Posterior
# For our model, we will be using an image model that consists of a raster of point sources,
# convolved with some pulse or kernel to make a `ContinuousImage` object with it `Comrade's.`
# generic image model. For this we will define a two argument function. The first argument
# is typically a named tuple with the model parameters. The second argument defines the metadata
# for the model that is typically constant. This allows us to explicitly pass arguments that
# are constant to the model, such as the image `cache` object that we will define below.
function sky(θ, metadata)
    (;fg, c, σimg) = θ
    (;meanpr, grid) = metadata
    ## Construct the image model we fix the flux to 0.6 Jy in this case
    cp = meanpr .+ σimg.*c.params
    rast = ((1-fg))*to_simplex(CenteredLR(), cp)
    m = ContinuousImage(rast, grid, BSplinePulse{3}())
    ## Add a large-scale gaussian to deal with the over-resolved mas flux
    g = modify(Gaussian(), Stretch(μas2rad(250.0), μas2rad(250.0)), Renormalize(fg))
    x0, y0 = centroid(m.img)
    return shifted(m, -x0, -y0) + g
end


# Now, let's set up our image model. The EHT's nominal resolution is 20-25 μas. Additionally,
# the EHT is not very sensitive to a larger field of views; typically, 60-80 μas is enough to
# describe the compact flux of M87. Given this, we only need to use a small number of pixels
# to describe our image.
npix = 32
fovxy = μas2rad(150.0)

# To define the image model we need to specify both the grid we will be using and the
# FT algorithm we will use, in this case the NFFT which is the most efficient.
grid = imagepixels(fovxy, fovxy, npix, npix)



# Now we need to specify our image prior. For this work we will use a Gaussian Markov
# Random field prior
using VLBIImagePriors, Distributions, DistributionsAD
# Since we are using a Gaussian Markov random field prior we need to first specify our `mean`
# image. For this work we will use a symmetric Gaussian with a FWHM of 50 μas
fwhmfac = 2*sqrt(2*log(2))
mpr = modify(Gaussian(), Stretch(μas2rad(40.0)./fwhmfac))
imgpr = intensitymap(mpr, grid)

# Now since we are actually modeling our image on the simplex we need to ensure that
# our mean image has unit flux before we transform
imgpr ./= flux(imgpr)
meanpr = to_real(CenteredLR(), baseimage(imgpr));

#
skymeta = (;meanpr, grid);

# In addition we want a reasonable guess for what the resolution of our image should be.
# For radio astronomy this is given by roughly the longest baseline in the image. To put this
# into pixel space we then divide by the pixel size.
beam = beamsize(dlcamp)
rat = (beam/(step(grid.X)))

# To make the Gaussian Markov random field efficient we first precompute a bunch of quantities
# that allow us to scale things linearly with the number of image pixels. This drastically improves
# the usual N^3 scaling you get from usual Gaussian Processes.
crcache = ConditionalMarkov(GMRF, grid; order=1)

# Now we can finally form our image prior. For this we use a heirarchical prior where the
# correlation length is given by a inverse gamma prior to prevent overfitting.
# Gaussian Markov random fields are extremly flexible models.
# To prevent overfitting it is common to use priors that penalize complexity. Therefore, we
# want to use priors that enforce similarity to our mean image, and prefer smoothness.
cprior = HierarchicalPrior(crcache, truncated(InverseGamma(1.0, -log(0.1)*rat); upper=2*npix))

prior = (c = cprior, σimg = Exponential(0.5), fg=Uniform(0.0, 1.0))

# Form the sky model
skym = SkyModel(sky, prior, grid; metadata=skymeta)
post = VLBIPosterior(skym, dlcamp, dcphase)

# ## Reconstructing the Image

# To sample from this posterior, it is convenient to first move from our constrained parameter space
# to an unconstrained one (i.e., the support of the transformed posterior is (-∞, ∞)). This is
# done using the `asflat` function.
tpost = asflat(post)

# We can now also find the dimension of our posterior or the number of parameters we will sample.
# !!! warning
#     This can often be different from what you would expect. This is especially true when using
#     angular variables, where we often artificially increase the dimension
#     of the parameter space to make sampling easier.
#-
ndim = dimension(tpost)


# Now we optimize using LBFGS
using ComradeOptimization
using OptimizationOptimJL
using Zygote
f = OptimizationFunction(tpost, Optimization.AutoZygote())
prob = Optimization.OptimizationProblem(f, prior_sample(rng, tpost), nothing)
sol = solve(prob, LBFGS(); maxiters=1000);


# Before we analyze our solution we first need to transform back to parameter space.
xopt = transform(tpost, sol)

# First we will evaluate our fit by plotting the residuals
using DisplayAs #hide
using Plots
p = residual(post, xopt)
DisplayAs.Text(DisplayAs.PNG(p)) #hide

# Now let's plot the MAP estimate.
import CairoMakie as CM
CM.activate!(type = "png", px_per_unit=1) #hide
g = imagepixels(μas2rad(150.0), μas2rad(150.0), 100, 100)
img = intensitymap(skymodel(post, xopt), g)
fig = imageviz(img);
DisplayAs.Text(DisplayAs.PNG(fig)) #hide


# That doesn't look great. The reason for this is that the MAP, from a Bayesian standpoint,
# is not representative of the posterior. In fact, the MAP often overfits the data and depends
# on the choice of parameterization of the posterior. Given this we will turn to sampling
# from the posterior instead.

# To sample from the posterior we will use HMC and more specifically the NUTS algorithm. For information about NUTS
# see Michael Betancourt's [notes](https://arxiv.org/abs/1701.02434).
# !!! note
#     For our `metric` we use a diagonal matrix due to easier tuning.
#-
using ComradeAHMC
using Zygote
metric = DiagEuclideanMetric(ndim)
chain = sample(post, AHMC(;metric, autodiff=Val(:Zygote)), 700; n_adapts=500, progress=false, initial_params=xopt);


# !!! warning
#     This should be run for longer!
#-
# Now that we have our posterior, we can assess which parts of the image are strongly inferred by the
# data. This is rather unique to `Comrade` where more traditional imaging algorithms like CLEAN and RML are inherently
# unable to assess uncertainty in their reconstructions.
#
# To explore our posterior let's first create images from a bunch of draws from the posterior
msamples = skymodel.(Ref(post), chain[501:2:end]);

# The mean image is then given by
using StatsBase
imgs = intensitymap.(msamples, Ref(g))
mimg = mean(imgs)
simg = std(imgs)
fig = CM.Figure(;resolution=(700, 700));
CM.image(fig[1,1], mimg,
                   axis=(xreversed=true, aspect=1, title="Mean Image"),
                   colormap=:afmhot);
CM.image(fig[1,2], simg./(max.(mimg, 1e-8)),
                   axis=(xreversed=true, aspect=1, title="1/SNR",), colorrange=(0.0, 2.0),
                   colormap=:afmhot);
CM.image(fig[2,1], imgs[1],
                   axis=(xreversed=true, aspect=1,title="Draw 1"),
                   colormap=:afmhot);
CM.image(fig[2,2], imgs[end],
                   axis=(xreversed=true, aspect=1,title="Draw 2"),
                   colormap=:afmhot);
CM.hidedecorations!.(fig.content)
DisplayAs.Text(DisplayAs.PNG(fig)) #hide

# Now let's see whether our residuals look better.
p = Plots.plot(layout=(2,1))
for s in sample(chain[501:end], 10)
    residual!(post, s)
end
p


# And viola, you have a quick and preliminary image of M87 fitting only closure products.
# For a publication-level version we would recommend
#    1. Running the chain longer and multiple times to properly assess things like ESS and R̂ (see [Geometric Modeling of EHT Data](@ref))
#    2. Fitting gains. Typically gain amplitudes are good to 10-20% for the EHT not the infinite uncertainty closures implicitly assume
