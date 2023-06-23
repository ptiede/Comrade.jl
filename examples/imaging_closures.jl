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


# To get started, we will load Comrade
using Comrade


using Pkg #hide
Pkg.activate(joinpath(dirname(pathof(Comrade)), "..", "examples")) #hide

using Pyehtim

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(123)


# ## Load the Data
# To download the data visit https://doi.org/10.25739/g85n-f134
# To load the eht-imaging obsdata object we do:
obs = ehtim.obsdata.load_uvfits(joinpath(dirname(pathof(Comrade)), "..", "examples", "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))

# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      are coherent.
#   - Add 1% systematic noise to deal with calibration issues that cause 1% non-closing errors.
obs = scan_average(obs).add_fractional_noise(0.01).flag_uvdist(uv_min=0.1e9)

# Now, we extract our closure quantities from the EHT data set.
dlcamp, dcphase  = extract_table(obs, LogClosureAmplitudes(;snrcut=3), ClosurePhases(;snrcut=3))

# ## Build the Model/Posterior
# For our model, we will be using an image model that consists of a raster of point sources,
# convolved with some pulse or kernel to make a `ContinuousImage` object with it `Comrade's.`
# generic image model. Note that `ContinuousImage(img, cache)` actually creates a [`Comrade.ModelImage`](@ref)
# object that allows `Comrade` to numerically compute the Fourier transform of the image.
function model(θ, metadata)
    (;c) = θ
    (;grid, cache) = metadata
    ## Construct the image model
    c = to_simplex(CenteredLR(), c.params)
    img = IntensityMap(c, grid)
    return  ContinuousImage(img, cache)
end


# Now, let's set up our image model. The EHT's nominal resolution is 20-25 μas. Additionally,
# the EHT is not very sensitive to a larger field of views; typically, 60-80 μas is enough to
# describe the compact flux of M87. Given this, we only need to use a small number of pixels
# to describe our image.
npix = 24
fovxy = μas2rad(100.0)

# Now, we can feed in the array information to form the cache
grid = imagepixels(fovxy, fovxy, npix, npix)
buffer = IntensityMap(zeros(npix,npix), grid)
cache = create_cache(NFFTAlg(dlcamp), buffer, BSplinePulse{3}())
metadata = (;K, grid, cache)



# Now we need to specify our image prior. For this work we will use a Gaussian Markov
# Random field prior
using VLBIImagePriors, Distributions, DistributionsAD
# Since we are using a Gaussian Markov random field prior we need to first specify our `mean`
# image. For this work we will use a symmetric Gaussian with a FWHM of 40 μas
fwhmfac = 2*sqrt(2*log(2))
mpr = modify(Gaussian(), Stretch(μas2rad(40.0)./fwhmfac))
imgpr = intensitymap(mpr, grid)

# Now since we are actually modeling our image on the simplex we need to ensure that
# our mean image has unit flux
imgpr ./= flux(imgpr)

meanpr = to_real(CenteredLR(), Comrade.baseimage(imgpr))

# In addition we want a reasonable guess for what the resolution of our image should be.
# For radio astronomy this is given by roughly the longest baseline in the image. To put this
# into pixel space we then divide by the pixel size.
hh(x) = hypot(x...)
beam = inv(maximum(hh.(uvpositions.(extract_table(obs, ComplexVisibilities()).data))))
rat = (beam/(step(grid.X)))

# To make the Gaussian Markov random field efficient we first precompute a bunch of quantities
# that allow us to scale things linearly with the number of image pixels. This drastically improves
# the usual N^3 scaling you get from usual Gaussian Processes.
crcache = MarkovRandomFieldCache(meanpr)

# One of the benefits of the Bayesian approach is that we can fit for the hyperparameters
# of our prior/regularizers unlike traditional RML appraoches. To construct this heirarchical
# prior we will first make a map that takes in our regularizer hyperparameters and returns
# the image prior given those hyperparameters.
fmap = let meanpr=meanpr, crcache=crcache, rat=rat
    x->GaussMarkovRandomField(meanpr, x.λ /rat, x.σ^2, crcache)
end

# Now we can finally form our image prior. For this we use a heirarchical prior where the
# inverse correlation length is given by a Half-Normal distribution whose peak is at zero and
# standard deviation is 1/3/rat. For the variance of the GP we use another
# half normal prior with standard deviation unity. The reason we use the half-normal priors is
# to prefer "simple" structures. Namely, Gaussian Markov random fields are extremly flexible models.
# To prevent overfitting it is common to use priors that penalize complexity. Therefore, we
# want to use priors that enforce similarity to our mean image, and prefer smoothness.
cprior = HierarchicalPrior(fmap, Comrade.NamedDist((;λ = truncated(Normal(0.0, 0.1); lower=0.0), σ=truncated(Normal(0.0, 1.0); lower=0.0))))

prior = (c = cprior, )

lklhd = RadioLikelihood(model, dlcamp, dcphase;
                        skymeta = metadata)
post = Posterior(lklhd, prior)

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
sol = solve(prob, LBFGS(); maxiters=5_000);


# Before we analyze our solution we first need to transform back to parameter space.
xopt = transform(tpost, sol)

# First we will evaluate our fit by plotting the residuals
using Plots
residual(skymodel(post, xopt), dlcamp, ylabel="Log Closure Amplitude Res.")
# and now closure phases
#-
residual(skymodel(post, xopt), dcphase, ylabel="|Closure Phase Res.|")

# Now these residuals look a bit high. However, it turns out this is because the MAP is typically
# not a great estimator and will not provide very predictive measurements of the data. We
# will show this below after sampling from the posterior.
img = intensitymap(skymodel(post, xopt), μas2rad(100.0), μas2rad(100.0), 100, 100)
plot(img, title="MAP Image")

# To sample from the posterior we will use HMC and more specifically the NUTS algorithm. For information about NUTS
# see Michael Betancourt's [notes](https://arxiv.org/abs/1701.02434).
# !!! note
#     For our `metric` we use a diagonal matrix due to easier tuning.
#-
using ComradeAHMC
using Zygote
metric = DiagEuclideanMetric(ndim)
chain, stats = sample(post, AHMC(;metric, autodiff=Val(:Zygote)), 2_000; nadapts=1_000, init_params=xopt)


# !!! warning
#     This should be run for longer!
#-
# Now that we have our posterior, we can assess which parts of the image are strongly inferred by the
# data. This is rather unique to `Comrade` where more traditional imaging algorithms like CLEAN and RML are inherently
# unable to assess uncertainty in their reconstructions.
#
# To explore our posterior let's first create images from a bunch of draws from the posterior
msamples = skymodel.(Ref(post), chain[1001:10:end]);

# The mean image is then given by
using StatsBase
imgs = intensitymap.(msamples, μas2rad(100.0), μas2rad(100.0), 128, 128)
mimg = mean(imgs)
simg = std(imgs)
p1 = plot(mimg, title="Mean Image")
p2 = plot(simg./mimg, title="1/SNR")
p3 = plot(imgs[1], title="Draw 1")
p4 = plot(imgs[end], title="Draw 2")
plot(p1, p2, p3, p4, size=(800,800), colorbar=:none)

# Now let's see whether our residuals look better.
p = plot();
for s in sample(chain[1001:end], 10)
    residual!(p, vlbimodel(post, s), dlcamp)
end
ylabel!("Log-Closure Amplitude Res.");
p
#-


p = plot();
for s in sample(chain[1001:end], 10)
    residual!(p, vlbimodel(post, s), dcphase)
end
ylabel!("Closure Phase Res.");
p

# And we see that the residuals are looking much better.



# And viola, you have a quick and preliminary image of M87 fitting only closure products.
# For a publication-level version we would recommend
#    1. Running the chain longer and multiple times to properly assess things like ESS and R̂ (see [Geometric Modeling of EHT Data](@ref))
#    2. Fitting gains. Typically gain amplitudes are good to 10-20% for the EHT not the infinite uncertainty closures implicitly assume
#    3. Making sure the posterior is unimodal (hint for this example it isn't!). The EHT image posteriors can be pretty complicated, so typically
#       you want to use a sampler that can deal with multi-modal posteriors. Check out the package [`Pigeons.jl`](https://github.com/Julia-Tempering/Pigeons.jl)
#       for an **in-development** package that should easily enable this type of sampling.
#-


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
