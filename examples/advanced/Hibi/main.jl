import Pkg #hide
__DIR = @__DIR__ #hide
pkg_io = open(joinpath(__DIR, "pkg.log"), "w") #hide
Pkg.activate(__DIR; io = pkg_io) #hide
Pkg.develop(; path = joinpath(__DIR, "..", "..", ".."), io = pkg_io) #hide
Pkg.instantiate(; io = pkg_io) #hide
Pkg.precompile(; io = pkg_io) #hide
close(pkg_io) #hide


# # Hierarchicial Interferometric Bayesian Imaging (HIBI)

# In this tutorial, we will demonstrate how to use `Comrade` to utilize
# Bayesian hierarchical modeling to reconstruct an image of M87.
# Regular imaging tries to make minimal assumptions about the source structure.
# We will demonstrate how a use can incorporate their knowledge of the source
# structure into the model to potentially improve the image reconstruction.

# Explicitly, we will take advantage of the fact that optically thin accretion flows
# from a black hole that is face on is expected to produce a ring-like structure.
# To take advantage of this we will use a hierarchical imaging approach where we
# can easily incorporate this domain model into the image reconstruction.

# This approach is called Hierarchical Interferometric Bayesian Imaging (HIBI) and is
# described in detail in the paper Hierarchical Interferometric Bayesian Imaging (in prep).
# For our image model we will use a raster image similar to
# [Stokes I Simultaneous Image and Instrument Modeling](@ref). We will decompose the image raster
# components $F_{ij}$ as follows
# ```math
# F_{ij} = F_0 \frac{\mu(x_i,y_j) \exp(σ r_{ij})}{\sum_{ij} \mu(x_i,y_j) \exp(σ r_{ij})}
# ```
# where $F_0$ is the total flux of the image, $μ(x_i,y_j)$ is a mean image and $r_{ij}$ is a
# is the multiplicative fluctuations about this mean image.

# However, before we continue let's load the packages and tutorials we will need.


# ## Loading the Data

# To get started we will load Comrade
using Comrade

# ## Load the Data


using Pyehtim

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(42)


# To download the data visit https://doi.org/10.25739/g85n-f134
# To load the eht-imaging obsdata object we do:
obs = ehtim.obsdata.load_uvfits(joinpath(__DIR, "..", "..", "Data", "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))

# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      coherent.
obs = scan_average(obs)

# For this tutorial we will once again fit complex visibilities since they
# provide the most information once the telescope/instrument model are taken
# into account. Additionally, we will preprocess the data to remove effects
# we don't want to model. Specifically we will remove the
#  - Baselines shorter than 0.1Gλ since this represents overresolved structure
#  - Add 1% systematic uncertainty to handle residual calibration errors such as leakage.
dvis = add_fractional_noise(flag(x -> uvdist(x) < 0.1e9, extract_table(obs, Visibilities())), 0.01)

# ## Building the Model/Posterior

# Now let's construct our model using the decomposition described above. For this we will need
# to define

function sky(θ, metadata)
    (; c, σimg, r, ain, aout) = θ
    (; ftot, grid) = metadata
    ## Form the image model
    ## First transform to simplex space first applying the non-centered transform
    mb = RingTemplate(RadialDblPower(ain, aout), AzimuthalUniform())
    mr = modify(mb, Stretch(r))
    mimg = intensitymap(mr, grid)
    rast = apply_fluctuations(CenteredLR(), mimg, σimg * c.params)
    rast .*= ftot
    mimg = ContinuousImage(rast, grid, BSplinePulse{3}())
    return mimg
end

# Unlike other imaging examples
# (e.g., [Imaging a Black Hole using only Closure Quantities](@ref)) we also need to include
# a model for the instrument, i.e., gains as well. The gains will be broken into two components
#   - Gain amplitudes which are typically known to 10-20%, except for LMT, which has amplitudes closer to 50-100%.
#   - Gain phases which are more difficult to constrain and can shift rapidly.

using VLBIImagePriors
using Distributions
fgain(x) = exp(x.lg + 1im * x.gp)
G = SingleStokesGain(fgain)

intpr = (
    lg = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.2)); LM = IIDSitePrior(ScanSeg(), Normal(0.0, 1.0))),
    gp = ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))); refant = SEFDReference(0.0)),
)
intmodel = InstrumentModel(G, intpr)


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
fovxy = μas2rad(150.0)
npix = 48
g = imagepixels(fovxy, fovxy, npix, npix)


# Part of hybrid imaging is to force a scale separation between
# the different model components to make them identifiable.
# To enforce this we will set the raster component to have a
# correlation length of 5 times the beam size.
cprior = corr_image_prior(g, dvis)


# For the other parameters we use a uniform priors for the ring fractional flux `f`
# ring radius `r`, ring width `σ`, and the flux fraction of the Gaussian component `fg`
# and the amplitude for the ring brightness modes. For the angular variables `ξτ` and `ξ`
# we use the von Mises prior with concentration parameter `inv(π^2)` which is essentially
# a uniform prior on the circle. Finally for the standard deviation of the MRF we use a
# half-normal distribution. This is to ensure that the MRF has small differences from the
# mean image.
skyprior = (
    c = cprior,
    σimg = truncated(Normal(0.0, 0.5); lower = 0.0),
    r = Uniform(μas2rad(10.0), μas2rad(40.0)),
    ain = Uniform(1.0, 20.0),
    aout = Uniform(1.0, 20.0),
)

# Now we form the metadata
skymetadata = (; ftot = 1.1, grid = g)
skym = SkyModel(sky, skyprior, g; metadata = skymetadata)

# This is everything we need to specify our posterior distribution, which our is the main
# object of interest in image reconstructions when using Bayesian inference.
using Enzyme
post = VLBIPosterior(skym, intmodel, dvis)

# We can sample from the prior to see what the model looks like
using DisplayAs #hide
using CairoMakie
xrand = prior_sample(rng, post)
gpl = refinespatial(g, 3)
fig = imageviz(intensitymap(skymodel(post, xrand), gpl));
fig |> DisplayAs.PNG |> DisplayAs.Text #hide


# ## Reconstructing the Image

# To find the image we will demonstrate two methods:
#  - Optimization to find the MAP (fast but often a poor estimator)
#  - Sampling to find the posterior (slow but provides a substantially better estimator)
# For optimization we will use the `Optimization.jl` package and the LBFGS optimizer.
# To use this we use the [`comrade_opt`](@ref) function
using Optimization, OptimizationLBFGSB
xopt, sol = comrade_opt(
    post, LBFGSB();
    initial_params = xrand, maxiters = 2000, g_tol = 1.0e0
);


# First we will evaluate our fit by plotting the residuals
res = residuals(post, xopt);
fig = plotfields(res[1], :uvdist, :res);
fig |> DisplayAs.PNG |> DisplayAs.Text #hide

# These residuals suggest that we are substantially overfitting the data. This is a common
# side effect of MAP imaging. As a result if we plot the image we see that there
# is substantial high-frequency structure in the image that isn't supported by the data.
fig = imageviz(intensitymap(skymodel(post, xopt), gpl), figure = (; resolution = (500, 400)));
fig |> DisplayAs.PNG |> DisplayAs.Text #hide


# To improve our results we will now move to Posterior sampling. This is the main method
# we recommend for all inference problems in `Comrade`. While it is slower the results are
# often substantially better. To sample we will use the `AdvancedHMC` package.
using AdvancedHMC
chain = sample(rng, post, NUTS(0.8), 700; n_adapts = 500, progress = false, initial_params = xopt);
# chain = load_samples(out)
# We then remove the adaptation/warmup phase from our chain
chain = chain[501:end]

# !!! warning
#     This should be run for 4-5x more steps to properly estimate expectations of the posterior
#-

# Now lets plot the mean image and standard deviation images.
# To do this we first clip the first 250 MCMC steps since that is during tuning and
# so the posterior is not sampling from the correct sitesary distribution.

using StatsBase
msamples = skymodel.(Ref(post), chain[begin:5:end]);

# The mean image is then given by
imgs = intensitymap.(msamples, Ref(gpl))
fig = imageviz(mean(imgs), size = (400, 300));
fig |> DisplayAs.PNG |> DisplayAs.Text #hide
#-
fig = imageviz(std(imgs), colormap = :batlow, size = (400, 300));
fig |> DisplayAs.PNG |> DisplayAs.Text #hide

#-
#


# Finally, let's take a look at some of the ring parameters

figd = Figure(; resolution = (650, 200));
p1 = density(figd[1, 1], rad2μas(chain.sky.r) * 2, axis = (xlabel = "Ring Pwr Law Transition (μas)",))
p2 = density(figd[1, 2], chain.sky.ain, axis = (xlabel = "Inner Radial Index",))
p3 = density(figd[1, 3], chain.sky.aout, axis = (xlabel = "Outer Radial Index",))
figd |> DisplayAs.PNG |> DisplayAs.Text #hide

# Now let's check the residuals using draws from the posterior
fig = Figure(; size = (600, 400))
resch = residuals(post, chain[end])
ax, = plotfields!(fig[1, 1], resch[1], :uvdist, :res, scatter_kwargs = (; label = "MAP", color = :blue, colorim = :red, marker = :circle), legend = false)
for s in sample(chain, 10)
    baselineplot!(ax, residuals(post, s)[1], :uvdist, :res, alpha = 0.2, label = "Draw")
end
axislegend(ax, merge = true)
fig |> DisplayAs.PNG |> DisplayAs.Text #hide
