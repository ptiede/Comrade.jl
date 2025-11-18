import Pkg #hide
__DIR = @__DIR__ #hide
pkg_io = open(joinpath(__DIR, "pkg.log"), "w") #hide
Pkg.activate(__DIR; io = pkg_io) #hide
Pkg.develop(; path = joinpath(__DIR, "..", "..", ".."), io = pkg_io) #hide
Pkg.instantiate(; io = pkg_io) #hide
Pkg.precompile(; io = pkg_io) #hide
close(pkg_io) #hide

# # Modeling the Power Spectrum of an AGN with Markov Random Field Expansion

# In this tutorial, we will reconstruct VLBA data observations of the AGN XYZ using closures.
# However we will also model the power spectrum of the AGN using a Markov Random Field expansion
# Where we can fit multiple scales of the power spectrum simultaneously up to some order.


# In this tutorial, we will do closure-only modeling of the AGN DA 193 observed with the VLBA
# at 15 GHz with the Mojave AGN project. Unlike the previous tutorials, we will constrain the
# power spectrum slope and using a Markov random field expansion. This will allow us to model
# more complex and multi-scale processes in AGN, which is expected to be common in black hole
# jets.


# To get started, we will load Comrade
using Comrade
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)
# Pyehtim loads eht-imaging using PythonCall this is necessary to load uvfits files
# currently.
using Pyehtim
using NonuniformFFTs

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(123)


# ## Load the Data
# For this tutorial we will image publicly available VLBA data of the AGN
# 1308+326 observed on 2021/03/19 at 43 GHz as part of the Boston University blazar monitoring program.
file = Base.download("https://www.bu.edu/blazars/VLBA_GLAST/1308/1308+326Q.2021-03-19.UVP.gz")
obs0 = ehtim.obsdata.load_uvfits(file)

# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      are coherent.
#   - Add 0.5% systematic noise to deal with calibration issues such as leakage.
obs = scan_average(obs0).add_fractional_noise(0.005)

# For this tutorial we will only use closure quantities to reconstruct the image however, polarized
# or complex visibilities can also be used with instrumental models following the other tutorials.
dlcamp, dcphase = extract_table(obs, LogClosureAmplitudes(; snrcut = 3), ClosurePhases(; snrcut = 3))


# ## Build the Model/Posterior
# Most of the model building here will look very similar to the previous [Imaging a Black Hole using only Closure Quantities](@ref)
# tutorial. However, we will be utilizing a more complex image prior. Specifically, [`VLBIImagePriors`](https://ptiede.github.io/VLBIImagePriors.jl/stable/)
# provides a basic framework for building stationary Gaussian random fields with cyclic boundary conditions.
# To define the random field we just need to define a spectral model. For this work we will use a
# Markovian spectral model. Namely, our power spectrum will be modeled as
# ```math
#   P(k) \propto \frac{\sigma}{1 + \sum_s \rho_s k^{2s}}
# ```
# where `σ` is the marginal variance of the image, `ρs` are the coefficients of the Markovian expansion,
# and `k` is the norm of the spatial wavenumber.
using VLBIImagePriors ## Defines the `MarkovPS` power spectrum model and `StationaryRandomField`
function sky(θ, metadata)
    (; fb, c, ρs, σimg) = θ
    (; mimg, pl) = metadata
    ## Apply the GMRF fluctuations to the image
    x = genfield(StationaryRandomField(MarkovPS(ρs .^ 2), pl), c)
    x .= σimg .* x
    fbn = fb / length(mimg)
    mb = mimg .* (1 - fb) .+ fbn
    rast = apply_fluctuations(CenteredLR(), mb, x)
    m = ContinuousImage(rast, BSplinePulse{3}())
    return m
end

# For this tutorial we decided to image a very compact AGN. Thus, we will use a small FOV for a 15 GHz
# observation. Namely, we will use a 5000 μas FOV with 64x64 pixels.
nx = 64
ny = 64
fovx = μas2rad(1_000)
fovy = fovx * ny / nx
grid = imagepixels(fovx, fovy, nx, ny, μas2rad(150.0), -μas2rad(150.0))

# Now we need to specify our image prior. For this work we will use a Gaussian Markov
# Random field prior

# Since we are using a Gaussian Markov random field prior we need to first specify our `mean`
# image. For this work we will use a symmetric Gaussian with a FWHM equal to the approximate
# beamsize of the array. This models the fact that we expect the AGN core to be compact.
fwhmfac = 2 * sqrt(2 * log(2))
mpr = modify(Gaussian(), Stretch(beamsize(dlcamp) / 4 / fwhmfac))
imgpr = intensitymap(mpr, grid)
# To momdel the power spectrum we also need to construct our execution plan for the given grid.
# This will be used to construct the actual correlated realization of the RF given some initial
# white noise.
pl = StationaryRandomFieldPlan(grid)
skymeta = (; mimg = imgpr ./ sum(imgpr), pl);


# For the stationary random field prior we also need to define the *noise* prior. Luckily
# VLBIImagePriors provides a helper function to do this for us.
cprior = std_dist(pl)

# For the coefficients of the spectral expansion we will use a uniform prior between 0.1 and
# 4 times the maximum dimension of the image. This prior is rather uninformative and
# allows for a wide range of power spectra. Additionally, we truncate the expansion at order 3
# for simplicity in this tutorial.
using Distributions
ρs = ntuple(Returns(Uniform(0.1, 2 * max(size(grid)...))), 3)

# Putting everything together the total prior is then our image prior, a prior on the
# standard deviation of the MRF, and a prior on the fractional flux of the Gaussian component.
#
prior = (;
    c = cprior,
    ρs = ρs,
    σimg = Exponential(2.0),
    fb = Uniform(0.0, 1.0),
)

# We can then define our sky model.
skym = SkyModel(sky, prior, grid; metadata = skymeta)

# Since we are fitting closures we do not need to include an instrument model, since
# the closure likelihood is approximately independent of gains in the high SNR limit.
using Enzyme
post = VLBIPosterior(skym, dlcamp, dcphase)

# ## Reconstructing the Image

# To reconstruct the image we will first use the MAP estimate. This is approach is basically
# a re-implentation of regularized maximum likelihood (RML) imaging. However, unlike traditional
# RML imaging we also fit the regularizer hyperparameters, thanks to our interpretation of
# as our imaging prior as a hierarchical model.

# To optimize our posterior `Comrade` provides the `comrade_opt` function. To use this
# functionality a user first needs to import `Optimization.jl` and the optimizer of choice.
# In this tutorial we will use the Adam optimizer.
# We also need to import Enzyme to allow for automatic differentiation.
using Optimization, OptimizationOptimisers
# tpost = asflat(post)
xopt, sol = comrade_opt(post, Adam(); maxiters = 5000)

using CairoMakie
using DisplayAs #hide
# The image we actually fit is a continuous object so we can easily refine the image
# to produce a higher resolution rendering.
# Here we refine the image by a factor of 3 in each dimension.
g = refinespatial(grid, 3)
# Now to produce the intensity map we just do
imgmap = intensitymap(skymodel(post, xopt), g)
fig = imageviz(imgmap, colorscale = log10, colorrange = (1.0e-8, 1.0e-4), size = (650, 500));
DisplayAs.Text(DisplayAs.PNG(fig)) #hide


# To see how well the MAP estimate fits the data we can plot the residuals.
res = Comrade.residuals(post, xopt)
fig = Figure(; size = (800, 300))
plotfields!(fig[1, 1], res[1], :uvdist, :res);
plotfields!(fig[1, 2], res[2], :uvdist, :res);
fig |> DisplayAs.PNG |> DisplayAs.Text


# Overall, the image looks reasonable. However, the MAP is not
# a robust estimator of the image morphology. For high dimensional problems the MAP is often
# not representative of the entire image posterior. For this reason Comrade's main goal is to sample
# the posterior of the image given the data.


# To sample from the posterior we will use HMC and more specifically the NUTS algorithm similar to
# the other imaging tutorials. For this tutorial we will also show how to use the `DiskStore`
# functionality that save the chain to disk to reduce memory usage.
# This is especially useful for high-dimensional imaging problems where the chain can easily
# reach multiple GBs in size. This also allows us to restart sampling from a previous chain if needed.
# by using the keyword argument `restart=true` in the `sample` function.
using AdvancedHMC
mc = sample(
    rng, post, AdvancedHMC.NUTS(0.8), 300 + 400, n_adapts = 400,
    initial_params = xopt, saveto = DiskStore(; stride = 10, name = "VLBA_2025")
);
chain = load_samples(mc)
# !!! warning
#     This should be run for longer!
#-
# Now that we have our posterior, we can assess which parts of the image are strongly inferred by the
# data. This is rather unique to `Comrade` where more traditional imaging algorithms like CLEAN and RML are inherently
# unable to assess uncertainty in their reconstructions.
#
# To explore our posterior let's first create images from a bunch of draws from the posterior
msamples = skymodel.(Ref(post), chain[501:5:end]);

k = range(1 / size(grid)[1], π / 2, length = 512)
fig = Figure()
ax = Axis(fig[1, 1], xscale = log10, yscale = log10)
for i in 501:10:length(chain)
    lines!(ax, k, VLBIImagePriors.ampspectrum.(Ref(MarkovPS(chain.sky.ρs[i] .^ 2)), tuple.(k, 0)))
end
fig

# The mean image is then given by
using StatsBase
gpl = refinespatial(grid, 3)
imgs = intensitymap.(msamples, Ref(gpl))
mimg = mean(imgs)
simg = std(imgs)
fig = Figure(; size = (500, 300));
crange = (5.0e-6, 5.0e-2)
axs = [Axis(fig[i, j], xreversed = true, aspect = DataAspect()) for i in 1:1, j in 1:2]
image!(axs[1, 1], mimg, colormap = :afmhot, colorscale = log10, colorrange = crange); axs[1, 1].title = "Mean"
image!(axs[1, 2], simg ./ (max.(mimg, 1.0e-12)), colormap = :afmhot);axs[1, 2].title = "Fractional Uncertainty"
hidedecorations!.(axs)
fig |> DisplayAs.PNG |> DisplayAs.Text

# We can also compare the Comrade reconstruction to the CLEAN reconstruction of the same data.
cleanf = Base.download("https://www.bu.edu/blazars/VLBA_GLAST/1308/1308+326Q.2021-03-19.IMAP.gz")
# By default this will load the clean components with the beam defined in the FITS header.
mcl = load_clean_components(cleanf)
# We can also choose the load the clean components with a user-defined beam.
mcl_25 = load_clean_components(cleanf, modify(Gaussian(), Stretch(beamsize(dlcamp) / 4 / fwhmfac)))

# Now we can produce the CLEAN images on the same grid as our Comrade reconstruction.
cleanimg = intensitymap(mcl, gpl)
cleanimg25 = intensitymap(mcl_25, gpl)

fig = Figure(; size = (900, 350));
axs = [Axis(fig[1, j], xreversed = true, aspect = DataAspect()) for j in 1:3]
image!(axs[1], mimg, colormap = :afmhot, colorscale = log10, colorrange = crange); axs[1].title = "Comrade Mean"
image!(axs[2], max.(cleanimg, 1.0e-20), colormap = :afmhot, colorscale = log10, colorrange = crange); axs[2].title = "CLEAN (Nominal beam)"
image!(axs[3], max.(cleanimg25, 1.0e-20), colormap = :afmhot, colorscale = log10, colorrange = crange); axs[3].title = "CLEAN (25% beam)"
hidedecorations!.(axs)
fig |> DisplayAs.PNG |> DisplayAs.Text

# From the plot you can see that the Comrade reconstruction is significantly superresolved compared
# to the CLEAN reconstruction with the nominal beam. If we use a smaller beam for CLEAN we see
# a reconstruction that is more similar to Comrade. However, unlike CLEAN Comrade automatically
# infers the effective resolution from the data itself and does not require a restoring beam.

# Additionally, Comrade allows us to fully explore the distributions of images that are consistent
# with the data. For example, we can plot a few random samples from the posterior to see the
# variety of images that are consistent with the data.
fig = Figure(; resolution = (800, 450))
axs = [Axis(fig[i, j], xreversed = true, aspect = DataAspect()) for i in 1:2, j in 1:3]
map(enumerate(axs)) do (i, ax)
    hidedecorations!(ax)
    image!(ax, sample(imgs), colormap = :afmhot, colorscale = log10, colorrange = crange)
    text!(ax, 0.05, 0.9, text = "χ²= $(round(mean(chi2(post, chain[51:5:end][i]; reduce = true)); digits = 2))", space = :relative, color = :white)
end
axcl = Axis(fig[1:2, 4], xreversed = true, aspect = DataAspect())
hidedecorations!(axcl)
image!(axcl, max.(cleanimg25, 1.0e-20), colormap = :afmhot, colorscale = log10, colorrange = crange)
axcl.title = "CLEAN (25% beam)"
Label(fig[0, 1:3], "Comrade Post. Samples", tellheight = true)
rowgap!(fig.layout, 1, 0.0)
fig

# In summary, we have demonstrated how to use Comrade to reconstruct VLBA data of an AGN
# using only closure quantities. Additionally, we have shown how to use a Markov Random Field
# expansion to model the power spectrum of the AGN. This allows us to model more complex
# structures in the AGN jet and infer the power spectrum directly from the data.
