import Pkg #hide
__DIR = @__DIR__ #hide
pkg_io = open(joinpath(__DIR, "pkg.log"), "w") #hide
Pkg.activate(__DIR; io = pkg_io) #hide
Pkg.develop(; path = joinpath(__DIR, "..", "..", ".."), io = pkg_io) #hide
Pkg.instantiate(; io = pkg_io) #hide
Pkg.precompile(; io = pkg_io) #hide
close(pkg_io) #hide

# # Modeling a Power Spectrum of AGN with Markov Random Field Expansion

# In this tutorial, we will reconstruct VLBA data observations of the AGN XYZ using closures.
# However we will also model the power spectrum of the AGN using a Markov Random Field expansion
# Where we can fit multiple scales of the power spectrum simultaneously up to some order. 


# ## Introduction to Closure Imaging
#
# In this tutorial, we will do closure-only modeling of the AGN XYZ observed with the VLBA at XYZ GHz.


# To get started, we will load Comrade
using Comrade
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)
VLBISkyModels.FFTW.set_num_threads(Threads.nthreads())
# Pyehtim loads eht-imaging using PythonCall this is necessary to load uvfits files
# currently.
using Pyehtim
using NonuniformFFTs

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(123)


# ## Load the Data
# To download the data visit https://doi.org/10.25739/g85n-f134
# To load the eht-imaging obsdata object we do:
# obs = ehtim.obsdata.load_uvfits("/home/ptiede/Smithsonian External Dropbox/Paul Tiede/MixedPolPaper/data/VLBA/BU/0954+658Q.2025-03-23.UVP.gz")
obs0 = ehtim.obsdata.load_uvfits(joinpath(__DIR, "..", "..", "Data", "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))

# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      are coherent.
#   - Add 2% systematic noise to deal with calibration issues such as leakage.
obs = scan_average(obs)

# Now, we extract our closure quantities from the EHT data set. We flag now SNR points since
# the closure likelihood we use is only applicable to high SNR data.
dlcamp, dcphase = extract_table(obs, LogClosureAmplitudes(; snrcut = 3), ClosurePhases(; snrcut = 3))

# !!! note
#     Fitting low SNR closure data is complicated and requires a more sophisticated likelihood.
#     If low-SNR data is very important we recommend fitting visibilties with a instrumental model.

function fast_centroid(img::IntensityMap{<:Real, 2})
    x0 = zero(eltype(img))
    y0 = zero(eltype(img))
    dp = domainpoints(img)
    fs = Comrade._fastsum(img)
    @inbounds for i in CartesianIndices(img)
        x0 += dp[i].X * img[i]
        y0 += dp[i].Y * img[i]
    end
    return x0 / fs, y0 / fs
end



# ## Build the Model/Posterior
# For our model, we will be using an image model that consists of a raster of point sources,
# convolved with some pulse or kernel to make a `ContinuousImage`.
# To define this model we define the standard two argument function `sky` that defines the
# sky model we want to fit. The first argument are the model parameters, and are typically
# a NamedTuple. The second argument defines the metadata
# for the model that is typically constant. For our model the constant `metdata` will just
# by the mean or prior image.
function sky(θ, metadata)
    (; fb, c, ρs, τ, ξτ, σimg) = θ
    (; mimg, pl) = metadata
    ## Apply the GMRF fluctuations to the image
    x = genfield(StationaryRandomField(MarkovPS(ρs), pl), c)
    x .= σimg .* x
    fbn = fb/length(mimg)
    mb = mimg.*(1 - fb) .+ fbn
    rast = apply_fluctuations(UnitFluxMap(exp), mb, x)
    x0, y0 = fast_centroid(rast)
    m = ContinuousImage(rast, DeltaPulse())
    ## Force the image centroid to be at the origin
    ## Add a large-scale gaussian to deal with the over-resolved mas flux
    return shifted(m, -x0, -y0)
end


# Now, let's set up our image model. The EHT's nominal resolution is 20-25 μas. Additionally,
# the EHT is not very sensitive to a larger field of views; typically, 60-80 μas is enough to
# describe the compact flux of M87. Given this, we only need to use a small number of pixels
# to describe our image.
npix = 64
fovxy = μas2rad(2000.0)

# To define the image model we need to specify both the grid we will be using and the
# FT algorithm we will use, in this case the NFFT which is the most efficient.
grid = imagepixels(fovxy, fovxy, npix, npix)


# Now we need to specify our image prior. For this work we will use a Gaussian Markov
# Random field prior
using VLBIImagePriors, Distributions

# Since we are using a Gaussian Markov random field prior we need to first specify our `mean`
# image. For this work we will use a symmetric Gaussian with a FWHM of 50 μas
fwhmfac = 2 * sqrt(2 * log(2))
mpr = modify(Gaussian(), Stretch(3*beamsize(dlcamp) / fwhmfac))
imgpr = intensitymap(mpr, grid)
# To momdel the power spectrum we also need to construct our execution plan for the given grid.
# This will be used to construct the actual correlated realization of the RF given some initial
# white noise.
pl = StationaryRandomFieldPlan(grid)
skymeta = (; mimg=imgpr./sum(imgpr), pl);


# Now we can finally form our image prior. For this we use a heirarchical prior where the
# direct log-ratio image prior is a Gaussian Markov Random Field. The correlation length
# of the GMRF is a hyperparameter that is fit during imaging. We pass the data to the prior
# to estimate what the maximumal resolutoin of the array is and prevent the prior from allowing
# correlation lengths that are much small than the telescope beam size. Note that this GMRF prior
# has unit variance. For more information on the GMRF prior see the [`corr_image_prior`](@ref) doc string.
cprior = std_dist(pl)

# Putting everything together the total prior is then our image prior, a prior on the
# standard deviation of the MRF, and a prior on the fractional flux of the Gaussian component.
dρ = Uniform(0.1, 4 * max(size(grid)...))
prior = (
    c = cprior,
    ρs = ntuple(Returns(dρ), 2),
    σimg = truncated(Normal(0.0, 0.5); lower=0.0),
    fb = Uniform(0.0, 1.0),
)

# We can then define our sky model.
skym = SkyModel(sky, prior, grid; metadata = skymeta, algorithm = NonuniformFFTsAlg())

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
# In this tutorial we will use Optiizations LBFGS optimizer.
# We also need to import Enzyme to allow for automatic differentiation.
using Optimization
xopt, sol = comrade_opt(
    post, Optimization.LBFGS();
    maxiters = 2000, initial_params = prior_sample(rng, post)
);

using CairoMakie
using DisplayAs #hide
g = refinespatial(grid, 2)
# img = intensitymap(skymodel(post, xopt), g)
img = skymodel(post, xopt).model.img
fig = imageviz(img, size = (600, 500));
DisplayAs.Text(DisplayAs.PNG(fig)) #hide


# First we will evaluate our fit by plotting the residuals
res = Comrade.residuals(post, xopt)
fig = Figure(; size = (800, 300))
plotfields!(fig[1, 1], res[1], :uvdist, :res);
plotfields!(fig[1, 2], res[2], :uvdist, :res);
fig |> DisplayAs.PNG |> DisplayAs.Text

# Now let's plot the MAP estimate.


# That doesn't look great. This is pretty common for the sparse EHT data. In this case the
# MAP often drastically overfits the data, producing a image filled with artifacts. In addition,
# we note that the MAP itself is not invariant to the model parameterization. Namely, if we
# changed our prior to use a fully centered parameterization we would get a very different image.
# Fortunately, these issues go away when we sample from the posterior, and construct expectations
# of the posterior, like the mean image.


# To sample from the posterior we will use HMC and more specifically the NUTS algorithm. For information about NUTS
# see Michael Betancourt's [notes](https://arxiv.org/abs/1701.02434).
# !!! note
#     For our `metric` we use a diagonal matrix due to easier tuning.
#-
using AdvancedHMC
out = sample(rng, post, AdvancedHMC.NUTS(0.8), 1000 + 1000, n_adapts = 1000, 
            saveto = DiskStore(name = joinpath(@__DIR__, "gausstest"), stride = 10), 
            initial_params = xopt, restart=false);
chain = load_samples(joinpath(@__DIR__, "gausstest"), 1:69*10)
# !!! warning
#     This should be run for longer!
#-
# Now that we have our posterior, we can assess which parts of the image are strongly inferred by the
# data. This is rather unique to `Comrade` where more traditional imaging algorithms like CLEAN and RML are inherently
# unable to assess uncertainty in their reconstructions.
#
# To explore our posterior let's first create images from a bunch of draws from the posterior
msamples = skymodel.(Ref(post), chain[301:5:end]);

k = range(1 / size(grid)[1], π/2, length = 512)
fig = Figure()
ax = Axis(fig[1, 1], xscale = log10, yscale = log10)
for i in 301:5:length(chain)
    lines!(ax, k, VLBIImagePriors.ampspectrum.(Ref(MarkovPS(chain.sky.ρs[i])), tuple.(k, 0)))
end
fig

# The mean image is then given by
using StatsBase
imgs = center_image.(parent.(VLBISkyModels.unmodified.(msamples)))
mimg = mean(imgs)
simg = std(imgs)
fig = Figure(; resolution = (700, 400));
axs = [Axis(fig[i, j], xreversed = true, aspect = DataAspect()) for i in 1:2, j in 1:2]
image!(axs[1, 1], mimg, colormap = :afmhot, ); axs[1, 1].title = "Mean"
image!(axs[1, 2], simg ./ (max.(mimg, 1.0e-8)), colormap = :afmhot);axs[1, 2].title = "Std"
image!(axs[2, 1], sample(imgs), colormap = :afmhot, );
image!(axs[2, 2], sample(imgs), colormap = :afmhot, );
hidedecorations!.(axs)
fig |> DisplayAs.PNG |> DisplayAs.Text

gpl = imagepixels(μas2rad(100.0), μas2rad(100.0), 128, 128)
pimgs = regrid.(imgs, Ref(gpl))

fig = Figure(;resolution=(600, 400))
axs = [Axis(fig[i, j], xreversed = true, aspect = DataAspect()) for i in 1:2, j in 1:3]
map(enumerate(axs)) do (i, ax)
    hidedecorations!(ax)
    image!(ax, pimgs[i], colormap = :afmhot)
    text!(ax, 0.05, 0.9, text="χ²= $(round(mean(chi2(post, chain[300:5:end][i]; reduce=true)); digits=2))", space=:relative, color=:white)
end
fig


# Now let's see whether our residuals look better.
fig = Figure(; size = (800, 300))
res = Comrade.residuals(post, chain[end])
ax1, = baselineplot(fig[1, 1], res[1], :uvdist, :res, label = "MAP residuals", axis = (ylabel = "LCA Normalized Residuals", xlabel = "uvdist (Gλ)"))
ax2, = baselineplot(fig[1, 2], res[2], :uvdist, :res, label = "MAP residuals", axis = (ylabel = "CP Normalized Residuals", xlabel = "uvdist (Gλ)"))
ax1.title = "χ²ᵣ = $(chi2(post, chain[end]; reduce=true)[1])"
ax2.title = "χ²ᵣ = $(chi2(post, chain[end]; reduce=true)[2])"
for s in sample(chain[201:end], 10)
    rs = Comrade.residuals(post, s)
    baselineplot!(ax1, rs[1], :uvdist, :res, color = :grey, alpha = 0.2, label = "Posterior Draw")
    baselineplot!(ax2, rs[2], :uvdist, :res, color = :grey, alpha = 0.2, label = "Posterior Draw")
end
axislegend(ax1, merge = true)
fig |> DisplayAs.PNG |> DisplayAs.Text


# And viola, you have a quick and preliminary image of M87 fitting only closure products.
# For a publication-level version we would recommend
#    1. Running the chain longer and multiple times to properly assess things like ESS and R̂ (see [Geometric Modeling of EHT Data](@ref))
#    2. Fitting gains. Typically gain amplitudes are good to 10-20% for the EHT not the infinite uncertainty closures implicitly assume
