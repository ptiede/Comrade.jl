import Pkg #hide
__DIR = @__DIR__ #hide
pkg_io = open(joinpath(__DIR, "pkg.log"), "w") #hide
Pkg.activate(__DIR; io = pkg_io) #hide
Pkg.develop(; path = joinpath(__DIR, "..", "..", ".."), io = pkg_io) #hide
Pkg.instantiate(; io = pkg_io) #hide
Pkg.precompile(; io = pkg_io) #hide
close(pkg_io) #hide

# # Neural Field Imaging with Lux in Comrade

# In this tutorial, we will do closure-only modeling of the AGN DA 193 observed with the VLBA
# at 15 GHz with the Mojave AGN project. Unlike other tutorials, we will not attempt a Bayesian 
# analysis and instead utilize a **neural field** to obtain a MAP estimate of the image.
# This is mainly due to the challenges of characterizing posterior uncertainty with neural fields.
# The goal here is not to provide a optimal analysis of this source, or neural fields, but rather
# to demonstrate how Comrade's flexible modeling framework can be used to implement neural fields
# with only minor modifications to existing code.


# To get started, we will load Comrade
using Comrade
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(16)
# Pyehtim loads eht-imaging using PythonCall this is necessary to load uvfits files
# currently.
using Pyehtim
using NonuniformFFTs
using Lux

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
# For our model we will use a neural field to model the log image fluctuations. That is our model
# will look like
#
#   δ(x) = NN(γ(x); θ)
#
# where NN is a neural network with weights and biases θ and $x$ are the spatial coordinates. γ denotes
# a Fourier feature embedding of the spatial coordinates, which we include to potentially help with
# learning high frequency structure in the image. This kind of model for VLBI is described in more detail in 
# the excellent paper [kine](@citet). 
#
# Let's first define the Fourier feature embedding function.
function fourierfeature(grid::RectiGrid; m = 4)
    nx, ny = size(grid)
    freqs = ntuple(n -> 2^(n-1), m)
    feats = zeros(Float32, 4 * m, nx, ny)
    x = range(-1f0*π, 1f0*π, length = nx+1)[begin:end-1] #clip the endpoint
    y = range(-1f0*π, 1f0*π, length = ny+1)[begin:end-1] #clip the endpoint
    for i in 1:nx, j in 1:ny, (ik, k) in enumerate(1:4:4*m)
        s, c = sincos(freqs[ik] * x[j])
        feats[k, j, i] = s
        feats[k+1, j, i] = c
        s, c = sincos(freqs[ik] * y[i])
        feats[k + 2, j, i] = s
        feats[k + 3, j, i] = c 
    end
    return reshape(feats, 4*m, nx*ny)
end

# For this tutorial we decided to image a very compact AGN. Thus, we will use a small FOV for a 15 GHz
# observation. Namely, we will use a 1000 μas FOV with 64x64 pixels.
nx = 128
ny = 128
fovx = μas2rad(1_000)
fovy = fovx * ny / nx
grid = imagepixels(fovx, fovy, nx, ny, μas2rad(150.0), -μas2rad(150.0))


# Next we define neural network using the Lux.jl package. 
using Lux

# We only use 1 Fourier feature embedding here for simplicity and because VLBI practioners typically
# prefer very smooth images. However, more Fourier features can be used to model
# finer scale structures.
ff = fourierfeature(grid; m = 2)


# For the neural field we use a very simple MLP with 3 hidden layers of 64 units each
# and a final output layer with a single unit. We use the `tanh_fast` activation for speed
# and to ensure smoothness in the image. Also be disable the bias in the last layer since
# in some sense is prefers zero mean output which is desirable for fluctuation modeling.
nnmodel = Chain(
    Dense(size(ff,1) => 128, Lux.tanh_fast),
    Dense(128 => 128, Lux.tanh_fast),
    Dense(128 => 1, use_bias=false),
)



# We can now define our sky model function. This function takes in the neural network parameters `nn`
# and hyperparameters `fb, σ` and returns a `ContinuousImage` representing the sky brightness distribution.
# This should look similar to essentially every other imaging model in the Comrade tutorials.
function sky(θ, metadata)
    (;nn, σ, fb) = θ
    (; nnmodel, mimg, ff, st) = metadata
    fbn = fb / length(mimg)
    mb = mimg .* (1 - fb) .+ fbn

    ## Compute the neural field
    δ = reshape(first(Lux.apply(nnmodel, ff, nn, st)), size(mimg))
    δ .*= σ 
    rast = apply_fluctuations(CenteredLR(), mb, δ)
    m = ContinuousImage(rast, BSplinePulse{3}())
    return m
end




# We assume that the prior is a IID standard normal on the neural network weights and biases. To set
# this we use `Lux.setup` to get the parameter layout of the network and then use `Comrade.rmap` to
# map each parameter to a standard normal prior.
using Random
ps, st = Lux.setup(Random.default_rng(), nnmodel)
nnprior = Comrade.rmap(x->VLBIImagePriors.StdNormal(size(x)), ps)

# Similar to the RF models we will use a so-called `mean image`. Note that in this case is isn't
# actually clean what the mean image is, but this will help with centroid drift so we include it. 
# We will use a symmetric Gaussian with a FWHM equal to the approximate
# beamsize of the array. This models the fact that we expect the AGN core to be compact.
fwhmfac = 2 * sqrt(2 * log(2))
mpr = modify(TBlob(4.0), Stretch(beamsize(dlcamp) / 3 / fwhmfac))
imgpr = intensitymap(mpr, grid)
skymeta = (; mimg = imgpr ./ sum(imgpr), nnmodel, ff, st);

using Distributions
# The only other hyperparameters we need to set priors for are the marginal standard deviation
# and the high-frequency tapering parameter. For these we will use an Exponential prior with shape 1.0
# for the standard deviation and a Uniform(0, 1) prior for constant background flux ratio.
prior = (;
    nn = nnprior,
    fb = Uniform(0.0, 1.0),
    σ = Exponential(0.25)
)

# We can then define our sky model.
skym = SkyModel(sky, prior, grid; metadata = skymeta)

# Since we are fitting closures we do not need to include an instrument model, since
# the closure likelihood is approximately independent of gains in the high SNR limit.
using Enzyme
post = VLBIPosterior(skym, dlcamp, dcphase)

# ## Reconstructing the Image

# To optimize our posterior `Comrade` provides the `comrade_opt` function. To use this
# functionality a user first needs to import `Optimization.jl` and the optimizer of choice.
# In this tutorial we will use the AdamW optimizer although others can be used as well.
# We also need to import Enzyme to allow for automatic differentiation.
using Optimization, OptimizationOptimisers
xinit = prior_sample(rng, post) ## draw initial parameters from the prior to reproducible results
imageviz(parent(skymodel(post, xinit)); colorscale = log10, colorrange = (1.0e-8, 1.0e-4), size = (500, 400)) |> DisplayAs.PNG |> DisplayAs.Text
xopt, sol = comrade_opt(post, AdamW(3e-4); maxiters = 5000, initial_params = xinit)

# !!! note
#     We are currently working on making Comrade utilize MLIR and XLA through Reactant. This 
#     will greatly speed up these neural models.

using CairoMakie
using DisplayAs #hide
# We can plot the MAP image using `imageviz`.
imgmap = parent(skymodel(post, xopt))
crange = (1.0e-8, 1.0e-4)
fig = imageviz(imgmap, colorscale = log10, colorrange = crange, size = (650, 500));
DisplayAs.Text(DisplayAs.PNG(fig)) #hide


# To see how well the MAP estimate fits the data we can plot the residuals.
res = Comrade.residuals(post, xopt)
fig = Figure(; size = (800, 300))
plotfields!(fig[1, 1], res[1], :uvdist, :res);
plotfields!(fig[1, 2], res[2], :uvdist, :res);
fig |> DisplayAs.PNG |> DisplayAs.Text


# Overall, the image looks reasonable. Unlike other models in `Comrade`, estimating posterior uncertainty
# with neural fields is highly challenging so we skip that step here.

# We can also compare the Comrade neural field reconstruction to the CLEAN reconstruction of the same data.
cleanf = Base.download("https://www.bu.edu/blazars/VLBA_GLAST/1308/1308+326Q.2021-03-19.IMAP.gz")
# By default this will load the clean components with the beam defined in the FITS header.
mcl = load_clean_components(cleanf)
# We can also choose the load the clean components with a user-defined beam.
mcl_25 = load_clean_components(cleanf, modify(Gaussian(), Stretch(beamsize(dlcamp) / 4 / fwhmfac)))

# Now we can produce the CLEAN images on the same grid as our Comrade reconstruction.
cleanimg = intensitymap(mcl, grid)
cleanimg25 = intensitymap(mcl_25, grid)

fig = Figure(; size = (900, 350));
axs = [Axis(fig[1, j], xreversed = true, aspect = DataAspect()) for j in 1:3]
crange = (1.0e-6, 1.0e-1)
image!(axs[1], imgmap, colormap = :afmhot, colorscale = log10, colorrange = crange); axs[1].title = "Comrade Neural Field"
image!(axs[2], max.(cleanimg, 1.0e-20), colormap = :afmhot, colorscale = log10, colorrange = crange); axs[2].title = "CLEAN (Nominal beam)"
image!(axs[3], max.(cleanimg25, 1.0e-20), colormap = :afmhot, colorscale = log10, colorrange = crange); axs[3].title = "CLEAN (25% beam)"
hidedecorations!.(axs)
fig |> DisplayAs.PNG |> DisplayAs.Text

# From the plot you can see that the neural field reconstruction appears to capture a lot of fine detail.
# Interpreting the reliability of these features is rather challenging without either 1 posterior uncertainty
# estimates or 2 knowledege of what these sources look like. However, this example demonstrates the flexibility of neural fields
# for VLBI imaging. 
#
# A disclaimer is that we have not made any attempt to optimize the neural network architecture or hyperparameters here.
# For a more thorough and insightful analysis a paper by Foschi et al is in preparation.

