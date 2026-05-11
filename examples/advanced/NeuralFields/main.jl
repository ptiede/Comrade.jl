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
using VLBIFiles
using Lux

# We will accelerate the optimization by tracing the forward+gradient through MLIR/XLA
# with [Reactant.jl](https://github.com/EnzymeAD/Reactant.jl). Reactant precompiles a single
# loss-and-gradient step into an XLA program, which is then re-used for every iteration of
# our optimizer.
using Reactant
using VLBISkyModels: ReactantNUFFTAlg
using Comrade: ReactantEx

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(123)


# ## Load the Data
# For this tutorial we will image publicly available VLBA data of the AGN
# 1308+326 observed on 2021/03/19 at 43 GHz as part of the Boston University blazar monitoring program.
file = Base.download("https://www.bu.edu/blazars/VLBA_GLAST/1308/1308+326Q.2021-03-19.UVP.gz")
uvd = VLBIFiles.load(
    VLBIFiles.UVData,
    file
)


# For this tutorial we will only use closure quantities to reconstruct the image however, polarized
# or complex visibilities can also be used with instrumental models following the other tutorials.
dvis = extract_table(
    uvd,
    Visibilities(; time_average = VLBI.GapBasedScans()),
)
add_fractional_noise!(dvis, 0.005)



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
    freqs = ntuple(n -> 2^(n - 1), m)
    feats = zeros(Float32, 4 * m, nx, ny)
    x = range(-1.0f0 * π, 1.0f0 * π, length = nx + 1)[begin:(end - 1)] #clip the endpoint
    y = range(-1.0f0 * π, 1.0f0 * π, length = ny + 1)[begin:(end - 1)] #clip the endpoint
    for i in 1:nx, j in 1:ny, (ik, k) in enumerate(1:4:(4 * m))
        s, c = sincos(freqs[ik] * x[j])
        feats[k, j, i] = s
        feats[k + 1, j, i] = c
        s, c = sincos(freqs[ik] * y[i])
        feats[k + 2, j, i] = s
        feats[k + 3, j, i] = c
    end
    return reshape(feats, 4 * m, nx * ny)
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
    Dense(size(ff, 1) => 64, Lux.tanh_fast),
    Dense(64 => 64, Lux.tanh_fast),
    Dense(64 => 64, Lux.tanh_fast),
    Dense(64 => 1, use_bias = false),
)


# We can now define our sky model function. This function takes in the neural network parameters `nn`
# and hyperparameters `fb, σ` and returns a `ContinuousImage` representing the sky brightness distribution.
# This should look similar to essentially every other imaging model in the Comrade tutorials.
function sky(θ, metadata)
    (; nn, σ, fb, ftot) = θ
    (; nnmodel, mimg, ff, st) = metadata
    fbn = fb / length(mimg)
    bimg = baseimage(mimg)
    mb = bimg .* (1 - fb) .+ fbn

    ## Compute the neural field
    δ = reshape(first(Lux.apply(nnmodel, ff, nn, st)), size(mimg))
    δ .*= σ 
    mδ = maximum(δ)
    rast = @. exp(δ - mδ) * mb
    rast .*= ftot/sum(rast)
    return ContinuousImage(rast, axisdims(mimg), DeltaPulse{Float64}())
end


# For the prior on the network parameters we use the so-called *NNGP parameterization* of
# [Neal (1996)](@cite) and [Lee et al. (2018)](@cite). For a `Dense` layer with weight
# matrix `W ∈ ℝ^{out × in}` and bias `b ∈ ℝ^{out}` we use
#
#   W_{ij} ∼ 𝒩(0, σ_w² / fan_in),   b_i ∼ 𝒩(0, σ_b²)
#
# i.e. the weight variance is scaled by the inverse fan-in. With this scaling the pre-activation
# variance is independent of the layer width, and as the hidden widths are taken to infinity the
# network output converges to a well-defined Gaussian process (the NNGP). With unit-variance
# `StdNormal` priors on every parameter the pre-activation variance grows linearly with width and
# no proper GP limit exists, so this scaling is the principled choice for a Bayesian neural field.
# We use `Lux.setup` to get the parameter layout of the network and then use `Comrade.rmap` to
# walk the parameter tree, dispatching on parameter rank: rank-2 arrays are weight matrices (with
# fan-in equal to `size(W, 2)`) and rank-1 arrays are biases.
using Random
ps, st = Lux.setup(Random.default_rng(), nnmodel) |> reactant_device()
function nngp_prior(x::AbstractArray)
    if ndims(x) == 2
        fan_in = size(x, 2)
        return VLBIGaussian(zero(eltype(x)), one(eltype(x)), size(x))
    else
        return VLBIImagePriors.StdNormal(size(x))
    end
end
nnprior = Comrade.rmap(nngp_prior, ps)

# Similar to the RF models we will use a so-called `mean image`. Note that in this case is isn't
# actually clean what the mean image is, but this will help with centroid drift so we include it.
# We will use a symmetric Gaussian with a FWHM equal to the approximate
# beamsize of the array. This models the fact that we expect the AGN core to be compact.
fwhmfac = 2 * sqrt(2 * log(2))
mpr = modify(TBlob(4.0), Stretch(beamsize(dvis) / 3 / fwhmfac))
imgpr = intensitymap(mpr, grid)
skymeta = (; mimg = imgpr ./ sum(imgpr), nnmodel, ff, st);

using Distributions
# The only other hyperparameters we need to set priors for are the marginal standard deviation
# and the high-frequency tapering parameter. For these we will use an Exponential prior with shape 1.0
# for the standard deviation and a Uniform(0, 1) prior for constant background flux ratio.
prior = (;
    nn = nnprior,
    fb = VLBIUniform(0.0, 1.0),
    σ = VLBIExponential(0.5),
    ftot = VLBIUniform(0.1, 10.0)
)

# We can then define our sky model. We pass `algorithm = ReactantAlg()` so the non-uniform
# FFT used to evaluate visibilities runs through Reactant's XLA-traced NUFT plan rather
# than the default CPU `NFFTAlg()`.
skym = SkyModel(sky, prior, grid; algorithm = VLBISkyModels.ReactantNUFFTAlg(), metadata = skymeta)


# Now we will fit gains since our data product are complex visibilities
g(x) = exp(complex(x.lg, x.gp))
G = SingleStokesGain(g)
intpr = (
    lg=ArrayPrior(
        IIDSitePrior(IntegSeg(), VLBIGaussian(0.0, 0.5));
    ),
    gp=ArrayPrior(
        IIDSitePrior(IntegSeg(), DiagonalVonMises(0.0, inv(0.5^2)));
        refant=SEFDReference(0.0),
        phase=true,
    ),
)
intmodel = InstrumentModel(G, intpr)



# Since we are fitting closures we do not need to include an instrument model, since
# the closure likelihood is approximately independent of gains in the high SNR limit.
using Enzyme
post_cpu = VLBIPosterior(skym, intmodel, dvis)

# Move every array Reactant cares about (sky-model grid + metadata, data tables, etc.)
# onto the Reactant device. `prepare_device` walks the posterior tree and replaces the
# concrete arrays with `Reactant.RArray`s, leaving priors and the (empty) instrument
# model untouched.
post = Comrade.prepare_device(post_cpu, ComradeBase.ReactantEx())

# ## Reconstructing the Image

# Because the optimizer needs to traverse the entire parameter space (network weights +
# `σ`, `fb`) we work on the unconstrained-flat representation produced by `asflat`. This
# wraps the posterior in a `TransformedVLBIPosterior` whose log-density takes a flat
# `Vector{Float64}` and includes the log-Jacobian of the bijection.
tpost = asflat(post)

# Draw initial parameters from the prior, push them through the inverse transform to get a
# flat vector, and visualize the corresponding initial image. `skymodel` traces through
# Reactant because `post` lives on the device, so we `@jit` the call and then materialize
# the underlying raster back onto the host before handing it to Makie.
using CairoMakie, DisplayAs
xinit_nt = prior_sample(rng, post) ## NamedTuple in the original parameter space
init_img_r = @jit skymodel(post, Reactant.to_rarray(xinit_nt))
init_img = Comrade.Adapt.adapt(Array, init_img_r.img)
imageviz(init_img, size = (500, 400), colorscale=log10, colorrange=(1e-6, 1e-4)) |> DisplayAs.PNG |> DisplayAs.Text

# Move the parameter vector onto the Reactant device.
xinit_r = Reactant.to_rarray(Comrade.inverse(tpost, xinit_nt))

# We use [Optimisers.jl](https://github.com/FluxML/Optimisers.jl) directly (rather than the
# `comrade_opt` wrapper) because Reactant compiles the whole step — including the optimizer
# update — into a single XLA program. AdamW with a learning rate of `3.0e-4` matches the
# Enzyme/CPU version of this tutorial.
using Optimisers
opt = Optimisers.AdamW(3.0e-3)
opt_state = @jit Optimisers.setup(opt, xinit_r)

# A single training step: compute the negative log-density and its Enzyme reverse-mode
# gradient, then apply one Optimisers.jl update. Reactant traces this whole function once
# and reuses the compiled XLA program for every iteration.
function step(tpost, x, opt_state)
    dx = Enzyme.make_zero(x)
    _, primal = Enzyme.autodiff(
        Enzyme.WithPrimal(Enzyme.set_strong_zero(Enzyme.Reverse)),
        Comrade.logdensityof, Active, Const(tpost), Duplicated(x, dx)
    )
    ## Optimisers expects gradients of the loss we are *minimising*; logdensityof returns
    ## the log-posterior so we negate.
    grads = -1 .* dx
    new_state, new_x = Optimisers.update(opt_state, x, grads)
    return new_state, new_x, -primal
end

# Compile the step. The first call is the (slow) MLIR/XLA compile; every subsequent call is
# a single dispatched program.
step_jit = @compile sync = true step(tpost, xinit_r, opt_state)

# Run the optimizer. We log the loss every 250 steps so the docs build still shows
# convergence without flooding stdout.
xcur = xinit_r
state_cur = opt_state
maxiters = 10_000
for i in 1:maxiters
    state_cur, xcur, loss = step_jit(tpost, xcur, state_cur)
    if i == 1 || i % 250 == 0 || i == maxiters
        @info "iter $i  -logdensity = $(Reactant.@allowscalar Float64(loss))"
    end
end

# Pull the optimized parameters back to the host as a NamedTuple in the original space.
xopt = @jit Comrade.transform(tpost, xcur);

using CairoMakie
using DisplayAs #hide
# We can plot the MAP image using `imageviz`.
imgmap = parent(@jit skymodel(post, xopt));
crange = (1.0e-7, 1.0e-3)
fig = imageviz(Comrade.Adapt.adapt(Array, imgmap), colorscale = log10, colorrange = crange, size = (650, 500));
DisplayAs.Text(DisplayAs.PNG(fig)) #hide


# To see how well the MAP estimate fits the data we can plot the residuals.
xopt_cpu = Comrade.transform(asflat(post_cpu), Array(xcur))
res = @jit(Comrade.residuals(post, xopt))[1]
res_cpu = Comrade.Adapt.adapt(Array, res)
fig = Figure(; size = (800, 300))
plotfields!(fig[1, 1], res_cpu, uvdist, :res);
fig |> DisplayAs.PNG |> DisplayAs.Text


# Overall, the image looks reasonable. Unlike other models in `Comrade`, estimating posterior uncertainty
# with neural fields is highly challenging so we skip that step here.

# We can also compare the Comrade neural field reconstruction to the CLEAN reconstruction of the same data.
cleanf = Base.download("https://www.bu.edu/blazars/VLBA_GLAST/1308/1308+326Q.2021-03-19.IMAP.gz")
# By default this will load the clean components with the beam defined in the FITS header.
mcl = load_clean_components(cleanf)
# We can also choose the load the clean components with a user-defined beam.
mcl_25 = load_clean_components(cleanf, modify(Gaussian(), Stretch(beamsize(dvis) / 4 / fwhmfac)))

# Now we can produce the CLEAN images on the same grid as our Comrade reconstruction.
cleanimg = intensitymap(mcl, grid)
cleanimg25 = intensitymap(mcl_25, grid)

fig = Figure(; size = (900, 350));
axs = [Axis(fig[1, j], xreversed = true, aspect = DataAspect()) for j in 1:3]
crange = (1.0e-4, 1.0e-1)
image!(axs[1], Comrade.Adapt.adapt(Array, imgmap), colormap = :afmhot, colorscale = log10, colorrange = crange); axs[1].title = "Comrade Neural Field"
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










