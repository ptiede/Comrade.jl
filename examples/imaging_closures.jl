# # Imaging a Black Hole using only Closure Quantities
# The EHT is the highest resolution telescope ever created. Its resolution is equivalent
# to roughly tracking a hockey puck on the moon when viewing it from earth. However,
# the EHT is also a unique interferometer. For one the data is produces is incredible sparse.
# The telescope it only from 8 geographic locations around the planet, each with its unique
# telescope. Additionally, the EHT observes at a much higher frequency than typical interferometers.
# As a result, the EHT data is often poorly calibrated. Meaning there can be large instrumental effects
# often called *gains* that can corrupt our signal. One way to deal with this is to fit quantities
# that are independent of gains, these are often called **closure quantities**. There are two
# type of closure quantities for total intensity observations:
#  1. Closure phases
#  2. Log-closure amplitudes
#  In this tutorial we will explore imaging M87 using just closure quantities. Note that while
# this simplifies the model, we are throwing out additional information we know about the gains
# so often the images are of poorer quality than fitting more directly the visibilities measured
# by the telescope.

using Pkg; Pkg.activate(@__DIR__)

# To get started we will load Comrade
using Comrade


# To download the data visit https://doi.org/10.25739/g85n-f134
# To load the eht-imaging obsdata object we do:
obs = load_ehtim_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))

# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      coherent.
obs = scan_average(obs).add_fractional_noise(0.01)

# For this tutorial we will stick to fitting closure only data, although we can
# get better results by also modeling gains, since closure only modeling is equivalent
# to assuming infinite gain priors.
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs; cut_trivial=true)

# Build the Model. Here we we a struct to hold some caches
# which will speed up imaging
struct Model{C}
    cache::C
    fovx::Float64
    fovy::Float64
    nx::Int
    ny::Int
    function Model(obs::Comrade.EHTObservation, fovx, fovy, nx, ny)
        buffer = IntensityMap(zeros(ny, nx), fovx, fovy, BSplinePulse{3}())
        cache = create_cache(NFFTAlg(obs), buffer)
        return new{typeof(cache)}(cache, fovx, fovy, nx, ny)
    end
end

# For our model we will be using a rasterized image. This can be viewed as something like a
# non-parametric model. As a result of this we will need to use a `modelimage` object to
# store cache information we will need to compute the numerical FT.
function model(θ, metadata)
    (;c) = θ
    (; grid, cache) = metadata
    # Construct the image model
    img = IntensityMap(c, grid)
    return  ContinuousImage(img, cache)
end


npix = 12
fovxy = μas2rad(70.0)
# Now we can feed in the array information to form the cache. We will be using a DFT since
# it is efficient for so few pixels
# We will use a Dirichlet prior to enforce that the flux sums to unity since closures are
# degenerate to total flux.
grid = imagepixels(fovxy, fovxy, npix, npix)
buffer = IntensityMap(zeros(npix,npix), grid)
cache = create_cache(DFTAlg(dlcamp), buffer, BSplinePulse{3}())
metadata = (;grid, cache)

prior = (c = ImageDirichlet(1.0, npix, npix), )

lklhd = RadioLikelihood(model, metadata, dlcamp, dcphase)
post = Posterior(lklhd, prior)

# Transform from simplex space to the unconstrained
tpost = asflat(post)
ℓ = logdensityof(tpost)

# Let's run an optimizer to get a nice starting location
# It turns out that gradients are really helpful here
ndim = dimension(tpost)
using Zygote
# Creates optimization function using Optimization.jl
f = OptimizationFunction(tpost, Optimization.AutoZygote())
# randn(ndim) is a random initialization guess
# nothing just says there are no additional arguments to the optimization function.
prob = OptimizationProblem(f, rand(ndim) .- 0.5, nothing)

ℓ = logdensityof(tpost)
# Find the best fit image! Using LBFGS optimizaer.
sol = solve(prob, LBFGS(); maxiters=2_000, callback=(x,p)->(@info ℓ(x); false), g_tol=1e-1)


# Move from unconstrained space to physical parameter space
xopt = transform(tpost, sol)

# Let's see how the fit looks
img = intensitymap(model(xopt, metadata), fovxy, fovxy, 512, 512)
plot(img)
residual(model(xopt, metadata), dlcamp)
residual(model(xopt, metadata), dcphase)

# Finally we sample with HMC
metric = DiagEuclideanMetric(ndim)
hchain, stats = sample(post, AHMC(;metric, autodiff=AD.ZygoteBackend()), 5000; nadapts=4000, init_params=xopt)


# Plot the mean image and standard deviation image
using StatsBase
samples = mms.(sample(hchain, 500))
imgs = intensitymap.(samples, fovx, fovy, 256, 256)

mimg, simg = mean_and_std(imgs)

 p1 = plot(mimg, title="Mean", clims=(0.0, maximum(mimg)))
p2 = plot(simg,  title="Std. Dev.", clims=(0.0, maximum(mimg)))
p2 = plot(simg./mimg,  title="Fractional Error", xlims=(-25.0,25.0), ylims=(-25.0,25.0))

# Computing information
# ```
# Julia Version 1.8.0
# Commit 742b9abb4d (2022-05-06 12:58 UTC)
# Platform Info:
#   OS: Linux (x86_64-pc-linux-gnu)
#   CPU: 11th Gen Intel(R) Core(TM) i7-1185G7 @ 3.00GHz
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-12.0.1 (ORCJIT, tigerlake)
# ```
