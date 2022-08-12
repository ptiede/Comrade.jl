using Pkg; Pkg.activate(@__DIR__)
#Pkg.add(url="https://github.com/ptiede/RadioImagePriors.jl")
using Comrade
using Distributions
using ComradeOptimization
using ComradeAHMC
using OptimizationBBO
using OptimizationOptimJL

using Plots
using StatsBase
using RadioImagePriors

# load eht-imaging we use this to load eht data
load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obs.add_scans()
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obs = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true).add_fractional_noise(0.02)
# extract log closure amplitudes and closure phases
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)

# Build the Model. Here we we a struct to hold some caches
# which will speed up imaging
struct ImModel{C}
    cache::C
    fov::Float64
    npix::Int
    function ImModel(obs::Comrade.EHTObservation, fov::Real, npix::Int)
        buffer = IntensityMap(zeros(npix, npix), fov, fov, BSplinePulse{3}())
        cache = create_cache(DFTAlg(obs), buffer)
        return new{typeof(cache)}(cache, fov, npix)
    end

end

# For our model we will be using a rasterized image. This can be viewed as something like a
# non-parametric model. As a result of this we will need to use a `modelimage` object to
# store cache information we will need to compute the numerical FT.
function (model::ImModel)(θ)
    (;c) = θ
    #Construct the image model
    img = IntensityMap(c, model.fov, model.fov, BSplinePulse{3}())
    #Create the modelimage object that will use a cache to compute the DFT
    m = modelimage(img, model.cache)
end


npix = 16
fovxy = μas2rad(70.0)
# Now we can feed in the array information to form the cache. We will be using a DFT since
# it is efficient for so few pixels
mms = ImModel(dlcamp, fovxy, npix)
# We will use a Dirichlet prior to enforce that the flux sums to unity since closures are
# degenerate to total flux.
img = IntensityMap(zeros(npix,npix), fovxy, fovxy)
xitr, yitr = Comrade.imagepixels(img)
prior = (c = CenteredImage(xitr, yitr, μas2rad(1.0),
            ImageDirichlet(1.0, npix, npix)
            )
            ,)

lklhd = RadioLikelihood(mms, dlcamp, dcphase)
post = Posterior(lklhd, prior)

# Transform from simplex space to the unconstrained
tpost = asflat(post)

# Let's run an optimizer to get a nice starting location
# It turns out that gradients are really helpful here
ndim = dimension(tpost)
using Zygote
# Creates optimization function using Optimization.jl
f = OptimizationFunction(tpost, Optimization.AutoZygote())
# randn(ndim) is a random initialization guess
# nothing just says there are no additional arguments to the optimization function.
prob = OptimizationProblem(f, -rand(ndim), nothing)

ℓ = logdensityof(tpost)
# Find the best fit image! Using LBFGS optimizaer.
sol = solve(prob, LBFGS(); maxiters=2_000, callback=(x,p)->(@info ℓ(x); false), g_tol=1e-1)

# Move from unconstrained space to physical parameter space
xopt = transform(tpost, sol)

# Let's see how the fit looks
residual(mms(xopt), dlcamp)
residual(mms(xopt), dcphase)
plot(mms(xopt), fovx=fovxy, fovy=fovxy, title="MAP")


# now we sample using hmc
metric = DiagEuclideanMetric(ndim)
hchain, stats = sample(post, AHMC(;metric, autodiff=AD.ZygoteBackend()), 4000; nadapts=3000, init_params=xopt)

# This takes about 1.75 hours on my laptop. Which isn't bad for a 575 dimensional model!

# Plot the mean image and standard deviation image
 using StatsBase
samples = mms.(sample(hchain, 500))
imgs = intensitymap.(samples, fovxy, fovxy, 96, 96)

mimg, simg = mean_and_std(imgs)

p1 = plot(mimg, title="Mean", clims=(0.0, maximum(mimg)))
p2 = plot(simg,  title="Std. Dev.", clims=(0.0, maximum(mimg)))
p2 = plot(simg./mimg,  title="Fractional Error", xlims=(-25.0,25.0), ylims=(-25.0,25.0))

# Computing information
# ```
# Julia Version 1.7.3
# Commit 742b9abb4d (2022-05-06 12:58 UTC)
# Platform Info:
#   OS: Linux (x86_64-pc-linux-gnu)
#   CPU: 11th Gen Intel(R) Core(TM) i7-1185G7 @ 3.00GHz
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-12.0.1 (ORCJIT, tigerlake)
# ```
