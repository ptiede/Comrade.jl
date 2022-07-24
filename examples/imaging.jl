using Pkg; Pkg.activate(@__DIR__)
Pkg.add(url="https://github.com/ptiede/RadioImagePriors.jl")
using Comrade
using Distributions
using ComradeOptimization
using ComradeAHMC
using OptimizationBBO
using Plots
using StatsBase
using OptimizationOptimJL
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
damp = extract_amp(obs)
dcphase = extract_cphase(obs)
lklhd = RadioLikelihood(damp, dcphase)

# Build the Model. Here we we a struct to hold some caches
# This will be useful to hold precomputed caches

struct GModel{C,G}
    cache::C
    gcache::G
    fov::Float64
    npix::Int
    function GModel(obs::Comrade.EHTObservation, fov, npix)
        buffer = IntensityMap(zeros(npix, npix), fov, fov, BSplinePulse{3}())
        cache = create_cache(DFTAlg(obs), buffer)
        gcache = GainCache(scantable(obs))
        return new{typeof(cache), typeof(gcache)}(cache, gcache, fov, npix)
    end
end

function (model::GModel)(θ)
    (;c, f, lgamp) = θ
    # Construct the image model
    img = IntensityMap(f*c, model.fov, model.fov, BSplinePulse{3}())
    m = modelimage(img, model.cache)
    # Now corrupt the model with Gains
    g = exp.(lgamp)
    Comrade.GainModel(model.gcache, g, m)
end

# First we define the station gain priors
distamp = (AA = Normal(0.0, 0.1),
           AP = Normal(0.0, 0.1),
           LM = Normal(0.0, 0.9),
           AZ = Normal(0.0, 0.1),
           JC = Normal(0.0, 0.1),
           PV = Normal(0.0, 0.1),
           SM = Normal(0.0, 0.1)
           )

npix = 16
fovxy = μas2rad(70.0)
prior = (
          c = ImageDirichlet(0.5, npix, npix),
          f = Uniform(0.2, 0.9),
          lgamp = Comrade.GainPrior(distamp, scantable(damp)),
        )


mms = GModel(damp, fovxy, npix)

post = Posterior(lklhd, prior, mms)

tpost = asflat(post)

# We will use HMC to sample the posterior.
# First to get in the right ballpark we will use `BlackBoxOptim.jl`
ndim = dimension(tpost)
f = OptimizationFunction(tpost, Optimization.AutoZygote())
prob = OptimizationProblem(f, rand(ndim), nothing, lb=fill(-5.0, ndim), ub=fill(5.0, ndim))
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1000_000)

# Now lets zoom to the peak using LBFGS
ndim = dimension(tpost)
using Zygote
f = OptimizationFunction(tpost, Optimization.AutoZygote())
prob = OptimizationProblem(f, randn(ndim), nothing)
ℓ = logdensityof(tpost)
sol = solve(prob, LBFGS(); maxiters=2_000, callback=(x,p)->(@info ℓ(x); false), g_tol=1e-1)

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


ℓ = logdensityof(tpost)
prob = OptimizationProblem(f, randn(ndim), nothing)
sol = solve(prob, LBFGS(); maxiters=2_000, callback=(x,p)->(@info ℓ(x); false), g_tol=1e-1)

xopt = transform(tpost, sol)

residual(mms(xopt), damp)
residual(mms(xopt), dcphase)

plot(mms(xopt), fovx=fovxy, fovy=fovxy)

metric = DenseEuclideanMetric(ndim)
chain, stats = sample(post, AHMC(;metric, autodiff=AD.ZygoteBackend()), 4000; nadapts=3000, init_params=xopt)
