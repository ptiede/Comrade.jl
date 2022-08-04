using Pkg; Pkg.activate(@__DIR__)
#Pkg.add(url="https://github.com/ptiede/RadioImagePriors.jl")
using Comrade
using Distributions
using ComradeOptimization
using ComradeAHMC
using OptimizationBBO
using Plots
using StatsBase
using OptimizationOptimJL
using RadioImagePriors
using DistributionsAD

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
    fovx::Float64
    fovy::Float64
    nx::Int
    ny::Int
    function GModel(obs::Comrade.EHTObservation, fovx, fovy, nx, ny)
        buffer = IntensityMap(zeros(ny, nx), fovx, fovy, BSplinePulse{3}())
        cache = create_cache(DFTAlg(obs), buffer)
        gcache = GainCache(scantable(obs))
        return new{typeof(cache), typeof(gcache)}(cache, gcache, fovx, fovy, nx, ny)
    end
end

function (model::GModel)(θ)
    (;c, f, lgamp) = θ
    # Construct the image model
    img = IntensityMap(f*c, model.fovx, model.fovy, BSplinePulse{3}())
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

fovx = μas2rad(85.0)
fovy = μas2rad(75.0)
nx = 16
ny = floor(Int, fovy/fovx*nx)
prior = (
          c = ImageDirichlet(1.0, ny, nx),
          f = Uniform(0.2, 0.9),
          lgamp = Comrade.GainPrior(distamp, scantable(damp)),
        )


mms = GModel(damp, fovx, fovy, nx, ny)

post = Posterior(lklhd, prior, mms)

tpost = asflat(post)

# We will use HMC to sample the posterior.
# First to get in the right ballpark we will use `BlackBoxOptim.jl`
ndim = dimension(tpost)
using Zygote
f = OptimizationFunction(tpost, Optimization.AutoZygote())
#prob = OptimizationProblem(f, rand(ndim), nothing, lb=fill(-5.0, ndim), ub=fill(5.0, ndim))
#sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1000_000)

# Now lets zoom to the peak using LBFGS
ndim = dimension(tpost)
using Zygote
f = OptimizationFunction(tpost, Optimization.AutoZygote())
prob = OptimizationProblem(f, rand(ndim).-0.5, nothing)
ℓ = logdensityof(tpost)
sol = solve(prob, LBFGS(); maxiters=1000, callback=(x,p)->(@info ℓ(x); false), g_tol=1e-1)

xopt = transform(tpost, sol)

# Let's see how the fit looks
residual(mms(xopt), damp)
residual(mms(xopt), dcphase)
plot(mms(xopt), fovx=fovx, fovy=fovy, title="MAP")

# Let's also plot the gain curves
gt = Comrade.caltable(mms(xopt))
plot(gt, ylims=:none, layout=(3,3), size=(600,500))

using Measurements


# now we sample using hmc
metric = DiagEuclideanMetric(ndim)
hchain, stats = sample(post, AHMC(;metric, autodiff=AD.ZygoteBackend()), 500; nadapts=400, init_params=xopt)

# Now plot the gain table with error bars
gamps = exp.(hcat(hchain.lgamp...))
mga = mean(gamps, dims=2)
sga = std(gamps, dims=2)

using Measurements
gmeas = measurement.(mga, sga)
ctable = caltable(mms.gcache, vec(gmeas))
plot(ctable)

# This takes about 1.75 hours on my laptop. Which isn't bad for a 575 dimensional model!

# Plot the mean image and standard deviation image
using StatsBase
samples = mms.(sample(hchain, 50))
imgs = intensitymap.(samples, fovx*1.1, fovy*1.1, 96, 96)

mimg, simg = mean_and_std(imgs)

p1 = plot(mimg, title="Mean", clims=(0.0, maximum(mimg)))
p2 = plot(simg,  title="Std. Dev.", clims=(0.0, maximum(mimg)))
p2 = plot(simg./mimg,  title="Fractional Error")

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
