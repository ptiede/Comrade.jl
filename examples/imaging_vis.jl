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
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "0316+413.2013.08.26.uvfits"))
obs.add_scans()
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obsavg = obs.avg_coherent(60.0).add_fractional_noise(0.05)
# extract log closure amplitudes and closure phases
dvis = extract_vis(obsavg)

# Build the Model. Here we we a struct to hold some caches
# This will be useful to hold precomputed caches

struct GModel{C,G}
    cache::C
    gcache::G
    fovx::Float64
    fovy::Float64
    x0::Float64
    y0::Float64
    nx::Int
    ny::Int
    function GModel(obs::Comrade.EHTObservation, fovx, fovy, x0, y0, nx, ny)
        buffer = IntensityMap(zeros(ny, nx), (fovx, fovy), (x0, y0), BSplinePulse{3}())
        cache = create_cache(NFFTAlg(obs), buffer)
        gcache = GainCache(scantable(obs))
        return new{typeof(cache), typeof(gcache)}(cache, gcache, fovx, fovy, x0, y0, nx, ny)
    end
end

function (model::GModel)(θ)
    (;c, f, lgamp, gphase) = θ
    # Construct the image model
    img = IntensityMap(f*c, (model.fovx, model.fovy), (model.x0, model.y0), BSplinePulse{3}())
    m = modelimage(img, model.cache)
    # Now corrupt the model with Gains
    # Now corrupt the model with Gains
    g = exp.(lgamp)
    gp = cis.(gphase)
    return Comrade.GainModel(model.gcache, g.*gp, m)
    #return m
end

# First we define the station gain priors
sites = stations(dvis)
distamp = ( BR = Normal(0.0, 0.2),
            FD = Normal(0.0, 0.2),
            HN = Normal(0.0, 0.2),
            KP = Normal(0.0, 0.2),
            LA = Normal(-0.05, 0.2),
            MK = Normal(0.3, 0.6),
            NL = Normal(0.0, 0.2),
            OV = Normal(0.0, 0.2),
            PT = Normal(-0.1, 0.6),
            SC = Normal(0.3, 0.5))
distphase = NamedTuple{Tuple(sites[2:end])}(Tuple(Normal(0.0, 0.5) for _ in sites[2:end]))
distphase2 = merge(NamedTuple{(first(sites),)}((Normal(0.0, 0.001),)), distphase)


fovx = μas2rad(2800.0)
fovy = μas2rad(3750.0)
x0 = μas2rad(000.0)
y0 = μas2rad(1500.0)
nx = 76
ny = floor(Int, fovy/fovx*nx)+1
prior = (
          c = ImageDirichlet(1.0, ny, nx),
          f = Uniform(10.0, 20.0),
          lgamp = GainPrior(distamp, scantable(dvis)),
          gphase = GainPrior(distphase2, scantable(dvis))
        )


mms = GModel(dvis, fovx, fovy, x0, y0, nx, ny)
lklhd = RadioLikelihood(mms, dvis)

post = Posterior(lklhd, prior)

tpost = asflat(post)

# We will use HMC to sample the posterior.

ndim = dimension(tpost)
using Zygote
f = OptimizationFunction(tpost, Optimization.AutoZygote())
prob = OptimizationProblem(f, randn(ndim), nothing)
ℓ = logdensityof(tpost)
sol = solve(prob, LBFGS(); maxiters=1500, callback=(x,p)->(@info ℓ(x); false), g_tol=1e-1)

xopt = transform(tpost, sol)

# Let's see how the fit looks
plot(mms(xopt), fovx=fovx, fovy=fovy, title="MAP", phasecenter=(x0, y0), colorbar_scale=:log10, clims=(2e-5,1e-2), size=(600,600))
residual(mms(xopt), dvis)
#residual(mms(xopt), dcphase)

# Let's also plot the gain curves
gt = Comrade.caltable(mms.gcache, xopt.gphase)
plot(gt, layout=(4,3), size=(600,500))

gt = Comrade.caltable(mms.gcache, exp.(xopt.lgamp))
plot(gt, layout=(4,3), size=(600,500))


using Measurements


# now we sample using hmc
metric = DiagEuclideanMetric(ndim)
hchain, stats = sample(post, AHMC(;metric, autodiff=AD.ZygoteBackend()), 5000; nadapts=4000, init_params=xopt)

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
imgs = intensitymap.(samples, fovx, fovy, 128,  128; phasecenter=(x0, y0))

mimg, simg = mean_and_std(imgs)

p1 = plot(mimg, title="Mean", colorbar_scale=:log10, clims=(5e-4, 1e-1))
p2 = plot(simg,  title="Std. Dev.", colorbar_scale=:log10, clims=(5e-4, 1e-1))
p3 = plot(mimg./simg,  title="SNR", clims=(2.0, 17.0))
p4 = plot(simg./mimg,  title="Fractional Error")

plot(p1,p2,p3,p4, layout=(2,2), size=(800,800))
savefig("3c84_complex_vis_2min_avg.png")
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
