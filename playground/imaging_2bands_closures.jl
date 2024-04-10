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
# To download the data visit https://doi.org/10.25739/g85n-f134
obslo = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
obshi = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obslo.add_scans()
obshi.add_scans()
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obslo = obslo.avg_coherent(0.0, scan_avg=true).add_fractional_noise(0.01)
obshi = obshi.avg_coherent(0.0, scan_avg=true).add_fractional_noise(0.01)
# extract log closure amplitudes and closure phases

# Build the Model. Here we we a struct to hold some caches
# This will be useful to hold precomputed caches

struct GModel1{C}
    cache::C
    fovx::Float64
    fovy::Float64
    nx::Int
    ny::Int
    function GModel1(obs::Comrade.EHTObservation, fovx, fovy, nx, ny)
        buffer = IntensityMap(zeros(ny, nx), fovx, fovy, (0.0,0.0), BSplinePulse{3}())
        cache = create_cache(DFTAlg(obs), buffer)
        return new{typeof(cache)}(cache, fovx, fovy, nx, ny)
    end
end

function (model::GModel1)(θ)
    (;c, fg) = θ
    # Construct the image model
    img = IntensityMap((1-fg)*c, model.fovx, model.fovy, (0.0, 0.0), BSplinePulse{3}())
    m = modelimage(img, model.cache)
    gaussian = fg*stretched(Gaussian(), μas2rad(1000.0), μas2rad(1000.0))
    # Now corrupt the model with Gains
    return m+gaussian
end

struct GModel2{C}
    cache::C
    fovx::Float64
    fovy::Float64
    nx::Int
    ny::Int
    function GModel2(obs::Comrade.EHTObservation, fovx, fovy, nx, ny)
        buffer = IntensityMap(zeros(ny, nx), fovx, fovy, (0.0,0.0), BSplinePulse{3}())
        cache = create_cache(DFTAlg(obs), buffer)
        return new{typeof(cache)}(cache, fovx, fovy, nx, ny)
    end
end

function (model::GModel2)(θ)
    (;c, fg) = θ
    # Construct the image model
    img = IntensityMap((1-fg)*c, model.fovx, model.fovy, (0.0, 0.0), BSplinePulse{3}())
    m = modelimage(img, model.cache)
    gaussian = fg*stretched(Gaussian(), μas2rad(1000.0), μas2rad(1000.0))
    # Now corrupt the model with Gains
    return m+gaussian
end



# First we define the sites gain priors
fovx = μas2rad(75.0)
fovy = μas2rad(75.0)
nx = 12
ny = floor(Int, fovy/fovx*nx)
prior = (
          c = ImageDirichlet(1.0, ny, nx),
          fg = Uniform(0.2, 1.0),
        )


mms1 = GModel1(damplo, fovx, fovy, nx, ny)
mms2 = GModel2(damphi, fovx, fovy, nx, ny)
lklhd1 = RadioLikelihood(mms1, damplo, dcphaselo)
lklhd2 = RadioLikelihood(mms2, damphi, dcphasehi)
lklhd = MultiRadioLikelihood(lklhd1, lklhd2)

post = Posterior(lklhd, prior)

tpost = asflat(post)

# We will use HMC to sample the posterior.

ndim = dimension(tpost)
using Zygote
f = OptimizationFunction(tpost, Optimization.AutoZygote())
prob = OptimizationProblem(f, rand(ndim) .- 0.5, nothing)
ℓ = logdensityof(tpost)
sol = solve(prob, LBFGS(); maxiters=2000, callback=(x,p)->(@info ℓ(x); false), g_tol=1e-1)

xopt = transform(tpost, sol)

# Let's see how the fit looks
residual(mms1(xopt), damplo)
residual(mms1(xopt), dcphaselo)
residual(mms2(xopt), damphi)
residual(mms2(xopt), dcphasehi)
plot(mms2(xopt), fovx=fovx, fovy=fovy, title="MAP")


# now we sample using hmc
metric = DiagEuclideanMetric(ndim)
hchain, stats = sample(post, AHMC(;metric, autodiff=Val(:Zygote)), 2000; n_adapts=1000, initial_params=xopt)

# Now plot the gain table with error bars
gamps1 = exp.(hcat(hchain.lgamp1...))
mga1 = mean(gamps1, dims=2)
sga1 = std(gamps1, dims=2)

gamps2 = exp.(hcat(hchain.lgamp2...))
mga2 = mean(gamps2, dims=2)
sga2 = std(gamps2, dims=2)


using Measurements
gmeas1 = measurement.(mga1, sga1)
gmeas2 = measurement.(mga2, sga2)
ctable1 = caltable(mms1.gcache, vec(gmeas1))
ctable2 = caltable(mms2.gcache, vec(gmeas2))
plot(ctable1, layout=(3,3), size=(600,500))
plot!(ctable2, layout=(3,3), size=(600,500))

# This takes about 1.75 hours on my laptop. Which isn't bad for a 575 dimensional model!

# Plot the mean image and standard deviation image
using StatsBase
samples = mms1.(sample(hchain, 50))
imgs = intensitymap.(samples, fovx, fovy, 96, 96)

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
