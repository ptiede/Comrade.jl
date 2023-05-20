using Pkg; Pkg.activate(@__DIR__)
Pkg.add(url="https://github.com/ptiede/RadioImagePriors.jl")
using Comrade
using Distributions
using Plots
using ComradeOptimization
using OptimizationOptimJL
using OptimizationBBO
using ComradeAHMC
using RadioImagePriors

# load eht-imaging we use this to load eht data
# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obs.add_scans()
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obs = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true).add_fractional_noise(0.01)
# extract log closure amplitudes and closure phases
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)
lklhd = RadioLikelihood(dlcamp, dcphase)

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


npix = 24
fovxy = μas2rad(62.5)
# Now we can feed in the array information to form the cache. We will be using a DFT since
# it is efficient for so few pixels
mim = ImModel(dlcamp, fovxy, npix)
# We will use a Dirichlet prior to enforce that the flux sums to unity since closures are
# degenerate to total flux.
im_pr = (c = ImageDirichlet(0.5, npix, npix),)

function mring(θ)
    (;ma, mp, r, τ, σ, ξ) = θ
    s,c = sincos(mp)
    α = ma*c
    β = ma*s
    return rotated(smoothed(stretched(MRing(α, β), r, r*τ),σ), ξ)
end

mring_pr = ( r = Uniform(μas2rad(10.0), μas2rad(30.0)),
             σ = Uniform(μas2rad(0.5), μas2rad(20.0)),
             τ = Uniform(0.5, 1.0),
             ξ = Uniform(-π/2, π/2),
             ma = Uniform(0.0, 0.5),
             mp = Uniform(0.0, 2π),
            )


struct HModel{C,F}
    cache::C
    fov::F
    npix::Int
    function HModel(obs::Comrade.EHTObservation, fov::Real, npix::Int)
        buffer = IntensityMap(zeros(npix, npix), fov, fov, BSplinePulse{3}())
        cache = create_cache(DFTAlg(obs), buffer)
        return new{typeof(cache), typeof(fov)}(cache, fov, npix)
    end
end

# Now define the actual model creation using a Julia "functor".
function (model::HModel)(θ)
    (;c, f, r, σ, τ, ξ, ma, mp) = θ
    img = IntensityMap(f*c, model.fov, model.fov, BSplinePulse{3}())
    mimg = modelimage(img, model.cache)
    s,c = sincos(mp)
    α = ma*c
    β = ma*s
    ring = stretched(MRing(α, β), r, r*τ)
    return ring + mimg
end

mhb = HModel(dlcamp, fovxy, npix)


# Now form model we are going to fit

hb_pr = (
          c = ImageDirichlet(0.5, npix, npix),
          f = Uniform(0.0, 1.0),
          r = Uniform(μas2rad(10.0), μas2rad(30.0)),
          σ = Uniform(μas2rad(0.5), μas2rad(20.0)),
          τ = Uniform(0.5, 1.0),
          ξ = Uniform(-π/2, π/2),
          ma = Uniform(0.0, 0.5),
          mp = Uniform(0.0, 2π),
        )

# Now we can build the posterior
post_mr = Posterior(lklhd, mring_pr, mring)
post_im = Posterior(lklhd, im_pr, mim)
post_hb = Posterior(lklhd, hb_pr, mhb)
# We will use HMC to explore the posterior. However, we can really help sampling if we first
# get a good starting location. As such we will first optimize. First lets transform the
# posterior parameters to the unconstrained space.
tpost_mr = asflat(post_mr)
tpost_im = asflat(post_im)
tpost_hb = asflat(post_hb)

x0_mr = prior_sample(tpost_mr)
x0_im = prior_sample(tpost_im)
x0_hb = prior_sample(tpost_hb)

ℓmr = logdensityof(tpost_mr)
ℓim = logdensityof(tpost_im)
ℓhb = logdensityof(tpost_hb)

using BenchmarkTools
using Zygote
@benchmark ℓmr(x0_mr)
@benchmark ℓim(x0_im)
@benchmark ℓhb(x0_hb)


@benchmark ℓmr'($(x0_mr))
@benchmark ℓim'($(x0_im))
@benchmark ℓhb'($(x0_hb))


ndim = dimension(tpost)

# Now we optimize. First we will use BlackBoxOptim which is a genetic algorithm to get us
# in the region of the best fit model.
f = OptimizationFunction(tpost, Optimization.AutoForwardDiff())
prob = OptimizationProblem(f, randn(ndim), nothing, lb=fill(-5.0, ndim), ub=fill(5.0, ndim))
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=100_000)

# Alright now we can zoom to the peak! But we can even do better, we will also compute
# the laplace approximation to the posterior as a byproduct.
prob = OptimizationProblem(f, sol.u, nothing)
ℓ = logdensityof(tpost)
sol = solve(prob, LBFGS(), maxiters=1_000, callback=((x,p)->(@info ℓ(x);false)), g_tol=1e-1)

# tranform back to parameter space
xopt = transform(tpost, sol)

# plot the residuals and the intensity map
residual(mms(xopt), dlcamp)
residual(mms(xopt), dcphase)
img = intensitymap(mms(xopt), μas2rad(160.0), μas2rad(160.0), 512, 512)
plot(abs.(img) , xlims=(-60.0,60.0), ylims=(-60.0,60.0))

# now we sample using hmc
metric = DenseEuclideanMetric(ndim)
chain, stats = sample(post, AHMC(;metric, autodiff=Val(:ForwardDiff)), 4000; nadapts=3000, init_params=xopt)

# this took 25 minutes on my laptop which has a 11 gen core i7

# Plot the mean image and standard deviation
# Adaptation ruins detailed balance/reversibility in the chain!
# We will also split the uncertainty and mean maps into ring and raster
using StatsBase
msamples = Comrade.components.(mms.(sample(chain, 500)))
ring_samples = first.(msamples)
rast_samples = last.(msamples)
ring_imgs = intensitymap.(ring_samples, fovxy, fovxy, 256, 256)
rast_imgs = intensitymap.(rast_samples, fovxy, fovxy, 256, 256)

ring_mean, ring_std = mean_and_std(ring_imgs)
rast_mean, rast_std = mean_and_std(rast_imgs)
both_mean, both_std = mean_and_std(ring_imgs .+ rast_imgs)

p1 = plot(ring_mean, title="Ring Mean", clims=(0.0, maximum(ring_mean)))
p2 = plot(ring_std,  title="Ring Std. Dev.", clims=(0.0, maximum(ring_mean)))
p3 = plot(rast_mean, title="Raster Mean", clims=(0.0, maximum(rast_mean)))
p4 = plot(rast_std,  title="Raster Std. Dev.", clims=(0.0, maximum(rast_mean)))
p5 = plot(both_mean, title="Both Mean", clims=(0.0, maximum(both_mean)))
p6 = plot(both_std,  title="Both Std. Dev.", clims=(0.0, maximum(both_mean)))

plot(p1,p2,p3,p4,p5,p6, layout=(3,2), size=(950,1000))

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
