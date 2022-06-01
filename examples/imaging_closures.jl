using Pkg; Pkg.activate(@__DIR__)
using Comrade
using Distributions
using ComradeGalactic
using ComradeAHMC
using GalacticBBO
using Plots
using StatsBase
using GalacticOptimJL
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
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)
lklhd = RadioLikelihood(dlcamp, dcphase)

# Build the Model. Here we we a struct to hold some caches
# which will speed up imaging
struct ImModel{C}
    cache::C
    fov::Float64
    npix::Int
end

# For our model we will be using a rasterized image. This can be viewed as something like a
# non-parametric model. As a result of this we will need to use a `modelimage` object to
# store cache information we will need to compute the numerical FT.
function (model::ImModel)(θ)
    (;c) = θ
    #Construct the image model
    cmat = reshape(c, model.npix, model.npix)
    img = IntensityMap(cmat, model.fov, model.fov, BSplinePulse{3}())
    #Create the modelimage object that will use a cache to compute the DFT
    m = modelimage(img, cache)
end
const npix = 10
const fovxy = μas2rad(65.0)
# Now we can feed in the array information to form the cache. We will be using a DFT since
# it is efficient for so few pixels
cache = create_cache(Comrade.DFTAlg(dlcamp), IntensityMap(rand(npix,npix), fovxy, fovxy, BSplinePulse{3}()))
mms = ImModel(cache, fovxy, npix)
# We will use a Dirichlet prior to enforce that the flux sums to unity since closures are
# degenerate to total flux.
prior = (c = DistributionsAD.TuringDirichlet(npix*npix, 1.0),)

post = Posterior(lklhd, prior, mms)
tpost = asflat(post)

# Let's run an optimizer to get a nice starting location
# It turns out that gradients are really helpful here
ndim = dimension(tpost)
f = OptimizationFunction(tpost, GalacticOptim.AutoZygote())
x0 = randn(ndim)
prob = OptimizationProblem(f, x0, nothing)
sol = solve(prob, LBFGS(); maxiters=4_000, show_trace=true, show_every=50, g_tol=1e-2)

xopt = transform(tpost, sol)

# Let's see how the fit looks
residual(mms(xopt), dlcamp)
residual(mms(xopt), dcphase)
plot(mms(xopt), fovx=1.2*fovxy, fovy=1.2*fovxy)


# now we sample using hmc
metric = DenseEuclideanMetric(ndim)
hchain, stats = sample(post, AHMC(;metric, autodiff=AD.ZygoteBackend()), 3000; nadapts=2000, init_params=xopt)


foox = let dc=dc, dd=dd, pr=tpost.lpost.prior, tr=tpost.transform, plan = cache.plan, phases=cache.phases, dmat=dlcamp.config.designmat, dmatc=dcphase.config.designmat
    x->begin
        y = transform(tr, x)
        lp = logdensityof(pr, y)
        vis = conj.(plan*y.c).*phases
        logdensityof(dd, dmat*log.(abs.(vis))) + logdensityof(dc, dmatc*angle.(vis)) + lp
    end
end
