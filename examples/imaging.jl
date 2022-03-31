using Pkg; Pkg.activate(@__DIR__)
using Comrade
using Distributions
using Pathfinder
using AdvancedHMC
using Plots
using Zygote
using Parameters
using GalacticOptim
using Optim
using NLopt

# load eht-imaging we use this to load eht data
load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obs.add_scans()
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obs = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true).add_fractional_noise(0.02)
# extract log closure amplitudes and closure phases
vis = extract_vis(obs)
lklhd = RadioLikelihood(vis)

# Build the Model. Here we we a struct to hold some caches
# This will be useful to hold precomputed caches

struct Model{C,G}
    cache::C
    gcache::G
    fov::Float64
    npix::Int
end

function (model::Model)(θ)
    ;c, f, gamp, gphase = θ
    # Construct the image model
    cmat = f.*reshape(c, model.npix, model.npix)
    img = IntensityMap(cmat, model.fov, model.fov, BSplinePulse{3}())
    m = modelimage(img, cache)
    # Now corrupt the model with Gains
    g = gamp.*cis.(gphase)
    Comrade.GainModel(gcache, g, m)
end

# First we define the station gain priors
distamp = (AA = LogNormal(0.0, 0.1),
         AP = LogNormal(0.0, 0.1),
         LM = LogNormal(0.0, 0.5),
         AZ = LogNormal(0.0, 0.1),
         JC = LogNormal(0.0, 0.1),
         PV = LogNormal(0.0, 0.1),
         SM = LogNormal(0.0, 0.1))
# We will use AA as a reference array
distphase = (AA = VonMises(0.0, 10000.0),
             AP = VonMises(0.0, 1/π^2),
             LM = VonMises(0.0, 1/π^2),
             AZ = VonMises(0.0, 1/π^2),
             JC = VonMises(0.0, 1/π^2),
             PV = VonMises(0.0, 1/π^2),
             SM = VonMises(0.0, 1/π^2))

npix = 8
fov = μas2rad(70.0)
prior = (
          c = Dirichlet(npix*npix, 1.0),
          f = Uniform(0.2, 0.9),
          gamp = Comrade.GainPrior(distamp, scantable(vis)),
          gphase = Comrade.GainPrior(distphase, scantable(vis))
        )


# Construct cache's and design matrices
gcache = Comrade.GainCache(scantable(vis))
cache = create_cache(Comrade.DFTAlg(vis), IntensityMap(rand(npix,npix), μas2rad(fov), μas2rad(fov), BSplinePulse{3}()))

mms = Model(cache, gcache, fov, npix)

post = Posterior(lklhd, prior, mms)
tpost = asflat(post)
# We will use HMC to sample the posterior.
# First to reduce burn in we use pathfinder
ndim = dimension(tpost)
f, tr = GalacticOptim.OptimizationFunction(post, GalacticOptim.AutoReverseDiff())
x0 = inverse(tr, rand(post.prior))
prob = GalacticOptim.OptimizationProblem(f, x0, nothing, lb=fill(-15.0, dimension(tr)), ub=fill(15.0, dimension(tr)))
sol, xopt = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), tr; maxiters=1000_000)

prob = GalacticOptim.OptimizationProblem(f, sol.u, nothing)
sol, xopt = solve(prob, LBFGS(), tr)

residual(mms(xopt), vis)

q, phi, _ = multipathfinder(post, 100)
# now we sample using hmc
metric = DiagEuclideanMetric(ndim)
chain, stats = sample(post, HMC(;metric, autodiff=AD.ReverseDiffBackend()), 2000; nadapts=1000, init_params=xopt)
