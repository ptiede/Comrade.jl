using Pkg; Pkg.activate(@__DIR__)
using Comrade
using PyCall
using Plots

using Zygote
using Distributions, DistributionsAD
using ComradeAHMC
using ComradeGalactic
using GalacticBBO
using ForwardDiff
using GalacticOptimJL
using DistributionsAD

load_ehtim()
@pyimport eht_dmc as ed

# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obs.add_scans()

obsavg = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true).add_fractional_noise(0.02)

obslist = obsavg.split_obs()

dvis = extract_vis(ehtim.obsdata.merge_obs(obslist))
st = scantable(dvis)

gcache = GainCache(st)

struct Model{G,C}
    gcache::G
    cache::C
    fovx::Float64
    fovy::Float64
    nx::Int
    ny::Int
end

function (model::Model)(θ)
    (;c, gamp, gphase) = θ
    g = exp.(gamp).*cis.(gphase)
    cmat = reshape(c, model.ny, model.nx)
    img = IntensityMap(cmat, model.fovx, model.fovy, BSplinePulse{3}())

    m = modelimage(img, model.cache)
    return GainModel(model.gcache, g, m)
end

# define the priors
distamp = (AA = Normal(0.0, 1e-3),
         AP = Normal(0.0, 0.1),
         LM = Normal(0.0, 0.5),
         AZ = Normal(0.0, 0.1),
         JC = Normal(0.0, 0.1),
         PV = Normal(0.0, 0.1),
         SM = Normal(0.0, 0.1)
         )
distphase = (AA = Normal(0.0, 0.001),
             AP = Normal(0.0, 1π),
             LM = Normal(0.0, 1π),
             AZ = Normal(0.0, 1π),
             JC = Normal(0.0, 1π),
             PV = Normal(0.0, 1π),
             SM = Normal(0.0, 1π)
             )
# Now form the posterior
npix = 8
prior = (
          c = MvLogNormal(fill(-log(10.0*npix*npix), npix*npix), fill(1.0, npix*npix)),
          #c = DistributionsAD.TuringDirichlet(npix*npix, 1.0),
          #f = Uniform(0.0, 1.0),
          gamp = Comrade.GainPrior(distamp, st),
          gphase = Comrade.GainPrior(distphase, st)
        )

# Now form the posterior

fovxy = μas2rad(65.0)
gcache = GainCache(st)
cache = create_cache(Comrade.DFTAlg(dvis), IntensityMap(zeros(npix,npix), fovxy, fovxy, BSplinePulse{3}()))
mms = Model(gcache, cache, fovxy, fovxy, npix, npix)

#mms = NoGModel(cache, fovxy, fovxy, npix, npix)
lklhd = RadioLikelihood(dvis)
post = Posterior(lklhd, prior, mms)
tpost = asflat(post)
# We will use HMC to sample the posterior.
# First to reduce burn in we use pathfinder
ndim = dimension(tpost)
f = OptimizationFunction(tpost, GalacticOptim.AutoZygote())
x0 = Comrade.HypercubeTransform.inverse(tpost, prior_sample(post))

prob = OptimizationProblem(f, zeros(ndim), nothing, lb=fill(-5.0, ndim), ub=fill(5.0, ndim))
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1_000_000)

prob = GalacticOptim.OptimizationProblem(f, 1.5*randn(ndim), nothing)
prob = GalacticOptim.OptimizationProblem(f, sol.u.+0.5*randn(ndim), nothing)
lp = logdensityof(tpost)
function cb(x, args...)
    @info "lpost: $(lp(x))"
    return false
end
sol = solve(prob, LBFGS(); g_tol=1e-2, maxiters=3_000, callback=cb)

xopt = transform(tpost, sol)

residual(mms(xopt), dvis)

plot(mms(xopt), xlims=(-55.0,55.0), ylims=(-55.0,55.0))

# now we sample using hmc
metric = DenseEuclideanMetric(ndim)
chain, stats = sample(post, AHMC(;metric, autodiff=AD.ZygoteBackend()), 5000; nadapts=3000, init_params=xopt)
