using Pkg; Pkg.activate(@__DIR__)
using Comrade
using PyCall
using Plots

using Zygote
using Distributions, DistributionsAD
using ComradeAHMC
using ComradeOptimization
using GalacticBBO
using ForwardDiff
using GalacticOptimJL

load_ehtim()
@pyimport eht_dmc as ed

# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obs.add_scans()

obsavg = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true).add_fractional_noise(0.02)

obslist = obsavg.split_obs()
obss = ehtim.obsdata.merge_obs(obslist)
damp = extract_amp(obss)
dcp = extract_cphase(obss)
st = scantable(damp)

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
    (;c, gamp) = θ
    g = exp.(gamp)
    cmat = reshape(c, model.ny, model.nx)
    img = IntensityMap(cmat, model.fovx, model.fovy, BSplinePulse{3}())

    m = modelimage(img, model.cache)
    return GainModel(model.gcache, g, m)
end

# define the priors
distamp = (AA = Normal(0.0, 0.1),
         AP = Normal(0.0, 0.1),
         LM = Normal(0.0, 1.0),
         AZ = Normal(0.0, 0.1),
         JC = Normal(0.0, 0.1),
         PV = Normal(0.0, 0.1),
         SM = Normal(0.0, 0.1)
         )
distphase = (AA = Normal(0.0, 0.0001),
             AP = Normal(0.0, 1π),
             LM = Normal(0.0, 1π),
             AZ = Normal(0.0, 1π),
             JC = Normal(0.0, 1π),
             PV = Normal(0.0, 1π),
             SM = Normal(0.0, 1π)
             )
# Now form the posterior
npix = 12
prior = (
          c = MvLogNormal(fill(-log(4*npix*npix), npix*npix), fill(2.0, npix*npix)),
          #c = DistributionsAD.TuringDirichlet(npix*npix, 1.0),
          #f = Uniform(0.0, 1.0),
          gamp = Comrade.GainPrior(distamp, st),
        )

# Now form the posterior

fovxy = μas2rad(70.0)
cache = create_cache(Comrade.DFTAlg(damp), IntensityMap(zeros(npix,npix), fovxy, fovxy, BSplinePulse{3}()))
mms = Model(gcache, cache, fovxy, fovxy, npix, npix)

#mms = NoGModel(cache, fovxy, fovxy, npix, npix)
lklhd = RadioLikelihood(damp, dcp)
post = Posterior(lklhd, prior, mms)
tpost = asflat(post)
lp = logdensityof(tpost)

# We will use HMC to sample the posterior.
# First to reduce burn in we use pathfinder
ndim = dimension(tpost)
f = OptimizationFunction(tpost, GalacticOptim.AutoZygote())
x0 = Comrade.HypercubeTransform.inverse(tpost, prior_sample(post))

#prob = OptimizationProblem(f, zeros(ndim), nothing, lb=fill(-5.0, ndim), ub=fill(5.0, ndim))
#sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1_000_000)

prob = OptimizationProblem(f, 0.5*randn(ndim), nothing)
prob = OptimizationProblem(f, sol.u, nothing)
function cb(x, args...)
    @info "lpost: $(lp(x))"
    return false
end
sol = solve(prob, LBFGS(); g_tol=1e-2, maxiters=2_000, callback=cb)

xopt = transform(tpost, sol)

residual(mms(xopt), dcp)
residual(mms(xopt), damp)

plot(mms(xopt), xlims=(-40.0,40.0), ylims=(-40.0,40.0))

# now we sample using hm
metric = DenseEuclideanMetric(ndim)
chain, stats = sample(post, AHMC(;metric), 6000; nadapts=4000, init_params=xopt)
