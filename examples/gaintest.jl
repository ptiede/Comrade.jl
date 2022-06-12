using Pkg;Pkg.activate(@__DIR__)
using Comrade
using PyCall
using Plots

using ComradeAHMC

load_ehtim()
@pyimport eht_dmc as ed

# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obs.add_scans()

obsavg = obs.avg_coherent(0.0, scan_avg=true)

img = ehtim.image.make_empty(512, μas2rad(250.0), obs.ra, obs.dec, obs.rf, obs.source)
img = img.add_gauss(1.0, [μas2rad(40.0), μas2rad(30.0), π/3, 0.0, 0.0])

obslist = obsavg.split_obs()

gainoff = 0.0
gainp = Dict("AA"=> 0.01,
             "AP" => 0.1,
             "AZ" => 0.1,
             "JC" => 0.1,
             "LM" => 0.2,
             "PV" => 0.1,
             "SM" => 0.1,
             "SR" => 0.1)
obsim = img.observe_same(ehtim.obsdata.merge_obs(obslist[1:15]), ttype="fast",  gainp=gainp, gain_offset=gainoff, ampcal=false, phasecal=false, stabilize_scan_phase=true, stabilize_scan_amp=true)
# now get cal table
ctable = ehtim.calibrating.self_cal.self_cal(obsim, img, caltable=true)
ehtim.caltable.save_caltable(ctable, obsim, datadir=joinpath(@__DIR__ ,"CaltableTest"))


dvis = extract_vis(obsim)
dlcamp = extract_lcamp(obsim)
dcphase = extract_cphase(obsim)
damp = extract_amp(obsim)
st = scantable(dvis)

gcache = Comrade.GainCache(st)

struct Test{G}
    gcache::G
end

function Test(st::Comrade.ScanTable)
    gcache = Comrade.GainCache(st)
    return Test{typeof(gcache)}(gcache)
end

function (mod::Test)(θ)
    (;f, σ, τ, ξ,  gamp, gphase) = θ
    g = exp.(gamp).*cis.(gphase)
    m = f*Comrade.rotate(stretch(Gaussian(), σ*τ, σ), ξ)
    Comrade.GainModel(mod.gcache, g, m)
end

# define the priors
distamp = (AA = Normal(0.0, 0.1),
           AP = Normal(0.0, 0.1),
           LM = Normal(0.0, 0.2),
           AZ = Normal(0.0, 0.1),
           JC = Normal(0.0, 0.1),
           PV = Normal(0.0, 0.1),
           #SM = LogNormal(0.0, 0.1)
         )
distphase = (AA = Normal(0.0, 1e-4),
             AP = Normal(0.0, π),
             LM = Normal(0.0, π),
             AZ = Normal(0.0, 1.0*π),
             JC = Normal(0.0, 1.0*π),
             PV = Normal(0.0, π),
             #SM = Normal(0.0, 1.0*π),
             )
prior = (
          f = Uniform(0.9, 1.1),
          σ = Uniform(μas2rad(10.0), μas2rad(30.0)),
          τ = Uniform(0.1, 1.0),
          ξ = Uniform(-π/2, π/2),
          gamp = Comrade.GainPrior(distamp, st),
          gphase = Comrade.GainPrior(distphase, st)
        )
# Now form the posterior

mms  = Test(st)

cllklhd = RadioLikelihood(dlcamp, dcphase)
clpost = Posterior(cllklhd, prior, mms)

lklhd = RadioLikelihood(dvis)#damp, dcphase)
post = Posterior(lklhd, prior, mms)
tpost = asflat(post)
# We will use HMC to sample the posterior.
# First to reduce burn in we use pathfinder

ndim = dimension(tpost)
f = OptimizationFunction(tpost, GalacticOptim.AutoForwardDiff())
x0 = Comrade.HypercubeTransform.inverse(tpost, rand(post.prior))

prob = OptimizationProblem(f, x0, nothing, lb=fill(-5.0, ndim), ub=fill(5.0, ndim))
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1_000_000)

prob = GalacticOptim.OptimizationProblem(f, 0.1*randn(ndim), nothing)
sol = solve(prob, LBFGS(); g_tol=1e-2, maxiters=2_000)

@info sol.minimum

xopt = transform(tpost, sol)

residual(mms(xopt), dvis)


metric = DenseEuclideanMetric(ndim)
chain, stats = sample(post, AHMC(;metric), 7_000; nadapts=5_000, init_params=xopt)
