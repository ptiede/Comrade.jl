using Pkg;Pkg.activate(@__DIR__)
using Comrade
using PyCall
using Plots

using GalacticOptim, Optim
using Zygote
using Parameters
using Distributions
using SpeedMapping
using AdvancedHMC
using Dynesty

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
obsim = img.observe_same(ehtim.obsdata.merge_obs(obslist[10:15]), ttype="fast",  gainp=gainp, gain_offset=gainoff, ampcal=false, phasecal=false, stabilize_scan_phase=true, stabilize_scan_amp=true)
# now get cal table
ctable = ehtim.calibrating.self_cal.self_cal(obsim, img, caltable=true)
ehtim.caltable.save_caltable(ctable, obsim, datadir=joinpath(@__DIR__ ,"CaltableTest"))


dvis = extract_vis(obsim)
dlcamp = extract_lcamp(obsim; count="min")
dcphase = extract_cphase(obsim, count="min-cut0bl")
damp = extract_amp(obsim)
st = scantable(dvis)

gcache = Comrade.GainCache(st)

struct Model{G}
    gcache::G
end

function Model(st::Comrade.ScanTable)
    gcache = Comrade.GainCache(st)
    return Model{typeof(gcache)}(gcache)
end

function (mod::Model)(θ)
    @unpack f, σ, τ, ξ,  gamp, gphase = θ
    g = gamp.*cis.(gphase)
    m = f*rotated(stretched(Gaussian(), σ*τ, σ), ξ)
    Comrade.GainModel(mod.gcache, g, m)
end

# define the priors
distamp = (AA = LogNormal(0.0, 0.1),
           AP = LogNormal(0.0, 0.1),
           LM = LogNormal(0.0, 0.2),
           AZ = LogNormal(0.0, 0.1),
           JC = LogNormal(0.0, 0.1),
           PV = LogNormal(0.0, 0.1),
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

mms  = Model(st)

cllklhd = RadioLikelihood(dlcamp, dcphase)
clpost = Posterior(cllklhd, prior, mms)

lklhd = RadioLikelihood(dvis)#damp, dcphase)
post = Posterior(lklhd, prior, mms)
tpost = asflat(post)
# We will use HMC to sample the posterior.
# First to reduce burn in we use pathfinder
ndim = dimension(tpost)
f, tr = GalacticOptim.OptimizationFunction(post, GalacticOptim.AutoForwardDiff{15}())
#x0 = inverse(tr, rand(post.prior))

x0 = rand(post.prior, 50)
q, ϕ, inds = multipathfinder(clpost, 100; init_params=x0, ndraws_per_run=50)

res = map(1:150) do i
    #x0 = xx
    #x0.gphase .= randn(length(xx.gphase))
    #x0.gamp .= 1 .+ 0.2*randn(length(xx.gamp))
    prob = GalacticOptim.OptimizationProblem(f, randn(ndim), nothing)
    #sol, xopt = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), tr; maxiters=200_000)

    #xopt = rand(post.prior)
    #prob = GalacticOptim.OptimizationProblem(f, sol.u, nothing)
    sol, xopt = solve(prob, LBFGS(), tr; iterations=2_000)
    @info "$i/50 $(sol.minimum)"
    xopt, sol.minimum
end



metric = DenseEuclideanMetric(ndim)
chain, stats = sample(post, HMC(;metric), 14_000; nadapts=10_000, init_params=res[ind[1]][1])
