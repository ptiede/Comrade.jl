using Pkg;Pkg.activate(@__DIR__)
using Comrade
using Plots
using ComradeAHMC
using ComradeOptimization
using OptimizationBBO
using Distributions
using DistributionsAD


# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obs.add_scans()
obsavg = obs.avg_coherent(0.0, scan_avg=true)

# Lets make a gaussian image
img = ehtim.image.make_empty(512, μas2rad(250.0), obs.ra, obs.dec, obs.rf, obs.source)
img = img.add_gauss(1.0, [μas2rad(40.0), μas2rad(30.0), π/3, 0.0, 0.0])

# Set the gain priors
gainoff = 0.0
gainp = Dict("AA"=> 0.01,
             "AP" => 0.1,
             "AZ" => 0.1,
             "JC" => 0.1,
             "LM" => 0.2,
             "PV" => 0.1,
             "SM" => 0.1,
             "SP" => 0.1)
obsim = img.observe_same(obsavg, ttype="fast",  gainp=gainp, gain_offset=gainoff, ampcal=false, phasecal=false, stabilize_scan_phase=true, stabilize_scan_amp=true)
# now get cal table
ctable = ehtim.calibrating.self_cal.self_cal(obsim, img, caltable=true)
#ehtim.caltable.save_caltable(ctable, obsim, datadir=joinpath(@__DIR__ ,"CaltableTest"))


dcphase = extract_cphase(obsim)
damp = extract_amp(obsim)
st = timetable(damp)

gcache = Comrade.GainCache(st)

# Now create the Comrade Gain model
struct Model{G}
    gcache::G
end

function Model(st::Comrade.TimeTable)
    gcache = Comrade.GainCache(st)
    return Model{typeof(gcache)}(gcache)
end

function (mod::Model)(θ)
    (;f, σ, τ, ξ,  gamp) = θ
    g = exp.(gamp)
    m = f*rotated(stretched(Gaussian(), σ*τ, σ), ξ)
    return GainModel(mod.gcache, g, m)
end

# define the priors
distamp = (AA = Normal(0.0, 0.01),
           AP = Normal(0.0, 0.1),
           LM = Normal(0.0, 0.2),
           AZ = Normal(0.0, 0.1),
           JC = Normal(0.0, 0.1),
           PV = Normal(0.0, 0.1),
           SM = Normal(0.0, 0.1)
         )
prior = (
          f = Uniform(0.9, 1.1),
          σ = Uniform(μas2rad(10.0), μas2rad(30.0)),
          τ = Uniform(0.1, 1.0),
          ξ = Uniform(-π/2, π/2),
          gamp = Comrade.GainPrior(distamp, st),
        )
# Now form the posterior

mms  = Model(st)

lklhd = RadioLikelihood(mms, damp, dcphase)
post = Posterior(lklhd, prior)

tpost = asflat(post)
# We will use HMC to sample the posterior.
# First to reduce burn in we use pathfinder

ndim = dimension(tpost)
f = OptimizationFunction(tpost, Optimization.AutoForwardDiff())
x0 = Comrade.HypercubeTransform.inverse(tpost, rand(post.prior))

using OptimizationOptimJL
prob = OptimizationProblem(f, rand(ndim) .- 0.5, nothing)
sol = solve(prob, LBFGS(); g_tol=1e-5, maxiters=2_000)

@info sol.minimum

xopt = Comrade.transform(tpost, sol)

residual(mms(xopt), damp)

plot(caltable(mms(xopt)), layout=(3,3), size=(600,500))

# Compare the results for LMT
ctab = caltable(mms(xopt))
scatter(ctab[:time], inv.(ctab[:LM]),
        label="Comrade", size=(400,300),
        xlabel="Time (hr)",
        ylabel="LMT Gain Amp.",
        )
scatter!(pyconvert(Vector, ctable.data["LM"]["time"]), abs.(pyconvert(Vector, ctable.data["LM"]["rscale"])), label="eht-imaging")
