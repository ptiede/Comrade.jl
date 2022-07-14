using Pkg;Pkg.activate(@__DIR__)
using Comrade
using Distributions
using ComradeAHMC
using ComradeOptimization
using OptimizationBBO
using Plots
# load eht-imaging we use this to load eht data
load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
file = "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, file))
obs.add_scans()
# kill 0-baselines since we don't care about
# large scale structure and make scan-average data
obs = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true)
# extract log closure amplitudes and closure phases
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)
# form the likelihood
lklhd = RadioLikelihood(dlcamp, dcphase)
# build the model: here we fit a ring with a azimuthal
# brightness variation and a Gaussian
function model(params)
  (;rad, wid, a, b, f, sig, asy, pa, x, y) = params
  ring = f*smoothed(stretched(MRing((a,), (b,)), rad, rad), wid)
  g = (1-f)*shifted(rotated(stretched(Gaussian(), sig*asy, sig), pa), x, y)
  return ring + g
end
# define the priors
uas2rad = pi/180.0/3600/1e6
prior = (
          rad = Uniform(uas2rad*(10.0), uas2rad*(30.0)),
          wid = Uniform(uas2rad*(1.0), uas2rad*(10.0)),
          a = Uniform(-0.5, 0.5), b = Uniform(-0.5, 0.5),
          f = Uniform(0.0, 1.0),
          sig = Uniform(uas2rad*(1.0), uas2rad*(40.0)),
          asy = Uniform(0.0, 0.75),
          pa = Uniform(0.0, 1pi),
          x = Uniform(-uas2rad*(80.0), uas2rad*(80.0)),
          y = Uniform(-uas2rad*(80.0), uas2rad*(80.0))
        )
# Now form the posterior
post = Posterior(lklhd, prior, model)
# We will use HMC to sample the posterior.
# First we will find a reasonable starting location using Optimization
# to have nice bounds we first transform to the unit hypercube
tpost = ascube(post)
ndim = dimension(tpost)
f = OptimizationFunction(tpost)
prob = OptimizationProblem(
            f, rand(ndim), nothing;
            lb=fill(1e-2, ndim), ub = fill(0.99, ndim)
            )
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=50_000)
# transform the solution back to regular space
xopt = transform(tpost, sol.u)
# Comrade is all about uncertainty quantification so now let's find the posterior!
# To do this we will use the `AdvancedHMC` package or rather its interface to Comrade.
metric = DiagEuclideanMetric(ndim)
chain, stats = sample(post, AHMC(;metric), 3000; nadapts=2000, init_params=xopt)
