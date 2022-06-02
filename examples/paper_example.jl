using Pkg; Pkg.activate(@__DIR__)
using Comrade
using Distributions
using ComradeGalactic
using ComradeAHMC
using GalacticBBO
using Plots
using StatsBase

# load eht-imaging we use this to load eht data
load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
obs.add_scans()
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obs = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true)
# extract log closure amplitudes and closure phases
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)
# form the likelihood
lklhd = RadioLikelihood(dlcamp, dcphase)
# build the model here we fit a ring with a azimuthal
#brightness variation and a Gaussian
function model(θ)
    (;rad, wid, a, b, f, sig, asy, pa, x, y) = θ
    ring = f*smoothed(stretched(MRing((a,), (b,)), rad, rad), wid)
    g = (1-f)*shifted(rotated(stretched(Gaussian(), sig*asy, sig), pa), x, y)
    return ring + g
end
# define the priors
prior = (
          rad = Uniform(μas2rad(10.0), μas2rad(30.0)),
          wid = Uniform(μas2rad(1.0), μas2rad(10.0)),
          a = Uniform(-0.5, 0.5), b = Uniform(-0.5, 0.5),
          f = Uniform(0.0, 1.0),
          sig = Uniform(μas2rad(1.0), μas2rad(40.0)),
          asy = Uniform(0.0, 0.9),
          pa = Uniform(0.0, 1π),
          x = Uniform(-μas2rad(80.0), μas2rad(80.0)),
          y = Uniform(-μas2rad(80.0), μas2rad(80.0))
        )
# Now form the posterior
post = Posterior(lklhd, prior, model)
# We will use HMC to sample the posterior.
# First we will find a reasonable starting location using GalacticOptim
# For optimization we need to specify what transform to use. Here we will transform to
# the unit hypercube
tpost = ascube(post)
ndim = dimension(tpost)
f = OptimizationFunction(tpost)
prob = OptimizationProblem(f, rand(ndim), nothing, lb=fill(0.01, ndim), ub = fill(0.99, ndim))
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=50_000)

# transform the solution back to regular space
xopt = transform(tpost, sol)
# Let's see if the best fit looks reasonable by plotting normalized residuals for the
# log-closure amplitudes
residual(model(xopt), dlcamp)
# we can also plot the best fit model or maximum likelihood estimate (MLE)
plot(model(xopt), xlims=(-80.0,80.0),ylims=(-80.0,80.0), colorbar=nothing, title="MLE M87")

# Comrade is all about uncertainty quantification so now let's find the posterior!
# To do this we will use the `AdvancedHMC` package or rather its interface to Comrade.
metric = DiagEuclideanMetric(ndim)
chain, stats = sample(post, AHMC(;metric), 2000; nadapts=1000, init_params=xopt)
# Now let's find the mean image
images = intensitymap.(model.(sample(chain, 200)), μas2rad(160.0), μas2rad(160.0), 256, 256)
plot(mean(images), xlims=(-80.0, 80.0), ylims=(-80.0,80.0), colorbar=nothing, title="Mean M87")
