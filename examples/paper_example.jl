using Pkg; Pkg.activate(@__DIR__)
using Comrade
using Distributions
using ComradeOptimization
using ComradeAHMC
using OptimizationBBO
using OptimizationOptimJL
using Plots
using StatsBase

# load eht-imaging we use this to load eht data
load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obs = scan_average(obs.flag_uvdist(uv_min=0.1e9).add_fractional_noise(0.02))
# extract log closure amplitudes and closure phases
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)
# form the likelihood
# build the model here we fit a ring with a azimuthal
#brightness variation and a Gaussian

function model(θ)
    (;rad, wid, a, b, f, sig, asy, pa, x, y) = θ
    ring = f*smoothed(stretched(MRing((a,), (b,)), μas2rad(rad), μas2rad(rad)), μas2rad(wid))
    g = (1-f)*shifted(rotated(stretched(Gaussian(), μas2rad(sig)*asy, μas2rad(sig)), pa), μas2rad(x), μas2rad(y))
    return ring + g
end
# define the priors
prior = (
          rad = Uniform(10.0, 30.0),
          wid = Uniform(1.0, 10.0),
          a = Uniform(-0.5, 0.5), b = Uniform(-0.5, 0.5),
          f = Uniform(0.0, 1.0),
          sig = Uniform((1.0), (60.0)),
          asy = Uniform(0.25, 1.0),
          pa = Uniform(0.0, 1π),
          x = Uniform(-(80.0), (80.0)),
          y = Uniform(-(80.0), (80.0))
        )
# Now form the posterior
lklhd = RadioLikelihood(model,dlcamp, dcphase)
post = Posterior(lklhd, prior)
# We will use HMC to sample the posterior.
# First we will find a reasonable starting location using GalacticOptim
# For optimization we need to specify what transform to use. Here we will transform to
# the unit hypercube
tpost = asflat(post)
ndim = dimension(tpost)
f = OptimizationFunction(tpost, Optimization.AutoForwardDiff())
prob = OptimizationProblem(f, rand(ndim), nothing, lb=fill(-5.0, ndim), ub = fill(5.0, ndim))
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=50_000)

# Now let's get the Laplace approximation since it is cheap!
prob = OptimizationProblem(f, sol.u, nothing)
ldist = laplace(prob, LBFGS(); show_trace=true)

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
chain, stats = sample(post, AHMC(;metric), 4000; nadapts=2000, init_params=xopt)
# chain has the MCMC chain and stats includes ancilliary information
# Now we should check that the chain acutally mixed well. To do that we can compute the ESS
using MCMCDiagnostics
using Tables
ess = map(effective_sample_size, Tables.columns(chain))
# We can also calculate the split-rhat or potential scale reduction. For this we should actually
# use at least 4 chains. However for demonstation purposes we will use one chain that we split in two
rhats = map(Tables.columns(chain)) do c
    c1 = @view c[2001:3000]
    c2 = @view c[3001:4000]
    return potential_scale_reduction(c1, c2)
end
# Ok we have a split-rhat < 1.01 on all parameters so we have success (in reality run more chains!).

# Now let's find the mean image
images = intensitymap.(model.(sample(chain, 200)), μas2rad(160.0), μas2rad(160.0), 256, 256)
plot(mean(images), xlims=(-80.0, 80.0), ylims=(-80.0,80.0), colorbar=nothing, title="Mean M87")

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
