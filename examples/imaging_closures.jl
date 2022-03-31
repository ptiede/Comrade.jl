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
using BlackBoxOptim

# load eht-imaging we use this to load eht data
load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obs.add_scans()
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obs = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true).add_fractional_noise(0.01)
# extract log closure amplitudes and closure phases
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)
lklhd = RadioLikelihood(dlcamp, dcphase)

# build the model here we fit a ring with a azimuthal
#brightness variation and a Gaussian

struct Model{C,F}
    cache::C
    fov::F
    npix::Int
end


function (model::Model)(θ)
    c  = θ.c
    ϵ = θ.ϵ
    cmat = reshape(c, model.npix, model.npix)
    img = IntensityMap(cmat, model.fov, model.fov, SqExpPulse(ϵ))
    return modelimage(img, model.cache)
end

fovxy = μas2rad(70.0)
npix = 9
prior = (
          c = Dirichlet(npix*npix, 1.0),
          ϵ = Uniform(0.5, 3.0)
        )
# Now form the posterior
cache = create_cache(Comrade.DFTAlg(dlcamp), IntensityMap(rand(npix,npix), fovxy, fovxy, SqExpPulse(1.0)))
mms = Model(cache, fovxy, npix)

post = Posterior(lklhd, prior, mms)
tpost = asflat(post)
# We will use HMC to sample the posterior.
# First to reduce burn in we use pathfinder
ndim = dimension(tpost)
f, tr = GalacticOptim.OptimizationFunction(post, GalacticOptim.AutoForwardDiff{25}())
x0 = inverse(tr, rand(post.prior))
prob = GalacticOptim.OptimizationProblem(f, x0, nothing, lb=fill(-15.0, dimension(tr)), ub=fill(15.0, dimension(tr)))
sol, xopt = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), tr; maxiters=1000_000)

prob = GalacticOptim.OptimizationProblem(f, 3*randn(ndim), nothing)
sol, xopt = solve(prob, LBFGS(), tr; show_trace=true, iterations=2000, g_tol=1e-2)

residual(mms(xopt), dlcamp)
residual(mms(xopt), dcphase)
plot(mms(xopt))
# now we sample using hmc
metric = DenseEuclideanMetric(ndim)
chain, stats = sample(post, HMC(;metric), 4500; nadapts=2500, init_params=xopt)

using Serialization
serialize(joinpath(@__DIR__, "m87_closure_fits.jls"), Dict(:chain => chain, :stats => stats, :model=>model, :post=>post))
