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
using Dynesty
using AdaptiveMCMC

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
    (;c, f, r, σ, τ, ξ, ma, mp, x, y) = θ
    #cmat = f*reshape(c, model.npix, model.npix)
    #img = IntensityMap(cmat, model.fov, model.fov, BSplinePulse{3}())
    #mimg = modelimage(img, model.cache)
    s,c = sincos(mp)
    α = ma*c
    β = ma*s
    #ring = (1-f)*rotated(smoothed(stretched(MRing((α,), (β,)), r, r*τ),σ), ξ)
    return MRing((α,), (β,))#stretched(MRing((α,), (β,)), r, r)#shifted(ring, x, y) #+ mimg
end

fovxy = μas2rad(90.0)
npix = 7
prior = (
          c = Dirichlet(npix*npix, 0.5),
          f = Uniform(0.0, 1.0),
          r = Uniform(μas2rad(10.0), μas2rad(30.0)),
          σ = truncated(Normal(0.0, μas2rad(0.5)), μas2rad(0.01), μas2rad(20.0)),
          τ = Uniform(0.5, 1.0),
          ξ = Uniform(-π/2, π/2),
          ma = Uniform(0.0, 0.5),
          mp = Uniform(-π, π),
          x = Uniform(-fovxy/2, fovxy/2),
          y = Uniform(-fovxy/2, fovxy/2)
        )
# Now form the posterior
cache = create_cache(Comrade.DFTAlg(dlcamp), IntensityMap(rand(npix,npix), fovxy, fovxy, BSplinePulse{3}()))
mms = Model(cache, fovxy, npix)

post = Posterior(lklhd, prior, mms)
tpost = asflat(post)
# We will use HMC to sample the posterior.
# First to reduce burn in we use pathfinder
ndim = dimension(tpost)
f, tr = GalacticOptim.OptimizationFunction(post, GalacticOptim.AutoForwardDiff{25}())
x0 = inverse(tr, rand(post.prior))
prob = GalacticOptim.OptimizationProblem(f, x0, nothing, lb=fill(-5.0, dimension(tr)), ub=fill(5.0, dimension(tr)))
sol, xopt = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), tr; maxiters=1000_000)

prob = GalacticOptim.OptimizationProblem(f, 3*randn(ndim), nothing)
sol, xopt = solve(prob, LBFGS(), tr; show_trace=true, iterations=4000, g_tol=1e-2)



start = [transform(tr, 3*randn(ndim)) for _ in 1:50]
q, ϕ, compinds = multipathfinder(post, 50; init_params=start, importance=true, iterations=2000)

residual(mms(xopt), dlcamp)
residual(mms(xopt), dcphase)
img = intensitymap(mms(xopt), μas2rad(160.0), μas2rad(160.0), 512, 512)
plot(abs.(img) , xlims=(-60.0,60.0), ylims=(-60.0,60.0), colorbar_scale=:log10, clims=(maximum(img)/100, maximum(img)),)

# now we sample using hmc
metric = DenseEuclideanMetric(ndim)
chain, stats = sample(post, HMC(;metric, autodiff=Comrade.AD.ForwardDiffBackend(;chunksize=Val(25))), 10000; nadapts=7000, init_params=xopt)

serialize("m87hybrid_fits.jls", Dict(:chain=>chain, :stats=>stats, :fov=>fovxy, :npix=>npix, :prior=>priors))
