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
using DistributionsAD

# load eht-imaging we use this to load eht data
load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obs.add_scans()
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obs = obs.avg_coherent(0.0, scan_avg=true).add_fractional_noise(0.01)
# extract log closure amplitudes and closure phases
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs; cut_trivial=true)
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
    cmat = reshape(c, model.npix, model.npix)
    img = IntensityMap(cmat, model.fov, model.fov, BSplinePulse{3}())
    return modelimage(img, model.cache) + stretched(Gaussian(), μas2rad(250.0), μas2rad(250.0))
end

fovxy = μas2rad(80.0)
npix = 15
prior = (
          #c = DistributionsAD.TuringDirichlet(npix*npix, 0.5),
          c = MvLogNormal(fill(log(inv(4*npix^2)), npix*npix), 3.0),
        )
# Now form the posterior
u, v = getuv(dlcamp.config)
cache = create_cache(Comrade.DFTAlg(dlcamp), IntensityMap(rand(npix,npix), fovxy, fovxy, BSplinePulse{3}()))
mms = Model(cache, fovxy, npix)

post = Posterior(lklhd, prior, mms)
tpost = asflat(post)
# We will use HMC to sample the posterior.
# First to reduce burn in we use pathfinder
ndim = dimension(tpost)
f, tr = GalacticOptim.OptimizationFunction(post, GalacticOptim.AutoZygote())
x0 = HypercubeTransform.inverse(tr, rand(post.prior))
prob = GalacticOptim.OptimizationProblem(f, x0, nothing, lb=fill(-5.0, dimension(tr)), ub=fill(5.0, dimension(tr)))
sol, xopt = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), tr; maxiters=1000_000)



prob = GalacticOptim.OptimizationProblem(f, randn(ndim), nothing)
sol, xopt = solve(prob, LBFGS(), tr; show_trace=true, show_every=50, iterations=4000, g_tol=1e-2)
@info sol.minimum

residual(mms(xopt), dlcamp)
residual(mms(xopt), dcphase)
plot(mms(xopt), fovx=1.2*fovxy, fovy=1.2*fovxy)

using AdaptiveMCMC
cacc, sacc = sample(post, Comrade.AdaptMCMC(ntemp=5), 500_000, 250_000; init_params=xopt, thin=50)

# now we sample using hmc
metric = DenseEuclideanMetric(ndim)
hchain, stats = sample(post, HMC(;metric, autodiff=AD.ZygoteBackend()), 7000; nadapts=5000, init_params=xopt)

using Serialization
serialize(joinpath(@__DIR__, "m87_closure_fits.jls"), Dict(:chain => chain, :stats => stats, :model=>model, :post=>post))

foox = let dc=dc, dd=dd, pr=tpost.lpost.prior, tr=tpost.transform, plan = cache.plan, phases=cache.phases, dmat=dlcamp.config.designmat, dmatc=dcphase.config.designmat
    x->begin
        y = transform(tr, x)
        lp = logdensityof(pr, y)
        vis = conj.(plan*y.c).*phases
        logdensityof(dd, dmat*log.(abs.(vis))) + logdensityof(dc, dmatc*angle.(vis)) + lp
    end
end
