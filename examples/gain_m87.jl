using Pkg; Pkg.activate(@__DIR__)
using Comrade
using PyCall
using Plots

using GalacticOptim, Optim
using Zygote
using Parameters
using Distributions, DistributionsAD
using AdvancedHMC
using AbstractDifferentiation
using ReverseDiff
using ForwardDiff

load_ehtim()
@pyimport eht_dmc as ed

# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obs.add_scans()

obsavg = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true).add_fractional_noise(0.01)

obslist = obsavg.split_obs()

dvis = extract_vis(ehtim.obsdata.merge_obs(obslist[1:20]))
st = scantable(dvis)

gcache = Comrade.GainCache(st)

function model(θ, cache, gcache, fovx, fovy, ny, nx)
    @unpack c, f, gamp, gphase = θ
    g = gamp.*cis.(gphase)
    cmat = reshape(c, ny, nx)
    img = IntensityMap(cmat, fovx, fovy, BSplinePulse{3}())
    m = f*modelimage(img, cache)
    Comrade.GainModel(gcache, g, m)
end
# define the priors
distamp = (AA = LogNormal(0.0, 0.01),
         AP = LogNormal(0.0, 0.1),
         LM = LogNormal(0.0, 0.5),
         AZ = LogNormal(0.0, 0.1),
         JC = LogNormal(0.0, 0.1),
         PV = LogNormal(0.0, 0.1),
         SM = LogNormal(0.0, 0.1)
         )
distphase = (AA = VonMises(0.0, 1e4),
             AP = VonMises(0.0, 1/π^2),
             LM = VonMises(0.0, 1/π^2),
             AZ = VonMises(0.0, 1/π^2),
             JC = VonMises(0.0, 1/π^2),
             PV = VonMises(0.0, 1/π^2),
             SM = VonMises(0.0, 1/π^2)
             )
# Now form the posterior

prior = (
          c = MvLogNormal(20*20, 1.0),
          f = Uniform(0.2, 0.9),
          gamp = Comrade.GainPrior(distamp, st),
          gphase = Comrade.GainPrior(distphase, st)
        )
# Now form the posterior


gcache = Comrade.GainCache(st)
cache = create_cache(Comrade.DFTAlg(dvis[:u], dvis[:v]), IntensityMap(rand(20,20), μas2rad(100.0), μas2rad(100.0), BSplinePulse{3}()))
mms = let cache=cache, gcache=gcache
    θ -> model(θ, cache, gcache, μas2rad(100.0), μas2rad(100.0), 20, 20)
end

lklhd = RadioLikelihood(dvis)
post = Posterior(lklhd, prior, mms)
tpost = asflat(post)
# We will use HMC to sample the posterior.
# First to reduce burn in we use pathfinder
ndim = dimension(tpost)
function make_gradient(ℓ)
    ftape = ReverseDiff.GradientTape(ℓ, rand(ndim))
    compile_ftape = ReverseDiff.compile(ftape)
    out = similar(rand(ndim))

    (y,x,p)->ReverseDiff.gradient!(y, compile_ftape, x)
end
f, tr = GalacticOptim.OptimizationFunction(post, grad=grad)
x0 = inverse(tr, rand(post.prior))
prob = GalacticOptim.OptimizationProblem(f, x0, nothing, lb=fill(-5.0, dimension(tr)), ub=fill(5.0, dimension(tr)))
sol, xopt = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), tr; maxiters=1000_000)

res = map(1:50) do _
    i = rand(1:5)
    #prob = GalacticOptim.OptimizationProblem(f, inverse(tr, xopts[i]) .+ 0.1*randn(ndim), nothing)
    #sol, xopt = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), tr; maxiters=1000_000)

    #xopt = rand(post.prior)
    prob = GalacticOptim.OptimizationProblem(f, 0.25*randn(ndim), nothing)
    sol, xopt = solve(prob, LBFGS(), tr, show_trace=true, g_tol=1e-2)
    @info -sol.minimum
    xopt, -sol.minimum
end

xopt = res[end][1]
residual(mms(xopt), dvis)

q, phi, _ = multipathfinder(post, 100)
# now we sample using hmc
metric = DenseEuclideanMetric(ndim)
chain, stats = sample(post, HMC(;metric, autodiff=AD.ReverseDiffBackend()), 2000; nadapts=1000, init_params=transform(tpost,randn(ndim)*0.05))
