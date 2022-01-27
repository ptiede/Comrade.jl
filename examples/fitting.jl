using Pkg; Pkg.activate(@__DIR__)
using BlackBoxOptim
using Comrade
using MeasureTheory
using ParameterHandling
using GalacticOptim
using AdvancedHMC
using AdvancedMH
import Distributions as Dists
using HypercubeTransform
using Optim
using Plots
using ForwardDiff
using Pathfinder
using SpeedMapping
using GCMAES
using Parameters
using Nonconvex
using NLopt



load_ehtim()
obs = ehtim.obsdata.load_uvfits("/run/media/ptiede/Samsung_T5/Research/Projects/JuKeBOX/Observations/eht2022/ma+0.5_r10_n1_tavg.uvfits")
obsavg = obs.avg_coherent(0.0, scan_avg=true)
obsnn = obsavg.add_fractional_noise(0.03)

damp = extract_amp(obsnn; debias=true)
dcphase = extract_cphase(obsnn; count="min-cut0bl")


function make_model_mring(θ)
    @unpack r, αamp, αpha, σ, f, rg, τg, ξg, xg, yg = θ
    α = αamp.*cos.(αpha)
    β = αamp.*sin.(αpha)
    xg = rg/sqrt(1-τg)
    yg = rg*sqrt(1-τg)
    #m = DImage(reshape(coeff, 6,6), SqExpPulse(ϵ))
    mg = shifted(rotated(stretched(Gaussian(), xg,yg), ξg), xg, yg)
    m = convolved(stretched(MRing(α, β), r, r), stretched(Gaussian(), σ, σ))
    return m*f + mg*(1-f)
end

prior_mring = (
         r = Dists.Uniform(μas2rad(5.0), μas2rad(40.0)),
         σ = Dists.Uniform(μas2rad(1.0), μas2rad(40.0)),
         αamp = (Dists.Uniform(0.0, 0.5), Dists.Uniform(0.0, 0.5)),
         αpha = (Dists.Uniform(-π, π), Dists.Uniform(-π, π)),
         rg = Dists.Uniform(μas2rad(20.0), μas2rad(200.0)),
         xg = Dists.Uniform(μas2rad(-50.0), μas2rad(50.0)),
         yg = Dists.Uniform(μas2rad(-50.0), μas2rad(50.0)),
         τg = Dists.Uniform(0.0, 0.9),
         ξg = Dists.Uniform(-π/2, π/2),
         f = Dists.Uniform(0.0, 1.0)
        )


function make_model_rimg(θ; r=μas2rad(60.0), nx=6, ny=6)
    @unpack coeff, f = θ
    m = stretched(DImage(reshape(coeff, ny,nx), BSplinePulse{3}()), r, r)
    return m*f
end

function make_model_trimg(θ; rimg=μas2rad(85.0), nx=7, ny=7)
    @unpack coeff, r, αamp, αpha, σ, f = θ
    m = stretched(DImage(reshape(coeff, ny,nx), BSplinePulse{3}()), rimg, rimg)
    α = αamp.*cos.(αpha)
    β = αamp.*sin.(αpha)
    #mr = MRing(α, β)
    mr = convolved(stretched(MRing(α, β), r, r), stretched(Gaussian(), σ, σ))
    return m*f + mr*(1-f)
end

nx = ny = 8
prior_trimg = (
         coeff = Dists.Dirichlet(nx*ny, 1.0),
         r = Dists.Uniform(μas2rad(5.0), μas2rad(40.0)),
         σ = Dists.Uniform(μas2rad(0.1), μas2rad(5.0)),
         αamp = ntuple(x->Dists.Uniform(0.0 ,0.5), 2),
         αpha = ntuple(x->Dists.Uniform(-π, π), 2),
         f = Dists.Uniform(0.0, 1.0),
        )
prior_rimg = (
            coeff = Dists.Dirichlet(nx*ny, 1.0),
            f = Dists.Uniform(0.0 ,3.0)
           )


function make_lpost(lklhd, prior, make_model)
    tr = asflat(prior)
    function lpost(x)
    y, logjac = transform_and_logjac(tr, x)
    return logdensity(lklhd, make_model(y))  +
           logjac
    end
end

function make_grad(lp, chunk)
    nd = dimension(lp.tr)
    out = zeros(nd)
    cfg = ForwardDiff.GradientConfig(lp, rand(nd), ForwardDiff.Chunk{chunk}())
    x->ForwardDiff.gradient!(out, lp, x, cfg)
end


lklhd = RadioLikelihood(damp, dcphase)

lp = make_lpost(lklhd, prior_mring, make_model_mring)
nd = dimension(lp.tr)
lp(rand(nd))

∂lp = make_grad(lp, 2)
∂lp(rand(nd))
@benchmark ∂lp($rand(nd))

using Pathfinder
θ₀s = collect.(eachcol(rand(nd, 50) .* 5 .- 2.5))
res = map(x->pathfinder(lp, ∂lp, rand(nd)*8.0, 10), 1:50);

xopts = transform.(Ref(lp.tr), eachcol(res[6][2]))

# fopt = OptimizationFunction((x,p)->-lp(x), GalacticOptim.AutoForwardDiff{25}())
# prob = GalacticOptim.OptimizationProblem(fopt, 2*randn(nd), nothing, lb = fill(-5.0, nd), ub = fill(5.0, nd))
# sol1 = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiter=200, maxtime=500.0)
# @info "First min" sol1.minimum

# prob = GalacticOptim.OptimizationProblem(fopt, 3*randn(nd), nothing, lb = fill(-10.0, nd), ub = fill(10.0, nd))
# sol = solve(prob, LBFGS())
# @show sol.minimum
# xopt = transform(lp.tr, sol1.u)

mopt = make_model_rimg(xopts[1]; nx, ny)
img = intensitymap(mopt, μas2rad(100.0), μas2rad(100.0), 1024, 1024)

x,y = imagepixels(img)
heatmap(rad2μas(x), rad2μas(y), img, xflip=true, aspect_ratio=:equal, size=(600,500))

metric = DiagEuclideanMetric(nd)
hamiltonian = Hamiltonian(metric, lp, x->(lp(x),∂lp(x)))
ϵ0 = find_good_stepsize(hamiltonian, ϕ[:,1])
integrator = Leapfrog(ϵ0)

proposal = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator)
adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator))

samples, stats = sample(hamiltonian, proposal, ϕ[:,1], 2000, adaptor, 1000; progress=true)


ac = arrayconfig(damp)
u,v = getuv(ac)
scatter(hypot.(u,v), (amplitude.(damp.data)))
scatter!(hypot.(u,v), amplitude.(Ref(mopt), u, v))
