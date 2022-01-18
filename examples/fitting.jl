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


load_ehtim()
obs = ehtim.obsdata.load_uvfits("/home/ptiede/Research/ngEHT/ngEHTChallenge2/data/Challenge_2/synthetic_data/M87_GRMHD/eht2022/230/M87_GRMHD_230_122_eht2022_04_allnoise_10s_230GHz_seed-1_netcal.uvfits")
obs.add_scans()
obsavg = obs.avg_coherent(0.0, scan_avg=true)
obsavg.add_fractional_noise(0.02)
dlcamp = extract_lcamp(obsavg; count="min")
dcphase = extract_cphase(obsavg; count="min-cut0bl")

lklhd = RadioLikelihood(dlcamp, dcphase)

function make_model_mring(θ; nx=6, ny=6)
    (;r, αamp, αpha, σ, f, rg, τg, ξg, xg, yg) = θ
    α = αamp.*cos.(αpha)
    β = αamp.*sin.(αpha)
    xg = rg/sqrt(1-τg)
    yg = rg*sqrt(1-τg)
    #m = RImage(reshape(coeff, 6,6), SqExpPulse(ϵ))
    mg = shifted(rotated(stretched(Gaussian(), xg,yg), ξg), xg, yg)
    m = convolved(stretched(MRing{length(αpha)}(α, β), r, r), stretched(Gaussian(), σ, σ))
    return m*f + mg*(1-f)
end

prior_mring = (
         r = Dists.Uniform(μas2rad(5.0), μas2rad(40.0)),
         σ = Dists.Uniform(μas2rad(1.0), μas2rad(40.0)),
         αamp = For(1:5) do _ Dists.Uniform(0.0, 0.5) end,
         αpha = For(1:5) do _ Dists.Uniform(-1π, 1π) end,
         rg = Dists.Uniform(μas2rad(20.0), μas2rad(200.0)),
         xg = Dists.Uniform(μas2rad(-500.0), μas2rad(500.0)),
         yg = Dists.Uniform(μas2rad(-500.0), μas2rad(500.0)),
         τg = Dists.Uniform(0.0, 0.9),
         ξg = Dists.Uniform(-π/2, π/2),
         f = Dists.Uniform(0.0, 1.0)

        )


function make_model_rimg(θ; r=μas2rad(50.0), nx=6, ny=6)
    (;ϵ, coeff, f, rg, τg, ξg, xg, yg) = θ
    m = stretched(RImage(reshape(coeff, ny,nx), SqExpPulse(ϵ)), r, r)
    xg = rg/sqrt(1-τg)
    yg = rg*sqrt(1-τg)
    mg = shifted(rotated(stretched(Gaussian(), xg,yg), ξg), xg, yg)
    return m*f + mg*(1-f)
end

function make_model_trimg(θ; rimg=μas2rad(80.0), nx=8, ny=8)
    (;ϵ, coeff, r, αamp, αpha, σ, f, x0, y0) = θ
    m = shifted(stretched(RImage(reshape(coeff, ny,nx), SqExpPulse(ϵ)), rimg, rimg), x0, y0)
    α = αamp.*cos.(αpha)
    β = αamp.*sin.(αpha)
    #mr = MRing(α, β)
    mr = convolved(stretched(MRing(α, β), r, r), stretched(Gaussian(), σ, σ))
    return m*f + mr*(1-f)
end

nx = ny = 8
prior_trimg = (
         ϵ = Dists.Uniform(0.2, 4.0),
         coeff = Dists.Dirichlet(nx*ny, 1.0),
         r = Dists.Uniform(μas2rad(5.0), μas2rad(40.0)),
         σ = Dists.Uniform(μas2rad(0.1), μas2rad(1.0)),
         αamp = ntuple(x->Dists.Uniform(0.0 ,0.5), 2),
         αpha = ntuple(x->Dists.Uniform(-π, π), 2),
         f = Dists.Uniform(0.0, 1.0),
         x0 = Dists.Uniform(-μas2rad(40.0), μas2rad(40.0)),
         y0 = Dists.Uniform(-μas2rad(40.0), μas2rad(40.0))
        )
prior_rimg = (
            ϵ = Dists.Uniform(1.0/3, 3.0),
            coeff = Dists.Dirichlet(nx*ny, 1.0),
            rg = Dists.Uniform(μas2rad(200.0), μas2rad(2000.0)),
            xg = Dists.Uniform(μas2rad(-500.0), μas2rad(500.0)),
            yg = Dists.Uniform(μas2rad(-500.0), μas2rad(500.0)),
            τg = Dists.Uniform(0.0, 0.9),
            ξg = Dists.Uniform(-π/2, π/2),
            f = Dists.Uniform(0.0, 1.0)

           )


tr = asflat(prior_trimg)
tt = asflat(prior_rimg)
tc = ascube(prior_mring)
nd = dimension(tr)

function make_lpost(lklhd, prior, make_model)
    tr = asflat(prior)
    function lpost(x)
    y, logjac = transform_and_logjac(tr, x)
    return logdensity(lklhd, make_model(y))  +
           logjac
    end
end

lp = make_lpost(lklhd, prior_trimg, make_model_trimg)
lp(rand(nd))
∇lpost(x) = ForwardDiff.gradient(lpost, x)
∇lpost(rand(nd))
@time ∇lpost(rand(nd));
f = OptimizationFunction((x,p)->-lpost(x), GalacticOptim.AutoForwardDiff())
res = bboptimize(x->-lpost(x);
                 SearchRange=(-7.0, 7.0),
                 NumDimensions=nd,
                 MaxFuncEvals=50_000
                )
@show best_fitness(res)
prob = GalacticOptim.OptimizationProblem(f, best_candidate(res)+randn(nd)*0.5, nothing)
prob = GalacticOptim.OptimizationProblem(f, sol.u+randn(nd)*0.5, nothing)
sol = solve(prob, LBFGS())
@show sol.minimum
xopt = transform(tr, sol.u)


mopt = make_model_trimg(xopt; nx, ny)
img = intensitymap(mopt, μas2rad(200.0), μas2rad(200.0), 1024, 1024)

x,y = imagepixels(img)
heatmap(rad2μas(x), rad2μas(y), img, xflip=true, aspect_ratio=:equal, size=(600,500))
