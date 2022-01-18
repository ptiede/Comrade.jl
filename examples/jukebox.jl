using Pkg; Pkg.activate(@__DIR__)
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
using JuKeBOX
const Jk = JuKeBOX
using BlackBoxOptim
using GCMAES


load_ehtim()
obs = ehtim.obsdata.load_uvfits("/home/ptiede/Research/ngEHT/ngEHTChallenge2/data/Challenge_2/synthetic_data/M87_GRMHD/eht2022/230/M87_GRMHD_230_122_eht2022_04_allnoise_10s_230GHz_seed-1_netcal.uvfits")
obs.add_scans()
obsflg = obs.flag_uvdist(uv_min=0.9e9)
obsavg = obsflg.avg_coherent(0.0, scan_avg=true)
obsnn = obsavg.add_fractional_noise(0.02)
dlcamp = extract_lcamp(obsnn; count="min")
dcphase = extract_cphase(obsnn; count="min-cut0bl")

lklhd = RadioLikelihood(dlcamp, dcphase)


function make_model(θ)
    (;spin, inc, rpeak, width, β, χ, ι, θg, ξ) = θ
    #(;spin, inc, θg, ξ) = θ
    o = Jk.Observer(1, inc)
    g, acc = Jk.bam(1, spin, 1.0, rpeak, width, β, χ, ι)
    macc = SingleStokes(Jk.SimpleModel(acc, g, o), :I)
    mimg = modelimage(macc; fovx=30.0, fovy=30.0, nx=64, ny=64)
    #fm = flux(mimg)
    return rotated(stretched(mimg, θg, θg), ξ)
end

function test()
    params = transform(tf, rand(dimension(tf)))
    model = make_model(params)
    #heatmap(intensitymap(model, 150.0, 150.0, 256, 256))
end

prior = (
         spin = Dists.Uniform(-0.98, 0.98),
         inc = Dists.Uniform(0.01, π/2),
         rpeak = Dists.Uniform(2.1, 10.0),
         width = Dists.Uniform(1.0, 10.0),
         β = Dists.Uniform(0.0, 0.98),
         χ = Dists.Uniform(0.0, π/2),
         ι = Dists.Uniform(0.0, π/2),
         θg = Dists.Uniform( μas2rad(0.1), μas2rad(10.0)),
         ξ = Dists.Uniform(-π, π),
         #f = Dists.Uniform(0.0, 1.0)
        )

function make_lpost(lklhd, prior, make_model)
    tr = asflat(prior)
    function lpost(x)
        y, logjac = transform_and_logjac(tr, x)
        return logdensity(lklhd, make_model(y))  +
               logjac
    end
end


lp = make_lpost(lklhd, prior, make_model)
nd = dimension(lp.tr)
p = make_model(transform(lp.tr, rand(nd)))
lp(rand(nd))
@btime lp($rand(nd))
∇lpost(x) = ForwardDiff.gradient(lp, x)
∇lpost(rand(nd))
@time ∇lpost(rand(nd));
f = OptimizationFunction((x,p)->-lp(x))
prob = GalacticOptim.OptimizationProblem(f, rand(nd), nothing, lb = fill(0.01, nd), ub = fill(0.99, nd))
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters=10_000, maxtime=1000.0)
f = OptimizationFunction((x,p)->-lp(x), GalacticOptim.AutoForwardDiff())
prob = GalacticOptim.OptimizationProblem(f, sol.u, nothing, lb = fill(0.01, nd), ub = fill(0.99, nd))
sol = solve(prob, GCMAESOpt())
@show sol.minimum
xopt = transform(lp.tr, sol.u)


mopt = make_model(xopt)
dvis = extract_vis(obsavg)
ac = arrayconfig(dvis)
u,v = getuv(ac)
scatter(hypot.(u,v), abs.(visibility.(dvis.data)))
scatter!(hypot.(u,v), amplitude.(Ref(mopt), u, v))
img = intensitymap(mopt, μas2rad(100.0), μas2rad(100.0), 64, 64)

x,y = imagepixels(img)
heatmap(rad2μas(x), rad2μas(y), img, xflip=true, aspect_ratio=:equal, size=(600,500))
