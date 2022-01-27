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
using Parameters


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
    @unpack spin, inc, a, b, β, ξ = θ
    θg = μas2rad(3.364)
    r0 = 4.5
    spin = -0.94
    inc = deg2rad(17.0)
    χ = deg2rad(-135.0)
    ι = π/3
    g = Jk.Kerr(spin)
    #(;spin, inc, θg, ξ) = θ
    prf = Jk.GaussianRing(r0, a, b)
    v = Jk.FluidVelocity(β, χ)
    bmag = Jk.MagneticField(ι, χ+1π)
    T = typeof(spin)
    o = Jk.Observer(1, inc)
    acc = BAM(1, 1.0, prf, v, bmag)
    macc = SingleStokes(Jk.SimpleModel(acc, g, o), :I)
    mimg = modelimage(macc; fovx=40.0, fovy=40.0, nx=128, ny=128)
    fm = flux(mimg)
    return rotated(stretched(mimg, θg, θg), ξ)/fm
end


prior = (
         spin = Dists.Uniform(-0.98, 0),
         r0 = Dists.Uniform(2.0, 8.0),
         inc = Dists.Uniform(deg2rad(5.0), deg2rad(40.0)),
         a = Dists.Uniform(-4.0, 4.0),
         b = Dists.Uniform(1.0, 5.0),
         β = Dists.Uniform(0.0, 0.98),
         χ = Dists.Uniform(0.0, π),
         ι = Dists.Uniform(0.0, π),
         θg = Dists.Uniform( μas2rad(0.1), μas2rad(10.0)),
         ξ = Dists.Uniform(-π, π),
         #f = Dists.Uniform(0.0, 1.0)
        )

function make_lpost(lklhd, prior, make_model; chunk=5)
    tr = asflat(prior)
    function lpost(x)
        y, logjac = transform_and_logjac(tr, x)
        return logdensity(lklhd, make_model(y))  +
               logjac
    end
end


lp = make_lpost(lklhd, prior, make_model);
nd = dimension(lp.tr)
@time lp(rand(nd))

fopt = OptimizationFunction((x,p)->-lp(x), GalacticOptim.AutoForwardDiff{5}())
prob = GalacticOptim.OptimizationProblem(fopt, rand(nd), nothing, lb = fill(-5.0, nd), ub = fill(5.0, nd))
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters=200, maxtime=500.0)
@info "First min" sol.minimum
prob = GalacticOptim.OptimizationProblem(fopt, sol.u, nothing, lb = fill(-10.0, nd), ub = fill(10.0, nd))
sol = solve(prob, GCMAESOpt(); maxiters=500)
@show sol.minimum
xopt = transform(lp.tr, sol.u)


mopt = make_model(xopt)
dvis = extract_vis(obsavg)
ac = arrayconfig(dvis)
u,v = getuv(ac)
scatter(hypot.(u,v), abs.(visibility.(dvis.data)))
scatter!(hypot.(u,v), amplitude.(Ref(mopt), u, v)./4.0)
img = intensitymap(mopt, μas2rad(100.0), μas2rad(100.0), 128, 128)

x,y = imagepixels(img)
heatmap(rad2μas(x), rad2μas(y), img, xflip=true, aspect_ratio=:equal, size=(600,500))
savefig(joinpath(@__DIR__,"m87_bestfit.png"))

scatter(abs2.((ll.x .-  ll.k.ops.μ(mopt/20.0)).*ll.k.ops.κ(2)))
