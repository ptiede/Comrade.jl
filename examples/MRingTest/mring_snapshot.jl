using ROSE
using Plots
using StatsPlots
import ForwardDiff
#import ReverseDiff
import DiffResults
using Optim
using DelimitedFiles
using StructArrays
using StaticArrays
using Soss
using BlackBoxOptim
import ROSE: nsamples
using AdvancedHMC
using VIDA
using DataFrames
using CSV
using Statistics
using NLopt
using PyCall
@pyimport ehtim

include(joinpath(@__DIR__, "utility.jl"))

obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "mring.uvfits"))
obs = obs.avg_coherent(60.0)

function process_eht_obs(obs)
    obs_scan = obs.copy()
    obs_fit = obs_scan.copy()
    obs_fit = obs_fit.add_leakage_noise(min_noise=0.02)
    noise = get(obs_fit.data,"sigma")
    noise = sqrt.(noise.^2 .+ 0.01^2)
    set!(obs_fit.data, "sigma", noise)
    vis = get(obs_fit.data, "vis")
    zbl = mean(get(obs.flag_uvdist(uv_max=0.1e9).unpack("amp", debias=true), "amp"))
    set!(obs_fit.data, "vis", vis./zbl)
    set!(obs_fit.data, "sigma", vis./zbl)
    obs_fit.add_cphase()
    obs_fit.add_amp()
    return obs_fit
end


#################################################################3
# Imaging model
#################################################################

model = Soss.@model uvamp, erramp, uvcphase, errcphase begin

    rad ~ truncated(Normal(25.0, 15.0), 0.0, Inf)
    σ ~ truncated(Normal(20.0, 20.0), 0.0, Inf)
    α1 ~ Uniform(-0.5, 0.5)
    α2 ~ Uniform(-0.5, 0.5)
    β1 ~ Uniform(-0.5, 0.5)
    β2 ~ Uniform(-0.5, 0.5)
    #f ~ Normal(0.1, 10.0)

    img = MRing(rad, SVector(α1, α2), SVector(β1, β2)) |>
                smoothed(σ) |>
                renormed(f)



    amp ~ For(eachindex(uvamp, erramp)) do i
        u,v = uvamp[i]
        mamp = ROSE.visibility_amplitude(img, u, v)
        Normal(mamp, erramp[i])
    end

    cphase ~ For(eachindex(uvcphase, errcphase)) do i
        u1,v1,u2,v2,u3,v3 = uvcphase[i]
        mphase = ROSE.closure_phase(mamp, u1, v1, u2, v2, u3, v3)

    end
end


###################################################################
obseht = process_eht_obs(obs)
cphases =
joint = model(u=data.u/rad2μas, v=data.v/rad2μas,
               err=data.error,
               nx=nx, ny=ny)
jimage = image(u=data.u/rad2μas, v=data.v/rad2μas,
               err=data.error,
               nx=nx, ny=ny)

# Set the data
observations = (visr=data.visr, visi=data.visi,
                α=1.0, scx=60.0, scy=60.0, x0=0.0, y0=0.0)
lj = LogJoint(observations, jimage)


function forward_lj(x)
    res = DiffResults.GradientResult(x)
    cfg = forwardchunk(lj, lj.transform.dimension,10)
    ForwardDiff.gradient!(res, lj, x, cfg)
    return (DiffResults.value(res), DiffResults.gradient(res))
end

function reverse_lj(x)
   tape = reversecache(lj, lj.transform.dimension)
   res = DiffResults.GradientResult(x)
  ReverseDiff.gradient!(res, tape, x)
  return (DiffResults.value(res), DiffResults.gradient(res))
end

function fopt(x, grad)
    if length(grad) > 0
        f,g = forward_lj(x)
        grad .= g
    else
        f = lj(x)
    end

    return f
end

# Optimize!
#range = [(-5.0,5.0) for i in 1:nx*ny]
#pushfirst!(srange,(-5.0, 5.0),(-3.0, 3.0))#, (-50.0,50.0),(-50.0,50.0), (1.0, 5.0), (1.0, 5.0))

srange = [(-5.0,5.0) for i in 1:lj.transform.dimension]
nlopt = NLopt.Opt(:LD_LBFGS, length(srange))
lower_bounds!(nlopt, first.(srange))
upper_bounds!(nlopt, last.(srange))
max_objective!(nlopt, fopt)
xtol_rel!(nlopt,1e-12)
maxeval!(nlopt, 50_000)
(minf,minx,ret) = NLopt.optimize(nlopt, zeros(lj.transform.dimension))


results, divs = threaded_nlopt(256, fopt, srange, 100000)


#dim = lj.transform.dimension
#srange = [(-5.0,5.0) for i in 1:dim]
#results2, divs2 = threaded_bbopt(4, lj, srange, 100_000)
#res = optimize(x->-lj(x), results2[1], LBFGS())
#results, divs = threaded_opt(8, lj, srange,1)

# Construct the image
start = results[1]#divs[1] < divs2[1] ? results[1] : results2[1]
opt = lj.transform(start)
rimage = shifted(renormed(stretched(RImage(reshape(opt[:coeffs], nx,ny), SqExpKernel(opt[:ϵ])),
                                    opt[:scx], opt[:scy]), opt[:f]), opt[:x0], opt[:y0])
optcres = opt[:ring]
router = optcres[:R]
rinner = router*(1-optcres[:ψ])
img = ConcordanceCrescent(router, rinner, 0.0, optcres[:s])
mc = shifted(renormed(rotated(img, optcres[:ξ]), optcres[:f]), optcres[:x0], optcres[:y0])

mopt = rimage + mc
simopt = ROSE.stokesimage(mopt, 512, 512, 160.0, 160.0)

#Find a chi-square because yes...
echi2,nres = chi2(mopt, data)
vmodel = visibility.(Ref(mopt), data.u/rad2μas, data.v/rad2μas)
rchi2 = echi2/(2*length(data)-lj.transform.dimension)

# Make EHTImage from VIDA and save the fits
img = EHTImage(512,512, -160/512.0, 160/512.0, "SGRA", obs.ra, obs.dec, 1.33/1000, float(obs.mjd), Matrix(simopt.im)[:,end:-1:1])
save_fits(img, "map_fit_nx:$(nx)_ny:$(ny)_adaptivegrid.fits")
plot(img)

#Plot some fit statistics stuff for MAP
p1=plot_vis_comp(obs.data, vmodel)
savefig(p1,joinpath(@__DIR__, "dataplot_riaf_nx:$(nx)_ny:$(ny)_adaptivegrid.png"))

p2=plot_res_uv(data, nres, rchi2, 1)
savefig(p2, joinpath(@__DIR__, "chi2_uv_nx:$(nx)_ny:$(ny)_adaptivegrid.png"))

p3=plot_res_time(data, nres, rchi2, 1)
savefig(p3, joinpath(@__DIR__, "chi2_time_dirichlet_nx:$(nx)_ny:$(ny)_adaptivegrid.png"))

p4=plot_res_density(nres)
savefig(p4, joinpath(@__DIR__, "density_residuals_dirichlet_nx:$(nx)_ny:$(ny)_adaptivegrid.png"))



samples, stats = nuts_sample(start, lj, forward_lj, 2250, 1750)

schain = lj.transform.(samples)
histogram(getindex.(getindex.(schain, :ring),:f))
sims = make_sims(schain, rand(1:length(schain), 50), nx, ny)


p5=plot_mean(sims, 1)
savefig(p5, joinpath(@__DIR__,"mean_nx:$(nx)_ny:$(ny)_adaptivegrid.png"))


anim = plot_samples(sims, rchi2, 1)
gif(anim, joinpath(@__DIR__, "chain_nx:$(nx)_ny:$(ny)_adaptivegrid.gif"), fps=5)

df = DataFrame(schain)
df |> CSV.write(joinpath(@__DIR__, "chain__nx:$(nx)_ny:$(ny)_adaptivegrid.csv"))
