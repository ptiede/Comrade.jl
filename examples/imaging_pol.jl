using Pkg; Pkg.activate(@__DIR__)
#Pkg.add(url="https://github.com/ptiede/RadioImagePriors.jl")
using Comrade
using Distributions
using ComradeOptimization
using ComradeAHMC
using OptimizationBBO
using Plots
using StatsBase
using OptimizationOptimJL
using RadioImagePriors
using DistributionsAD

# load eht-imaging we use this to load eht data
load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
# obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "0316+413.2013.08.26.uvfits"))
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "polarized_synthetic_data.uvfits"))
obs.add_scans()
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obsavg = scan_average(obs)
obs_split = obsavg.split_obs()
# extract log closure amplitudes and closure phases
dvis = extract_coherency(obsavg)

# Build the Model. Here we we a struct to hold some caches
# This will be useful to hold precomputed caches

# function polimg(stokesI, pol, sphere)
#     I = stokesI
#     Q = stokesI .* pol .* getindex.(sphere, 1)
#     U = stokesI .* pol .* getindex.(sphere, 2)
#     V = stokesI .* pol .* getindex.(sphere, 3)
#     return (;I, Q, U, V)
# end

# function ChainRulesCore.rrule(::typeof(polimg), stokesI, pol, sphere)
#     s1 = getindex.(sphere, 1)
#     s2 = getindex.(sphere, 2)
#     s3 = getindex.(sphere, 3)
#     I = stokesI
#     Q = stokesI .* pol .* s1
#     U = stokesI .* pol .* s2
#     V = stokesI .* pol .* s3

#     function _polimg_pullback(Δ)
#         ΔI = Δ.I
#         ΔQ = Δ.Q
#         ΔU = Δ.U
#         ΔV = Δ.V

#         ΔstokesI = ΔI .+ ΔQ.*pol.*s1 .+ ΔU.*pol.*s2 .+ ΔV.*pol.*s3
#         Δpol = stokesI.*(ΔQ.*s1 .+ ΔU.*s2 .+ ΔV.*s3)
#         Δsphere =
#     end
# end

using StructArrays
function model(θ, metadata)
    (;c, f, p, angparams, lgR, lgL, dRa, dRp, dLa, dLp) = θ
    (; fovx, fovy, cache, gcache, tcache) = metadata
    # Construct the image model
    # produce Stokes images from parameters
    csa = angparams
    imgI = f.*c
    pimgI = imgI.*p
    I = IntensityMap(imgI, fovx, fovy)
    Q = IntensityMap(pimgI .* csa[1], fovx, fovy)
    U = IntensityMap(pimgI .* csa[2], fovx, fovy)
    V = IntensityMap(pimgI .* csa[3], fovx, fovy)

    pimg = PolIntensityMap(I, Q, U, V)
    cimg = ContinuousImage(pimg, BSplinePulse{3}())

    m = modelimage(cimg, cache)
    jT = jonesT(tcache)
    # calibration parameters
    gR = exp.(lgR .+ 0.0im)
    gL = exp.(lgL .+ 0.0im)
    G = jonesG(gR, gL, gcache)
    D = jonesD(dRa.*cis.(dRp), dLa.*cis.(dLp), dcache)
    J = G*D*jT
    return JonesModel(J, m, CirBasis())
end



# First we define the station gain priors
distamp = (AA = Normal(0.0, 0.1),
           AP = Normal(0.0, 0.1),
           LM = Normal(0.0, 0.3),
           AZ = Normal(0.0, 0.1),
           JC = Normal(0.0, 0.1),
           PV = Normal(0.0, 0.1),
           SM = Normal(0.0, 0.1),
           )

distamp2 = (AA = Normal(0.0, 1.0),
           AP = Normal(0.0,  1.0),
           LM = Normal(0.0,  1.0),
           AZ = Normal(0.0,  1.0),
           JC = Normal(0.0,  1.0),
           PV = Normal(0.0,  1.0),
           SM = Normal(0.0,  1.0),
           )


distphase = (AA = DiagonalVonMises([0.0], [inv(π^2)]),
             AP = DiagonalVonMises([0.0], [inv(π^2)]),
             LM = DiagonalVonMises([0.0], [inv(π^2)]),
             AZ = DiagonalVonMises([0.0], [inv(π^2)]),
             JC = DiagonalVonMises([0.0], [inv(π^2)]),
             PV = DiagonalVonMises([0.0], [inv(π^2)]),
             SM = DiagonalVonMises([0.0], [inv(π^2)]),
           )


# Set up the cache structure
fovx = μas2rad(65.0)
fovy = μas2rad(65.0)
nx = 10
ny = floor(Int, fovy/fovx*nx)

buffer = IntensityMap(zeros(nx, ny), fovx, fovy)
cache = create_cache(DFTAlg(dvis), buffer, BSplinePulse{3}())
tcache = TransformCache(dvis)
gcache = JonesCache(dvis, ScanSeg())
dcache = JonesCache(dvis, TrackSeg())
metadata = (;cache, fovx, fovy, tcache, gcache, dcache)



prior = (
          c = ImageDirichlet(2.0, nx, ny),
          f = Uniform(0.1, 2.0),
          p = ImageUniform(nx, ny),
          angparams = ImageSphericalUniform(nx, ny),
          lgR = CalPrior(distamp, gcache),
          lgL = CalPrior(distamp, gcache),
          dRa = CalPrior(distamp2, dcache),
          dRp = CalPrior(distphase, dcache),
          dLa = CalPrior(distamp2, dcache),
          dLp = CalPrior(distphase, dcache),
          )



lklhd = RadioLikelihood(model, metadata, dvis)

post = Posterior(lklhd, prior)

tpost = asflat(post)
ndim = dimension(tpost)

ℓ = logdensityof(tpost)

# We will use HMC to sample the posterior.

using Zygote
f = OptimizationFunction(tpost, Optimization.AutoZygote())
prob = OptimizationProblem(f, randn(ndim)*0.2, nothing)
sol = solve(prob, LBFGS(); maxiters=3000, callback=(x,p)->(@info ℓ(x); false), g_tol=1e-1)
xopt = transform(tpost, sol)





# Let's see how the fit looks
plot(model(xopt, metadata), fovx=fovx, fovy=fovy)
timg, hdr = Comrade.load(joinpath(@__DIR__, "polarized_synthetic_data.fits"), StokesIntensityMap)
plot(timg, xlims=(-32.5, 32.5), ylims=(-32.5, 32.5))

Comrade.residuals(model(xopt, metadata), dvis)
#residual(mms(xopt), dcphase)

# Let's also plot the gain curves
gt = Comrade.caltable(gcache, xopt.lgL)
plot(gt, layout=(3,3), size=(600,500))

gt = Comrade.caltable(gcache, exp.(xopt.lgamp))
plot(gt, layout=(3,3), size=(600,500))


using Measurements

using Pathfinder
res = pathfinder(
        ℓ, ℓ';
        init=sol.u .+ 0.01*randn(ndim),
        dim = ndim,
        optimizer=LBFGS(m=6),
        g_tol=1e-1,
        maxiters=1000,
        callback = (x,p)->(l = ℓ(x); @info l; isnan(l))
)

vis = dvis[:measurement]
err = dvis[:error]
mvis = visibilities(model(xopt, metadata), arrayconfig(dvis))
uvdist = hypot.(values(getuv(dvis))...)

p1 = scatter(uvdist, real.(getindex.(vis, 1, 1)), yerr=getindex.(err, 1, 1), title="RR", color=:blue, marker=:circ, label="Data real")
scatter!(uvdist, imag.(getindex.(vis, 1, 1)), yerr=getindex.(err, 1, 1), color=:orange, marker=:circ, label="Data imag")
scatter!(uvdist, real.(getindex.(vis, 1, 1)), color=:cyan, alpha=0.5, marker=:square, label="Model real")
scatter!(uvdist, imag.(getindex.(vis, 1, 1)), color=:yellow, alpha=0.5, marker=:square, label="Model imag")

p2 = scatter(uvdist, real.(getindex.(vis, 2, 1)), yerr=getindex.(err, 2, 1), title="LR", color=:blue, marker=:circ, label="Data real RR")
scatter!(uvdist, imag.(getindex.(vis, 2, 1)), yerr=getindex.(err, 2, 1), color=:orange, marker=:circ, label="Data imag RR")
scatter!(uvdist, real.(getindex.(vis, 2, 1)), color=:cyan, alpha=0.5, marker=:square, label="Model real RR")
scatter!(uvdist, imag.(getindex.(vis, 2, 1)), color=:yellow, alpha=0.5, marker=:square, label="Model imag RR", legend=:false)

p3 = scatter(uvdist, real.(getindex.(vis, 1, 2)), yerr=getindex.(err, 1, 2), title="RL", color=:blue, marker=:circ, label="Data real RR")
scatter!(uvdist, imag.(getindex.(vis, 1, 2)), yerr=getindex.(err, 1, 2), color=:orange, marker=:circ, label="Data imag RR")
scatter!(uvdist, real.(getindex.(vis, 1, 2)), color=:cyan, alpha=0.5, marker=:square, label="Model real RR")
scatter!(uvdist, imag.(getindex.(vis, 1, 2)), color=:yellow, alpha=0.5, marker=:square, label="Model imag RR", legend=false)

p4 = scatter(uvdist, real.(getindex.(vis, 2, 2)), yerr=getindex.(err, 2, 2), title="LL", color=:blue, marker=:circ, label="Data real RR")
scatter!(uvdist, imag.(getindex.(vis, 2, 2)), yerr=getindex.(err, 2, 2), color=:orange, marker=:circ, label="Data imag RR")
scatter!(uvdist, real.(getindex.(vis, 2, 2)), color=:cyan, alpha=0.5, marker=:square, label="Model real RR")
scatter!(uvdist, imag.(getindex.(vis, 2, 2)), color=:yellow, alpha=0.5, marker=:square, label="Model imag RR", legend=false)

plot(p1, p2, p3, p4, layout=(2,2))

# now we sample using hmc
using LinearAlgebra
metric = DenseEuclideanMetric(res.fit_distribution.Σ)
hchain, stats = sample(post, AHMC(;metric, autodiff=AD.ZygoteBackend()), 2_000; nadapts=1_000, init_params=xopt)

# Now plot the gain table with error bars
gamps = (hcat(hchain.gphase...))
mga = mean(gamps, dims=2)
sga = std(gamps, dims=2)

using Measurements
gmeas = measurement.(mga, sga)
ctable = caltable(gcache, vec(gmeas))
plot(ctable, layout=(3,3), size=(600,500))

# This takes about 1.75 hours on my laptop. Which isn't bad for a 575 dimensional model!

# Plot the mean image and standard deviation image
using StatsBase
samples = model.(sample(hchain, 50), Ref(metadata))
imgs = intensitymap.(samples, fovx, fovy, 128,  128)

mimg, simg = mean_and_std(imgs)

p1 = plot(mimg, title="Mean")
p2 = plot(simg,  title="Std. Dev.")
p3 = plot(mimg./simg,  title="SNR")
p4 = plot(simg./mimg,  title="Fractional Error")

plot(p1,p2,p3,p4, layout=(2,2), size=(800,800))
savefig("3c84_complex_vis_2min_avg.png")
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
