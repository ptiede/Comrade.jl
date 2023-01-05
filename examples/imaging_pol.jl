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
# obs = load_ehtim_uvfits(joinpath(@__DIR__, "PolarizedData/convention_test1/convention_test1.uvfits"), joinpath(@__DIR__, "PolarizedData/convention_test1/template_array.txt"))
# obs = load_ehtim_uvfits(joinpath(@__DIR__, "polarized_synthetic_data.uvfits"))
# obs = load_ehtim_uvfits(joinpath(@__DIR__, "PolarizedData/hops_lo_3601_M87+ALMArot.uvfits"), joinpath(@__DIR__, "PolarizedData/array.txt"))
obs = load_ehtim_uvfits(joinpath(@__DIR__, "PolarizedExamples/polarized_gaussian_all_corruptions.uvfits"),
                                joinpath(@__DIR__, "PolarizedExamples/array.txt"))

obs.add_scans()
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obsavg = scan_average(obs)
# extract log closure amplitudes and closure phases
dvis = extract_coherency(obsavg)



function model(θ, metadata)
    (;c, f, p, angparams, dRx, dRy, dLx, dLy, lgR, lgL, gpR, gpL) = θ
    (; grid, cache, pulse, tcache, gcache, dcache) = metadata
    # Construct the image model
    # produce Stokes images from parameters
    imgI = f.*c
    pimg = PoincareSphere2Map(imgI, p, angparams, grid)
    cimg = ContinuousImage(pimg, pulse)

    m = modelimage(cimg, cache)
    jT = jonesT(tcache)
    # calibration parameters
    gR = exp.(lgR).*cis.(gpR)
    gL = exp.(lgL).*cis.(gpL)
    G = jonesG(gR, gL, gcache)
    D = jonesD(complex.(dRx, dRy), complex.(dLx, dLy), dcache)
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


distphase = (AA = DiagonalVonMises([0.0], [inv(1e-4)]),
             AP = DiagonalVonMises([0.0], [inv(π^2)]),
             LM = DiagonalVonMises([0.0], [inv(π^2)]),
             AZ = DiagonalVonMises([0.0], [inv(π^2)]),
             JC = DiagonalVonMises([0.0], [inv(π^2)]),
             PV = DiagonalVonMises([0.0], [inv(π^2)]),
             SM = DiagonalVonMises([0.0], [inv(π^2)]),
           )

distD = (   AA = Normal(0.0, 0.1),
             AP = Normal(0.0, 0.1),
             LM = Normal(0.0, 0.1),
             AZ = Normal(0.0, 0.1),
             JC = Normal(0.0, 0.1),
             PV = Normal(0.0, 0.1),
             SM = Normal(0.0, 0.1),
           )


# Set up the cache structure
fovx = μas2rad(50.0)
fovy = μas2rad(50.0)
nx = 24
ny = floor(Int, fovy/fovx*nx)

grid = imagepixels(fovx, fovy, nx, ny)
buffer = IntensityMap(zeros(nx, ny), grid)
pulse = BSplinePulse{3}()
cache = create_cache(NFFTAlg(dvis), buffer, pulse)
tcache = TransformCache(dvis; add_fr=true, ehtim_fr_convention=false)
gcache = JonesCache(dvis, ScanSeg())
dcache = JonesCache(dvis, TrackSeg())
metadata = (;cache, grid, pulse, tcache, gcache, dcache)


prior = (
          c = ImageDirichlet(1.0, nx, ny),
          f = Uniform(0.7, 1.3),
          p = ImageUniform(nx, ny),
          angparams = ImageSphericalUniform(nx, ny),
          lgR = CalPrior(distamp, gcache),
          lgL = CalPrior(distamp, gcache),
          gpR = CalPrior(distphase, gcache),
          gpL = CalPrior(distphase, gcache),
          dRx = CalPrior(distD, dcache),
          dRy = CalPrior(distD, dcache),
          dLx = CalPrior(distD, dcache),
          dLy = CalPrior(distD, dcache),
          )



lklhd1 = RadioLikelihood(model1, metadata, dvis)

post1 = Posterior(lklhd1, prior)

tpost1 = asflat(post1)
ndim1 = dimension(tpost1)


ℓ1 = logdensityof(tpost1)

lklhd2 = RadioLikelihood(model2, metadata, dvis)

post2 = Posterior(lklhd2, prior)

tpost = asflat(post)
ndim = dimension(tpost)

ℓ = logdensityof(tpost)


# We will use HMC to sample the posterior.

using Zygote
f = OptimizationFunction(tpost, Optimization.AutoZygote())
prob = OptimizationProblem(f, rand(ndim) .- 0.5, nothing)
sol = solve(prob, LBFGS(); maxiters=10_000, callback=(x,p)->(@info ℓ(x); false), g_tol=1e-1)
xopt = transform(tpost, sol)


plot(model(xopt, metadata), fovx=fovx, fovy=fovy)
# Let's see how the fit looks

img = intensitymap(model(xopt, metadata), fovx, fovy, 128, 128)
plot(img)
Comrade.save(joinpath(@__DIR__, "DomTest/comrade_fit_3601_lo_allcorr.fits"), img)

residual(model(xopt, metadata), dvis)
plot(model(xopt, metadata), dvis)

# Let's also plot the calibration tables for gains
gL = Comrade.caltable(gcache, exp.(xopt.lgL))
plot(gL, layout=(3,3), size=(600,500))

gR = Comrade.caltable(gcache, exp.(xopt.lgR))
plot!(gR, layout=(3,3), size=(600,500))

# Let's also plot the calibration tables for gains
gpL = Comrade.caltable(gcache, (xopt.gpL))
plot(gpL, layout=(3,3), size=(600,500))

gpR = Comrade.caltable(gcache, (xopt.gpR))
plot!(gpR, layout=(3,3), size=(600,500))

# And the calibration tables for d-terms
dR = caltable(dcache, complex.(xopt.dRx, xopt.dRy))
dL = caltable(dcache, complex.(xopt.dLx, xopt.dLy))


using Measurements

using Pathfinder
res = pathfinder(
        ℓ, ℓ';
        init=sol.u .+ 0.05*randn(ndim),
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
hchain, stats = sample(tpost, AHMC(;metric, autodiff=AD.ZygoteBackend()), 10_000; nadapts=9_000, init_params=res.draws[:,1])

# Now plot the gain table with error bars
gamps = (hcat(hchain.dLx...))
mga = mean(gamps, dims=2)
sga = std(gamps, dims=2)

using Measurements
gmeas = measurement.(mga, sga)
ctable = caltable(dcache, vec(gmeas))
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


gAA = CSV.read(joinpath(@__DIR__, "DomTest/gain"))
