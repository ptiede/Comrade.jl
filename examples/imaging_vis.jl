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
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "simdata_3644_hops-b3_point+disk.uvfits"))
obs.add_scans()
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obsavg = scan_average(obs).add_fractional_noise(0.01)
obs_split = obsavg.split_obs()
# extract log closure amplitudes and closure phases
dvis = extract_vis(obsavg)

# Build the Model. Here we we a struct to hold some caches
# This will be useful to hold precomputed caches

function model(θ, metadata)
    (;c, f, fg, lgamp, gphase) = θ
    (; fovx, fovy, cache, gcache) = metadata
    # Construct the image model
    img = IntensityMap(f*c, fovx, fovy)
    cimg = ContinuousImage(img, BSplinePulse{3}())
    m = modelimage(cimg, cache)
    gaussian = fg*stretched(Gaussian(), μas2rad(1000.0), μas2rad(1000.0))
    # Now corrupt the model with Gains
    ga = exp.(lgamp)
    gp = cis.(gphase)
    return Comrade.GainModel(gcache, ga.*gp, m+gaussian)
end



# First we define the station gain priors
# distamp = (AA = Normal(0.0, 0.1),
#            AP = Normal(0.0, 0.1),
#            LM = Normal(0.0, 0.3),
#            AZ = Normal(0.0, 0.1),
#            JC = Normal(0.0, 0.1),
#            PV = Normal(0.0, 0.1),
#            SM = Normal(0.0, 0.1),
#            GL = Normal(0.0, 0.5)
#            )

# distphase = (AA = DiagonalVonMises([0.0], [inv(1e-4)]),
#              AP = DiagonalVonMises([0.0], [inv(π^2)]),
#              LM = DiagonalVonMises([0.0], [inv(π^2)]),
#              AZ = DiagonalVonMises([0.0], [inv(π^2)]),
#              JC = DiagonalVonMises([0.0], [inv(π^2)]),
#              PV = DiagonalVonMises([0.0], [inv(π^2)]),
#              SM = DiagonalVonMises([0.0], [inv(π^2)]),
#              GL = DiagonalVonMises([0.0], [inv(π^2)])
#            )

distamp = (AA = Normal(0.0, 0.1),
           PV = Normal(0.0, 0.1),
           AX = Normal(0.0, 0.1),
           MG = Normal(0.0, 0.1),
           LM = Normal(0.0, 0.3),
           MM = Normal(0.0, 0.1),
           SW = Normal(0.0, 0.1),
           GL = Normal(0.0, 0.5)
           )

distphase = (AA = DiagonalVonMises([0.0], [inv(1e-4)]),
             PV = DiagonalVonMises([0.0], [inv(π^2)]),
             AX = DiagonalVonMises([0.0], [inv(π^2)]),
             MG = DiagonalVonMises([0.0], [inv(π^2)]),
             LM = DiagonalVonMises([0.0], [inv(π^2)]),
             MM = DiagonalVonMises([0.0], [inv(π^2)]),
             SW = DiagonalVonMises([0.0], [inv(π^2)]),
             GL = DiagonalVonMises([0.0], [inv(π^2)])
           )



10
fovx = μas2rad(100.0)
fovy = μas2rad(100.0)
nx = 10
ny = floor(Int, fovy/fovx*nx)
prior = (
          c = ImageDirichlet(2.0, nx, ny),
          f = truncated(Normal(0.6, 0.05), 0.0, 1.0),
          fg = Uniform(0.0, 1.0),
          lgamp = Comrade.GainPrior(distamp, scantable(dvis)),
          gphase = Comrade.GainPrior(distphase, scantable(dvis))
        )


buffer = IntensityMap(zeros(nx, ny), fovx, fovy)
cache = create_cache(DFTAlg(dvis), buffer, BSplinePulse{3}())
gcache = GainCache(scantable(dvis))
metadata = (;cache, fovx, fovy, gcache)

lklhd = RadioLikelihood(model, metadata, dvis)

post = Posterior(lklhd, prior)

tpost = asflat(post)

ℓ = logdensityof(tpost)

# We will use HMC to sample the posterior.

ndim = dimension(tpost)
using Zygote
f = OptimizationFunction(tpost, Optimization.AutoZygote())
prob = OptimizationProblem(f, rand(ndim) .- 0.5, nothing)
ℓ = logdensityof(tpost)
sol = solve(prob, LBFGS(); maxiters=3000, callback=(x,p)->(@info ℓ(x); false), g_tol=1e-1)
xopt = transform(tpost, sol)


# Let's see how the fit looks
plot(model(xopt, metadata), fovx=fovx, fovy=fovy, title="MAP")

residual(model(xopt, metadata), dvis)
#residual(mms(xopt), dcphase)

# Let's also plot the gain curves
gt = Comrade.caltable(gcache, xopt.gphase)
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
