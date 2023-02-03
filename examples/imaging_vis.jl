# # Stokes I simultaneous Image and Instrument Modeling

# In this tutorial we will create a preliminary reconstruction of the 2017 M87 data on April 6
# by simultaneously creating and image and model for the instrument. By instrument model we
# mean something akin to self-calibration in traditional VLBI imaging terminology. However,
# unlike traditional self-cal we will at each point in our parameter space effectively explore
# the possible self-cal solutions. This will allow us to constrain and marginalize over the
# instrument effects such as time variable gains.

# ## Introduction to Complex Visibility Fitting


using Pkg; Pkg.activate(@__DIR__)

# To get started we will load Comrade
using Comrade

# ## Load the Data
# To download the data visit https://doi.org/10.25739/g85n-f134
# To load the eht-imaging obsdata object we do:
obs = load_ehtim_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))

# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      coherent.
#   - Add 1% systematic noise to deal with calibration issues that cause 1% non-closing errors.
obs = scan_average(obs).add_fractional_noise(0.015).flag_uvdist(uv_min=0.1e9)

# Now we extract our complex visibilities.
dvis = extract_vis(obs)

# ##Building the Model/Posterior
# Unlike Here we we a struct to hold some caches
# This will be useful to hold precomputed caches

function model(θ, metadata)
    (;c, lgamp, gphase) = θ
    (; grid, cache, gcache) = metadata
    # Construct the image model we fix the flux to 0.6 Jy in this case
    img = IntensityMap(0.6*c, grid)
    cimg = ContinuousImage(img, BSplinePulse{3}())
    m = modelimage(cimg, cache)
    # Now corrupt the model with Gains
    j = @fastmath jonesStokes(exp.(lgamp).*cis.(gphase), gcache)
    return JonesModel(j, m)
end



# First we define the station gain priors
distamp = (AA = Normal(0.0, 0.1),
           AP = Normal(0.0, 0.1),
           LM = Normal(0.0, 1.0),
           AZ = Normal(0.0, 0.1),
           JC = Normal(0.0, 0.1),
           PV = Normal(0.0, 0.1),
           SM = Normal(0.0, 0.1),
           )

distphase = (AA = DiagonalVonMises(0.0, inv(1e-3)),
             AP = DiagonalVonMises(0.0, inv(π^2)),
             LM = DiagonalVonMises(0.0, inv(π^2)),
             AZ = DiagonalVonMises(0.0, inv(π^2)),
             JC = DiagonalVonMises(0.0, inv(π^2)),
             PV = DiagonalVonMises(0.0, inv(π^2)),
             SM = DiagonalVonMises(0.0, inv(π^2)),
           )



fovx = μas2rad(75.0)
fovy = μas2rad(75.0)
nx = 7
ny = floor(Int, fovy/fovx*nx)

buffer = IntensityMap(zeros(nx, ny), fovx, fovy)
cache = create_cache(NFFTAlg(dvis), buffer, BSplinePulse{3}())
gcache = JonesCache(dvis, ScanSeg())
metadata = (;cache, fovx, fovy, gcache)


X, Y = imagepixels(buffer)
prior = (
          c = ImageDirichlet(2.0, nx, ny),
          lgamp = CalPrior(distamp, gcache),
          gphase = CalPrior(distphase, gcache)
        )



lklhd = RadioLikelihood(model, metadata, dvis)

post = Posterior(lklhd, prior)

tpost = asflat(post)

ndim = dimension(tpost)
ℓ = logdensityof(tpost)

# We will use HMC to sample the posterior.

using Zygote
f = OptimizationFunction(tpost, Optimization.AutoZygote())
prob = OptimizationProblem(f, rand(ndim) .- 0.5, nothing)
sol = solve(prob, LBFGS(); maxiters=3000, callback=(x,p)->(@info ℓ(x); false), g_tol=1e-1)
xopt = transform(tpost, sol)


# Let's see how the fit looks
img = intensitymap(model(xopt, metadata), fovx, fovy, 128, 128)
plot(img)

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
        init=sol.u .+ 0.035*randn(ndim),
        dim = ndim,
        optimizer=LBFGS(m=6),
        g_tol=1e-1,
        maxiters=1500,
        callback = (x,p)->(l = ℓ(x); @info l; isnan(l))
)



# now we sample using hmc
using LinearAlgebra
metric = DenseEuclideanMetric((res.fit_distribution.Σ))
hchain, stats = sample(post, AHMC(;metric, autodiff=AD.ZygoteBackend()), 12_000; nadapts=11_000, init_params=transform(tpost, res.draws[:,1]))

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
imgs = intensitymap.(samples, μas2rad(115.0), μas2rad(115.0), 128,  128)

mimg, simg = mean_and_std(imgs)

using CairoMakie
function Makie.convert_arguments(::SurfaceLike, img::Comrade.IntensityMap)
    return rad2μas.(values(imagepixels(img)))..., Comrade.baseimage(img)
end

fig = Figure(;resolution=(400,400))
ax = Axis(fig[1,1], xreversed=true, aspect=DataAspect())
hidedecorations!(ax)
image!(ax, mimg, colormap=:afmhot)
lines!(ax, [15.0, 55.0], [-45.0, -45.0], color=:white, linewidth=3)
text!(ax, 35.0, -50.0, text=L"$40\,\mu$as", color=:white, align=(:center, :center), fontsize=20)
fig
save("test.png", fig)

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
