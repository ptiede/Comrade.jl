import Pkg; #hide
__DIR = @__DIR__ #hide
pkg_io = open(joinpath(__DIR, "pkg.log"), "w") #hide
Pkg.activate(__DIR; io = pkg_io) #hide
Pkg.develop(; path = joinpath(__DIR, "..", "..", ".."), io = pkg_io) #hide
Pkg.instantiate(; io = pkg_io) #hide
Pkg.precompile(; io = pkg_io) #hide
close(pkg_io) #hide


# # Geometric Modeling of EHT Data

# `Comrade` has been designed to work with the EHT and ngEHT.
# In this tutorial, we will show how to reproduce some of the results
# from [EHTC VI 2019](https://iopscience.iop.org/article/10.3847/2041-8213/ab1141).

# In EHTC VI, they considered fitting simple geometric models to the data
# to estimate the black hole's image size, shape, brightness profile, etc.
# In this tutorial, we will construct a similar model and fit it to the data in under
# 50 lines of code (sans comments). To start, we load Comrade and some other packages we need.

# To get started we load Comrade.
#-

using Comrade


# Currently we use eht-imaging for data management, however this will soon be replaced
# by a pure Julia solution.
using Pyehtim
# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(42)
#-
# The next step is to load the data. We will use the publically
# available M 87 data which can be downloaded
# from [cyverse](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/EHTC_FirstM87Results_Apr2019).
# For an introduction to data loading, see [Loading Data into Comrade](@ref).
obs = ehtim.obsdata.load_uvfits(joinpath(__DIR, "..", "..", "Data", "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))

# Now we will kill 0-baselines since we don't care about large-scale flux and
# since we know that the gains in this dataset are coherent across a scan, we make scan-average data
obs = Pyehtim.scan_average(obs.flag_uvdist(uv_min = 0.1e9)).add_fractional_noise(0.02)

# Now we extract the data products we want to fit
dlcamp, dcphase = extract_table(obs, LogClosureAmplitudes(; snrcut = 3.0), ClosurePhases(; snrcut = 3.0))

# !!!warn
#    We remove the low-snr closures since they are very non-gaussian. This can create rather
#    large biases in the model fitting since the likelihood has much heavier tails that the
#    usual Gaussian approximation.
#-

# For the image model, we will use a modified `MRing`, a
# infinitely thin delta ring with an azimuthal structure given by a Fourier expansion.
# To give the MRing some width, we will convolve the ring with a Gaussian and add an
# additional gaussian to the image to model any non-ring flux.
# Comrade expects that any model function must accept a named tuple and returns  must always
# return an object that implements the [VLBISkyModels Interface](https://ehtjulia.github.io/VLBISkyModels.jl/stable/interface/)
#-
function sky(θ, p)
    (; radius, width, ma, mp, τ, ξτ, f, σG, τG, ξG, xG, yG) = θ
    α = ma .* cos.(mp .- ξτ)
    β = ma .* sin.(mp .- ξτ)
    ring = f * smoothed(modify(MRing(α, β), Stretch(radius, radius * (1 + τ)), Rotate(ξτ)), width)
    g = (1 - f) * shifted(rotated(stretched(Gaussian(), σG, σG * (1 + τG)), ξG), xG, yG)
    return ring + g
end

# To construct our likelihood `p(V|M)` where `V` is our data and `M` is our model, we use the `RadioLikelihood` function.
# The first argument of `RadioLikelihood` is always a function that constructs our Comrade
# model from the set of parameters `θ`.

# We now need to specify the priors for our model. The easiest way to do this is to
# specify a NamedTuple of distributions:

using Distributions, VLBIImagePriors
prior = (
    radius = Uniform(μas2rad(10.0), μas2rad(30.0)),
    width = Uniform(μas2rad(1.0), μas2rad(10.0)),
    ma = (Uniform(0.0, 0.5), ),
    mp = (Uniform(0, 2π),),
    τ = Uniform(0.0, 1.0),
    ξτ = Uniform(0.0, π),
    f = Uniform(0.0, 1.0),
    σG = Uniform(μas2rad(1.0), μas2rad(100.0)),
    τG = Exponential(1.0),
    ξG = Uniform(0.0, 1π),
    xG = Uniform(-μas2rad(80.0), μas2rad(80.0)),
    yG = Uniform(-μas2rad(80.0), μas2rad(80.0)),
)

# Note that for `α` and `β` we use a product distribution to signify that we want to use a
# multivariate uniform for the mring components `α` and `β`. In general the structure of the
# variables is specified by the prior. Note that this structure must be compatible with the
# model definition `model(θ)`.

# We can now construct our Sky model, which typically takes a model, prior and the
# on sky grid. Note that since our model is analytic the grid is not directly used when
# computing visibilities.
skym = SkyModel(sky, prior, imagepixels(μas2rad(200.0), μas2rad(200.0), 128, 128))


# In this tutorial we will be using closure products as our data. As such we do not need to specify a
# instrument model, since for stokes I imaging, the likelihood is approximately invariant to the instrument
# model.
post = VLBIPosterior(skym, dlcamp, dcphase)

# !!! note
#     When fitting visibilities a instrument is required, and a reader can refer to
#     [Stokes I Simultaneous Image and Instrument Modeling](@ref).


# This constructs a posterior density that can be evaluated by calling `logdensityof`.
# For example,

logdensityof(
    post, (
        sky = (
            radius = μas2rad(20.0),
            width = μas2rad(10.0),
            ma = (0.3, ),
            mp = (π / 2, ),
            τ = 0.1,
            ξτ = π / 2,
            f = 0.6,
            σG = μas2rad(50.0),
            τG = 0.1,
            ξG = 0.5,
            xG = 0.0,
            yG = 0.0,
        ),
    )
)

# ## Reconstruction

# Now that we have fully specified our model, we now will try to find the optimal reconstruction
# of our model given our observed data.

# Currently, `post` is in **parameter** space. Often optimization and sampling algorithms
# want it in some modified space. For example, nested sampling algorithms want the
# parameters in the unit hypercube. To transform the posterior to the unit hypercube, we
# can use the `ascube` function

cpost = ascube(post)

# If we want to flatten the parameter space and move from constrained parameters to (-∞, ∞)
# support we can use the `asflat` function

fpost = asflat(post)

# These transformed posterior expect a vector of parameters. For example, we can draw from the
# prior in our usual parameter space
p = prior_sample(rng, post)

# and then transform it to transformed space using T
logdensityof(cpost, Comrade.inverse(cpost, p))
logdensityof(fpost, Comrade.inverse(fpost, p))

# note that the log densit is not the same since the transformation has causes a jacobian to ensure volume is preserved.

# ### Finding the Optimal Image

# Typically, most VLBI modeling codes only care about finding the optimal or best guess
# image of our posterior `post` To do this, we will use [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/) and
# specifically the [`BlackBoxOptim.jl`](https://github.com/robertfeldt/BlackBoxOptim.jl) package. For Comrade, this workflow is
# we use the [`comrade_opt`](@ref) function.

using Optimization
using OptimizationBBO
xopt, sol = comrade_opt(post, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters = 50_000);

# Given this we can now plot the optimal image or the *maximum a posteriori* (MAP) image.

using DisplayAs
using CairoMakie
g = imagepixels(μas2rad(200.0), μas2rad(200.0), 256, 256)
fig = imageviz(intensitymap(skymodel(post, xopt), g), colormap = :afmhot, size = (500, 400));
DisplayAs.Text(DisplayAs.PNG(fig))


# ### Quantifying the Uncertainty of the Reconstruction

# While finding the optimal image is often helpful, in science, the most important thing is to
# quantify the certainty of our inferences. This is the goal of Comrade. In the language
# of Bayesian statistics, we want to find a representation of the posterior of possible image
# reconstructions given our choice of model and the data.
#
# Comrade provides several sampling and other posterior approximation tools. To see the
# list, please see the Libraries section of the docs. For this example, we will be using
# [Pigeons.jl](https://github.com/Julia-Tempering/Pigeons.jl) which is a state-of-the-art
# parallel tempering sampler that enables global exploration of the posterior. For smaller dimension
# problems (< 100) we recommend using this sampler, especially if you have access to > 1 core.
using Pigeons
pt = pigeons(target = cpost, explorer = SliceSampler(), record = [traces, round_trip, log_sum_ratio], n_chains = 16, n_rounds = 8)


# That's it! To finish it up we can then plot some simple visual fit diagnostics.
# First we extract the MCMC chain for our posterior.
chain = sample_array(cpost, pt)

# First to plot the image we call
using DisplayAs #hide
imgs = intensitymap.(skymodel.(Ref(post), sample(chain, 100)), Ref(g))
fig = imageviz(imgs[end], colormap = :afmhot)
DisplayAs.Text(DisplayAs.PNG(fig))


# What about the mean image? Well let's grab 100 images from the chain, where we first remove the
# adaptation steps since they don't sample from the correct posterior distribution
meanimg = mean(imgs)
fig = imageviz(meanimg, colormap = :afmhot);
DisplayAs.Text(DisplayAs.PNG(fig))


# That looks similar to the EHTC VI, and it took us no time at all!. To see how well the
# model is fitting the data we can plot the model and data products. As of Comrade 0.11.7 Makie
# is the preferred plotting tool. For plotting data there are two classes of functions
#  - `baselineplot` which gives complete control of plotting
#  - `plotfields, plotfields!` which are more automated and limited but will automatically add
#     labels, legends, titles etc.
# A reasonable workflow is to use `plotfields` to set up the initial figure and axis labels and then
# then use `baselineplot!` to add additional plots to the axis. For example,
lcsim, cpsim = simulate_observation(post, xopt; add_thermal_noise = false)
fig, ax1 = plotfields(lcsim, uvdist, measwnoise, scatter_kwargs = (;marker=:circle, label="MAP"), figure_kwargs = (;size=(800,300)), legend=false);
baselineplot!(ax1, dlcamp, uvdist, measurement, marker = :+, color = :black, label = "Data")
ax2, = plotfields!(fig[1,2], cpsim, uvdist, mod2pi ∘ measwnoise, scatter_kwargs = (;marker=:circle, label="MAP"), axis_kwargs = (ylabel = "Closure Phase (rad)",))
baselineplot!(ax2, dcphase, uvdist, mod2pi ∘ measurement, marker = :+, color = :black, label = "Data")
axislegend(ax1, framevisible = false)
DisplayAs.Text(DisplayAs.PNG(fig))

# We can also plot random draws from the posterior predictive distribution.
# The posterior predictive distribution create a number of synthetic observations that
# are marginalized over the posterior.
fig = Figure(; size = (800, 300))
ax1 = Axis(fig[1, 1], xlabel = "√Quadrangle Area", ylabel = "Log Closure Amplitude")
ax2 = Axis(fig[1, 2], xlabel = "√Triangle Area", ylabel = "Closure Phase (rad)")
for i in 1:10
    mobs = simulate_observation(post, sample(chain, 1)[1])
    mlca = mobs[1]
    mcp = mobs[2]
    baselineplot!(ax1, mlca, uvdist, measurement, color = :grey, alpha = 0.2)
    baselineplot!(ax2, mcp, uvdist, mod2pi ∘ measurement, color = :grey, alpha = 0.2)
end
baselineplot!(ax1, dlcamp, uvdist, measurement, marker = :x)
baselineplot!(ax2, dcphase, uvdist, mod2pi ∘ measurement, marker = :x)
DisplayAs.Text(DisplayAs.PNG(fig))


# Finally, we can also put everything onto a common scale and plot the normalized residuals.
# The normalied residuals are the difference between the data
# and the model, divided by the data's error:
rd = residuals(post, chain[end])
fig, ax = plotfields(rd[1], uvdist, :res, axis_kwargs = (;ylabel = "Norm. Res. LCA"))
plotfields!(fig[2,1], rd[2], uvdist, :res, axis_kwargs = (;ylabel = "Norm. Res. CP"))
DisplayAs.Text(DisplayAs.PNG(fig))
