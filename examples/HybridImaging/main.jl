# # Hybrid Imaging of a Black Hole

# In this tutorial, we will use **hybrid imaging** to analyze the 2017 EHT data.
# By hybrid imaging, we mean decomposing the model into simple geometric models, e.g., rings
# and such, plus a rasterized image model to soak up the additional structure.
# This approach was first developed in [`BB20`](https://iopscience.iop.org/article/10.3847/1538-4357/ab9c1f)
# and applied to EHT 2017 data. We will use a similar model in this tutorial.

# ## Introduction to Hybrid modeling and imaging
# The benefit of using a hybrid-based modeling approach is the effective compression of
# information/parameters when fitting the data. Hybrid modeling requires the user to
# incorporate specific knowledge of how you expect the source to look like. For instance
# for M87, we expect the image to be dominated by a ring-like structure. Therefore, instead
# of using a high-dimensional raster to recover the ring, we can use a ring model plus
# a raster to soak up the additional degrees of freedom.
# This is the approach we will take in this tutorial to analyze the April 6 2017 EHT data
# of M87.

import Pkg #hide
__DIR = @__DIR__ #hide
pkg_io = open(joinpath(__DIR, "pkg.log"), "w") #hide
Pkg.activate(__DIR; io=pkg_io) #hide
Pkg.develop(; path=joinpath(__DIR, "..", ".."), io=pkg_io) #hide
Pkg.instantiate(; io=pkg_io) #hide
Pkg.precompile(; io=pkg_io) #hide
close(pkg_io) #hide

ENV["GKSwstype"] = "nul" #hide


# ## Loading the Data

# To get started we will load Comrade
using Comrade

# ## Load the Data


using Pyehtim

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(42)


# To download the data visit https://doi.org/10.25739/g85n-f134
# To load the eht-imaging obsdata object we do:
obs = ehtim.obsdata.load_uvfits(joinpath(__DIR, "../Data/SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))

# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      coherent.
obs = scan_average(obs).add_fractional_noise(0.02)

# For this tutorial we will once again fit complex visibilities since they
# provide the most information once the telescope/instrument model are taken
# into account.
dvis  = extract_table(obs, ComplexVisibilities())

# ## Building the Model/Posterior

# Now we build our intensity/visibility model. That is, the model that takes in a
# named tuple of parameters and perhaps some metadata required to construct the model.
# For our model, we will use a raster or `ContinuousImage` model, an `m-ring` model,
# and a large asymmetric Gaussian component to model the unresolved short-baseline flux.

function sky(θ, metadata)
    (;c, σimg, f, r, σ, τ, ξτ, ma, mp, fg) = θ
    (;ftot, grid) = metadata
    ## Form the image model
    ## First transform to simplex space first applying the non-centered transform
    rast = (ftot*f*(1-fg)).*to_simplex(CenteredLR(), σimg.*c)
    mimg = ContinuousImage(rast, grid, BSplinePulse{3}())
    ## Form the ring model
    α = ma.*cos.(mp .- ξτ)
    β = ma.*sin.(mp .- ξτ)
    ring = smoothed(modify(MRing(α, β), Stretch(r, r*(1+τ)), Rotate(ξτ), Renormalize((ftot*(1-f)*(1-fg)))), σ)
    gauss = modify(Gaussian(), Stretch(μas2rad(250.0)), Renormalize(ftot*f*fg))
    ## We group the geometric models together for improved efficiency. This will be
    ## automated in future versions.
    return mimg + (ring + gauss)
end

# Unlike other imaging examples
# (e.g., [Imaging a Black Hole using only Closure Quantities](@ref)) we also need to include
# a model for the instrument, i.e., gains as well. The gains will be broken into two components
#   - Gain amplitudes which are typically known to 10-20%, except for LMT, which has amplitudes closer to 50-100%.
#   - Gain phases which are more difficult to constrain and can shift rapidly.

using VLBIImagePriors
using Distributions, DistributionsAD
G = SingleStokesGain() do x
    lg = x.lg
    gp = x.gp
    return exp(lg + 1im*gp)
end

intpr = (
    lg= ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.2)); LM = IIDSitePrior(ScanSeg(), Normal(0.0, 1.0))),
    gp= ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))); refant=SEFDReference(0.0))
        )
intmodel = InstrumentModel(G, intpr)


# Before we move on, let's go into the `model` function a bit. This function takes two arguments
# `θ` and `metadata`. The `θ` argument is a named tuple of parameters that are fit
# to the data. The `metadata` argument is all the ancillary information we need to construct the model.
# For our hybrid model, we will need two variables for the metadata, a `grid` that specifies
# the locations of the image pixels and a `cache` that defines the algorithm used to calculate
# the visibilities given the image model. This is required since `ContinuousImage` is most easily
# computed using number Fourier transforms like the [`NFFT`](https://github.com/JuliaMath/NFFT.jl)
# or [FFT](https://github.com/JuliaMath/FFTW.jl).
# To combine the models, we use `Comrade`'s overloaded `+` operators, which will combine the
# images such that their intensities and visibilities are added pointwise.

# Now let's define our metadata. First we will define the cache for the image. This is
# required to compute the numerical Fourier transform.
fovxy  = μas2rad(200.0)
npix   = 32
g      = imagepixels(fovxy, fovxy, npix, npix)



# Part of hybrid imaging is to force a scale separation between
# the different model components to make them identifiable.
# To enforce this we will set the
# length scale of the raster component equal to the beam size of the telescope in units of
# pixel length, which is given by
beam = beamsize(dvis)
rat = (beam/(step(g.X)))
cprior = GaussMarkovRandomField(rat, size(g); order=2)

# additionlly we will fix the standard deviation of the field to unity and instead
# use a pseudo non-centered parameterization for the field.
# GaussMarkovRandomField(meanpr, 0.1*rat, 1.0, crcache)


# Finally we can put form the total model prior
skyprior = (
          c  = cprior,
          ## We use a strong smoothing prior since we want to limit the amount of high-frequency structure in the raster.
          σimg = truncated(Normal(0.0, 1.0); lower=0.01),
          f  = Uniform(0.0, 1.0),
          r  = Uniform(μas2rad(10.0), μas2rad(30.0)),
          σ  = Uniform(μas2rad(0.1), μas2rad(10.0)),
          τ  = truncated(Normal(0.0, 0.1); lower=0.0, upper=1.0),
          ξτ = Uniform(-π/2, π/2),
          ma = ntuple(_->Uniform(0.0, 0.5), 2),
          mp = ntuple(_->Uniform(0.0, 2π), 2),
          fg = Uniform(0.0, 1.0),
        )

# Now we form the metadata
skymetadata = (;ftot=1.1, grid = g)
skym = SkyModel(sky, skyprior, g; metadata=skymetadata)

# This is everything we need to specify our posterior distribution, which our is the main
# object of interest in image reconstructions when using Bayesian inference.
post = VLBIPosterior(skym, intmodel, dvis)

# To sample from our prior we can do
xrand = prior_sample(rng, post)

# and then plot the results
using DisplayAs #hide
import CairoMakie as CM
CM.activate!(type = "png", px_per_unit=1) #hide
gpl = imagepixels(μas2rad(200.0), μas2rad(200.0), 128, 128)
fig = imageviz(intensitymap(skymodel(post, xrand), gpl), size=(400, 400));
DisplayAs.Text(DisplayAs.PNG(fig)) #hide


# ## Reconstructing the Image

# To sample from this posterior, it is convenient to first move from our constrained parameter space
# to an unconstrained one (i.e., the support of the transformed posterior is (-∞, ∞)). This is
# done using the `asflat` function.
tpost = asflat(post)

# We can now also find the dimension of our posterior or the number of parameters we will sample.
# !!! warning
#     This can often be different from what you would expect. This is especially true when using
#     angular variables, where we often artificially increase the dimension
#     of the parameter space to make sampling easier.
#-
ndim = dimension(tpost)

# Now we optimize using LBFGS
using ComradeOptimization
using OptimizationOptimJL
using Zygote
f = OptimizationFunction(tpost, Optimization.AutoZygote())
prob = Optimization.OptimizationProblem(f, prior_sample(rng, tpost), nothing)
sol = solve(prob, LBFGS(); maxiters=2000, g_tol=1e-0);


# Before we analyze our solution we first need to transform back to parameter space.
xopt = transform(tpost, sol)

# First we will evaluate our fit by plotting the residuals
using Plots
fig = residual(post, xopt);
DisplayAs.Text(DisplayAs.PNG(fig)) #hide

# and now closure phases
#-

# Now these residuals look a bit high. However, it turns out this is because the MAP is typically
# not a great estimator and will not provide very predictive measurements of the data. We
# will show this below after sampling from the posterior.
imageviz(intensitymap(skymodel(post, xopt), gpl), figure=(;resolution=(500, 400),))


# We will now move directly to sampling at this point.
using ComradeAHMC
metric = DiagEuclideanMetric(ndim)
chain = sample(rng, post, AHMC(;metric, autodiff=Val(:Enzyme)), 700; n_adapts=500, progress=true);

# We then remove the adaptation/warmup phase from our chain
chain = chain[501:end]

# !!! warning
#     This should be run for 4-5x more steps to properly estimate expectations of the posterior
#-

# Now lets plot the mean image and standard deviation images.
# To do this we first clip the first 250 MCMC steps since that is during tuning and
# so the posterior is not sampling from the correct sitesary distribution.

using StatsBase
msamples = skymodel.(Ref(post), chain[begin:2:end]);

# The mean image is then given by
imgs = intensitymap.(msamples, Ref(gpl))
fig = imageviz(mean(imgs), colormap=:afmhot, size=(400, 300));
DisplayAs.Text(DisplayAs.PNG(fig)) #hide
#-
fig = imageviz(std(imgs), colormap=:batlow, size=(400, 300));
DisplayAs.Text(DisplayAs.PNG(fig)) #hide

#-
#
# We can also split up the model into its components and analyze each separately
comp = Comrade.components.(msamples)
ring_samples = getindex.(comp, 2)
rast_samples = first.(comp)
ring_imgs = intensitymap.(ring_samples, Ref(gpl));
rast_imgs = intensitymap.(rast_samples, Ref(gpl));

ring_mean, ring_std = mean_and_std(ring_imgs);
rast_mean, rast_std = mean_and_std(rast_imgs);

fig = CM.Figure(; resolution=(400, 400));
axes = [CM.Axis(fig[i, j], xreversed=true, aspect=CM.DataAspect()) for i in 1:2, j in 1:2]
CM.image!(axes[1,1], ring_mean, colormap=:afmhot); axes[1,1].title = "Ring Mean"
CM.image!(axes[1,2], ring_std, colormap=:afmhot); axes[1,2].title = "Ring Std. Dev."
CM.image!(axes[2,1], rast_mean, colormap=:afmhot); axes[2,1].title = "Rast Mean"
CM.image!(axes[2,2], rast_mean./rast_std, colormap=:afmhot, colorrange=(1.5, 3)); axes[2,2].title = "Rast SNR"
CM.hidedecorations!.(axes)
DisplayAs.Text(DisplayAs.PNG(fig)) #hide


# Finally, let's take a look at some of the ring parameters
figd = CM.Figure(;resolution=(650, 400));
p1 = CM.density(figd[1,1], rad2μas(chain.sky.r)*2, axis=(xlabel="Ring Diameter (μas)",))
p2 = CM.density(figd[1,2], rad2μas(chain.sky.σ)*2*sqrt(2*log(2)), axis=(xlabel="Ring FWHM (μas)",))
p3 = CM.density(figd[1,3], -rad2deg.(chain.sky.mp.:1) .+ 360.0, axis=(xlabel = "Ring PA (deg) E of N",))
p4 = CM.density(figd[2,1], 2*chain.sky.ma.:2, axis=(xlabel="Brightness asymmetry",))
p5 = CM.density(figd[2,2], 1 .- chain.sky.f, axis=(xlabel="Ring flux fraction",))
DisplayAs.Text(DisplayAs.PNG(figd)) #hide

# Now let's check the residuals using draws from the posterior
p = Plots.plot();
for s in sample(chain, 10)
    residual!(p, post, s, legend=false)
end
DisplayAs.Text(DisplayAs.PNG(p)) #hide

# And everything looks pretty good! Now comes the hard part: interpreting the results...
