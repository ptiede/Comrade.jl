# # Stokes I Simultaneous Image and Instrument Modeling

# In this tutorial, we will create a preliminary reconstruction of the 2017 M87 data on April 6
# by simultaneously creating an image and model for the instrument. By instrument model, we
# mean something akin to self-calibration in traditional VLBI imaging terminology. However,
# unlike traditional self-cal, we will at each point in our parameter space effectively explore
# the possible self-cal solutions. This will allow us to constrain and marginalize over the
# instrument effects, such as time variable gains.

# To get started we load Comrade.


using Pkg #hide
Pkg.activate(joinpath(@__DIR__, "../examples")) #hide
#-
using Comrade
using Pyehtim
using LinearAlgebra

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(42)

# Enzyme.Compiler.bitcode_replacement!(false)

# ## Load the Data


# To download the data visit https://doi.org/10.25739/g85n-f134
# First we will load our data:
obs = ehtim.obsdata.load_uvfits(joinpath(dirname(pathof(Comrade)), "..", "examples", "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))

obs = scan_average(obs.add_fractional_noise(0.01))

# Now we extract our complex visibilities.
dvis = extract_table(obs, ComplexVisibilities())



function sky(θ)
    (;radius, width, α1, β1, α2, β2, f, σG, τG, ξG, xG, yG) = θ
    α = (α1, α2)
    β = (β1, β2)
    ring = f*smoothed(stretched(MRing(α, β), radius, radius), width)
    g = (1-f)*shifted(rotated(stretched(Gaussian(), σG, σG*(1+τG)), ξG), xG, yG)
    return ring + g
end

function instrument(θ, metadata)
    (; lgamp, gphase) = θ
    (; gcache, gcachep) = metadata
    ## Now form our instrument model
    gvis = exp.(lgamp)
    gphase = exp.(1im.*gphase)
    jgamp = jonesStokes(gvis, gcache)
    jgphase = jonesStokes(gphase, gcachep)
    return JonesModel(jgamp*jgphase)
end

gcache = jonescache(dvis, ScanSeg())
gcachep = jonescache(dvis, ScanSeg(); autoref=SEFDReference((complex(1.0))))

using VLBIImagePriors
# Now we can form our metadata we need to fully define our model.
metadata = (;gcache, gcachep)

# We will also fix the total flux to be the observed value 1.1. This is because
# total flux is degenerate with a global shift in the gain amplitudes making the problem
# degenerate. To fix this we use the observed total flux as our value.

# Moving onto our prior, we first focus on the instrument model priors.
# Each sites requires its own prior on both the amplitudes and phases.
# For the amplitudes
# we assume that the gains are apriori well calibrated around unit gains (or 0 log gain amplitudes)
# which corresponds to no instrument corruption. The gain dispersion is then set to 10% for
# all sites except LMT, representing that we expect 10% deviations from scan-to-scan. For LMT
# we let the prior expand to 100% due to the known pointing issues LMT had in 2017.
using Distributions
using DistributionsAD
distamp = site_tuple(dvis, Normal(0.0, 0.1); LM = Normal(1.0))

distphase = site_tuple(dvis, DiagonalVonMises(0.0, inv(π^2)))



# We can now form our model parameter priors. Like our other imaging examples, we use a
# Dirichlet prior for our image pixels. For the log gain amplitudes, we use the `CalPrior`
# which automatically constructs the prior for the given jones cache `gcache`.
prior = NamedDist(
    radius = Uniform(μas2rad(10.0), μas2rad(30.0)),
    width = Uniform(μas2rad(1.0), μas2rad(10.0)),
    α1 = Uniform(-0.5, 0.5),
    β1 = Uniform(-0.5, 0.5),
    α2 = Uniform(-0.5, 0.5),
    β2 = Uniform(-0.5, 0.5),
    f = Uniform(0.0, 1.0),
    σG = Uniform(μas2rad(1.0), μas2rad(40.0)),
    τG = Uniform(0.0, 0.75),
    ξG = Uniform(0.0, 1π),
    xG = Uniform(-μas2rad(80.0), μas2rad(80.0)),
    yG = Uniform(-μas2rad(80.0), μas2rad(80.0)),
    lgamp = CalPrior(distamp, gcache),
    gphase = CalPrior(distphase, gcachep),
        )


# Putting it all together we form our likelihood and posterior objects for optimization and
# sampling.
lklhd = RadioLikelihood(sky, instrument, dvis; instrumentmeta=metadata)
post = Posterior(lklhd, prior)

# ## Reconstructing the Image and Instrument Effects

# To sample from this posterior, it is convenient to move from our constrained parameter space
# to an unconstrained one (i.e., the support of the transformed posterior is (-∞, ∞)). This is
# done using the `asflat` function.
tpost = asflat(post)
ndim = dimension(tpost)

# Our `Posterior` and `TransformedPosterior` objects satisfy the `LogDensityProblems` interface.
# This allows us to easily switch between different AD backends and many of Julia's statistical
# inference packages use this interface as well.
using Zygote
using Enzyme
# Enzyme.API.runtimeActivity!(true)

x0 = randn(rng, ndim)
ℓ = logdensityof(tpost)
gz, = Zygote.gradient(ℓ, x0)
dx0 = zero(x0)
autodiff(Reverse, Const(ℓ), Active, Duplicated(x0, dx0))
