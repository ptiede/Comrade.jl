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



# ## Load the Data


# To download the data visit https://doi.org/10.25739/g85n-f134
# First we will load our data:
obs = ehtim.obsdata.load_uvfits(joinpath(dirname(pathof(Comrade)), "..", "examples", "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))

obs = scan_average(obs.add_fractional_noise(0.01))

# Now we extract our complex visibilities.
dvis = extract_table(obs, ComplexVisibilities())



function sky(θ, metadata)
    c = θ
    (;grid, cache) = metadata
    img = IntensityMap(c, grid)
    m = ContinuousImage(img, cache)
    return m
end

npix = 12
fovx = μas2rad(80.0)
fovy = μas2rad(80.0)

grid = imagepixels(fovx, fovy, npix, npix)
buffer = IntensityMap(zeros(npix, npix), grid)
cache = create_cache(NFFTAlg(dvis), buffer, DeltaPulse())



using VLBIImagePriors
prior = ImageUniform(npix, npix)

skymetadata = (;grid, cache)
# instrumentmetadata = (;gcache, gcachep)
lklhd = RadioLikelihood(sky, dvis; skymeta=skymetadata)
post = Posterior(lklhd, prior)

tpost = asflat(post)
ndim = dimension(tpost)

using Enzyme
Enzyme.API.runtimeActivity!(true)
# Enzyme.API.printall!(false)
x0 = prior_sample(tpost)
dx0 = zero(x0)
lt=logdensityof(tpost)
ℓ = logdensityof(tpost, x0)
dtpost = deepcopy(tpost)
using BenchmarkTools
@time autodiff(Reverse, logdensityof, Const(tpost), Duplicated(x0, fill!(dx0, 0.0)))
@btime autodiff(Reverse, logdensityof, $(Duplicated(tpost, dtpost)), Duplicated($x0, fill!($dx0, 0.0)))
