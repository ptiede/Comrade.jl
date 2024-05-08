# # Stokes I Simultaneous Image and Instrument Modeling

# In this tutorial, we will create a preliminary reconstruction of the 2017 M87 data on April 6
# by simultaneously creating an image and model for the instrument. By instrument model, we
# mean something akin to self-calibration in traditional VLBI imaging terminology. However,
# unlike traditional self-cal, we will at each point in our parameter space effectively explore
# the possible self-cal solutions. This will allow us to constrain and marginalize over the
# instrument effects, such as time variable gains.

# To get started we load Comrade.


using Pkg #hide
Pkg.activate(joinpath(@__DIR__, ".")) #hide
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
obs = ehtim.obsdata.load_uvfits(joinpath(dirname(pathof(Comrade)), "..", "examples", "Data", "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))

obs = scan_average(obs.add_fractional_noise(0.01))

# Now we extract our complex visibilities.
dvis = extract_table(obs, ComplexVisibilities())



function sky(θ, metadata)
    c = θ.c
    (;cache) = metadata
    m = ContinuousImage(c, cache)
    return m
end

npix = 13
fovx = μas2rad(80.0)
fovy = μas2rad(80.0)

grid = imagepixels(fovx, fovy, npix, npix)
cache = create_cache(DFTAlg(dvis), grid, BSplinePulse{3}())



using VLBIImagePriors

skymetadata = (;grid, cache)

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

using Distributions
using DistributionsAD
distamp = site_tuple(dvis, Normal(0.0, 0.1); LM = Normal(1.0))

distphase = site_tuple(dvis, DiagonalVonMises(0.0, inv(π^2)))


prior = NamedDist(
            (c = ImageDirichlet(1.0, npix, npix),
            lgamp = CalPrior(distamp, gcache),
            gphase = CalPrior(distphase, gcachep),)
            )



lklhd = RadioLikelihood(sky, instrument, dvis; skymeta=skymetadata, instrumentmeta=metadata)
post = Posterior(lklhd, prior)

tpost = asflat(post)
ndim = dimension(tpost)

using Enzyme
using Zygote
Enzyme.API.runtimeActivity!(true)
# Enzyme.Compiler.bitcode_replacement!(false)

# Enzyme.API.printall!(false)
x0 = randn(ndim)
dx0 = zero(x0)
lt=logdensityof(tpost)
ℓ = logdensityof(tpost, x0)
gz,  = Zygote.gradient(lt, x0)
using BenchmarkTools
autodiff(Reverse, Const(lt), Active, Duplicated(x0, fill!(dx0, 0.0)))
@benchmark autodiff($Reverse, $logdensityof, $(Const(tpost)), Duplicated($x0, fill!($dx0, 0.0)))
# @benchmark Zygote.gradient($lt, $x0)
