using Pkg; Pkg.activate(@__DIR__)
#Make sure you are using the main branch of Comrade
Pkg.add("Comrade#main")
using Comrade
using Distributions
# Add the interface package
Pkg.add("https://github.com/ptiede/Comrade.jl/lib/ComradePigeons")
using ComradePigeons
using Pigeons
using Plots
using StatsBase
using VLBIImagePriors

# load eht-imaging we use this to load eht data
load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obs.add_scans()
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obs = scan_average(obs).add_fractional_noise(0.02)
# extract log closure amplitudes and closure phases
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)

# Build the Model. Here we we a struct to hold some caches
# which will speed up imaging
# For our model we will be using a rasterized image. This can be viewed as something like a
# non-parametric model. As a result of this we will need to use a `modelimage` object to
# store cache information we will need to compute the numerical FT.
function model(θ, metadata)
    (;c) = θ
    (; cache, grid, pulse) = metadata
    #Construct the image model
    imap = IntensityMap(c, grid)
    cimg = ContinuousImage(imap, pulse)
    #Create the modelimage object that will use a cache to compute the DFT
    return modelimage(cimg, cache)
end

# Set up the grid
npix = 8
fovxy = μas2rad(72.0)
# Now we can feed in the array information to form the cache. We will be using a DFT since
# it is efficient for so few pixels
grid = imagepixels(fovxy, fovxy, npix, npix)
# We will use a Dirichlet prior to enforce that the flux sums to unity since closures are
# degenerate to total flux.
pulse = BSplinePulse{3}()
img = IntensityMap(zeros(npix,npix), grid)
cache = create_cache(DFTAlg(damp), img, pulse)
metadata = (;cache, img, grid, pulse)
prior = (c = ImageDirichlet(0.5, npix, npix),)

lklhd = RadioLikelihood(model, metadata, dlcamp, dcphase)
post = Posterior(lklhd, prior)

# Transform from simplex space to the unconstrained
tpost = asflat(post)
ℓ = logdensityof(tpost)


using OnlineStats
input = Pigeons.Inputs(tpost; recorder_builders=[Pigeons.index_process], seed=40, n_rounds=12, n_chains=50)

out = pigeons(input)

#Plot the image
ind = sortperm(logdensityof.(Ref(tpost), getproperty.(out.replicas, :state)))
plot(model(transform(tpost, out.replicas[ind[end]].state), metadata), xlims=(-32.5, 32.5), ylims=(-32.5, 32.5))
