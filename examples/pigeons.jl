using Pkg; Pkg.activate(@__DIR__)
using Comrade
using Distributions
using Plots
using StatsBase
using VLBIImagePriors
using Pigeons

# load eht-imaging we use this to load eht data
load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
obs = load_ehtim_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obs.add_scans()
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obs = scan_average(obs).add_fractional_noise(0.02)
# extract log closure amplitudes and closure phases
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)
damp = extract_amp(obs)


# For our model we will be using a rasterized image. This can be viewed as something like a
# non-parametric model. As a result of this we will need to use a `modelimage` object to
# store cache information we will need to compute the numerical FT.
function model(θ, metadata)
    (;c) = θ
    (; grid, cache, pulse) = metadata
    #Construct the image model
    img = IntensityMap(c, grid)
    cimg = ContinuousImage(img, pulse)
    #Create the modelimage object that will use a cache to compute the DFT
    return modelimage(cimg, cache)
end


nx = 8
ny = 8
fovxy = μas2rad(70.0)
# Now we can feed in the array information to form the cache. We will be using a DFT since
# it is efficient for so few pixels
pulse = BSplinePulse{3}()
buffer = IntensityMap(zeros(nx, ny), fovxy, fovxy)
cache = create_cache(DFTAlg(damp), buffer, pulse)
metadata = (;cache, pulse, grid=axiskeys(buffer))

# We will use a Dirichlet prior to enforce that the flux sums to unity since closures are
# degenerate to total flux.
prior = (c = ImageDirichlet(1.0, ny, nx),)

lklhd = RadioLikelihood(model, metadata, dlcamp, dcphase)
post = Posterior(lklhd, prior)

# Transform from simplex space to the unconstrained
tpost = asflat(post)
ℓ = logdensityof(tpost)

# Let's run an optimizer to get a nice starting location
# It turns out that gradients are really helpful here
ndim = dimension(tpost)

ℓ = logdensityof(tpost)


struct TransformedPr{P,T}
    p::P
    t::T
end

struct PigeonPost{P}
    p::P
end

TransformedPr(p::Comrade.NamedDist) = TransformedPr(p, asflat(p))
function (p::TransformedPr)(x)
    y, lj = Comrade.transform_and_logjac(p.t, x)
    return logdensityof(p.p, y) + lj
end

Pigeons.@provides target PigeonPost(p::Comrade.TransformedPosterior) = PigeonPost{typeof(p)}(p)

function (p::PigeonPost)(x)
    y = transform(p.p, x)
    return logdensityof(p.p.lpost, y)
end



Pigeons.create_state_initializer(target::PigeonPost, ::Inputs) = target
Pigeons.initialization(target::PigeonPost, rng::Pigeons.SplittableRandom, _::Int) = prior_sample(target.p)
Pigeons.create_explorer(::PigeonPost, ::Inputs) = Pigeons.SliceSampler()
Pigeons.create_reference_log_potential(target::PigeonPost, ::Inputs) = TransformedPr(target.p.lpost.prior)
function Pigeons.sample_iid!(pot::PigeonPost, replica)
    replica.state = initialization(pot, replica.rng, replica.replica_index)
end

function Pigeons.sample_iid!(pot::TransformedPr, replica)
    replica.state = inverse(pot.t, rand(pot.p))
end


res = pigeons(target = PigeonPost(tpost), n_rounds=5, recorder_builders=[index_process], n_chains=80, )

# Plot the mean image and standard deviation image
using StatsBase
samples = mms.(sample(hchain, 500))
imgs = intensitymap.(samples, fovx, fovy, 256, 256)

mimg, simg = mean_and_std(imgs)

 p1 = plot(mimg, title="Mean", clims=(0.0, maximum(mimg)))
p2 = plot(simg,  title="Std. Dev.", clims=(0.0, maximum(mimg)))
p2 = plot(simg./mimg,  title="Fractional Error", xlims=(-25.0,25.0), ylims=(-25.0,25.0))

# Computing information
# ```
# Julia Version 1.8.0
# Commit 742b9abb4d (2022-05-06 12:58 UTC)
# Platform Info:
#   OS: Linux (x86_64-pc-linux-gnu)
#   CPU: 11th Gen Intel(R) Core(TM) i7-1185G7 @ 3.00GHz
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-12.0.1 (ORCJIT, tigerlake)
# ```
