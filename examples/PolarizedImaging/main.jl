# # Polarized Image and Instrumental Modeling

# In this tutorial, we will analyze a simulated simple polarized dataset to demonstrate
# Comrade's polarized imaging capabilities.

# ## Introduction to Polarized Imaging
# The EHT is a polarized interferometer. However, like all VLBI interferometers, it does not
# directly measure the Stokes parameters (I, Q, U, V). Instead, it measures components
# related to the electric field at the telescope along two *directions* using feeds.
# There are two types of feeds at telescopes: circular, which measure $R/L$ components of the
# electric field, and linear feeds, which measure $X/Y$ components of the electric field.
# Most sites in the
# EHT use circular feeds, meaning they measure the right (R) and left electric field (L) at each telescope.
# These circular electric field measurements are then correlated, producing **coherency matrices**,
#
# ```math
#  C_{ij} = \begin{pmatrix}
#        RR^* &  RL^*\\
#        LR^* &  LL^*
#      \end{pmatrix}.
# ```
#
# These coherency matrices are the fundamental object in interferometry and what
# the telescope observes. For a perfect interferometer, these coherency matrices are related to
# the usual Fourier transform of the stokes parameters by
#
# ```math
#   \begin{pmatrix}
#       \tilde{I}\\ \tilde{Q} \\ \tilde{U} \\ \tilde{V}
#   \end{pmatrix}
#   =\frac{1}{2}
#   \begin{pmatrix}
#      RR^* + LL^* \\
#      RL^* + LR^* \\
#      i(LR^* - RL^*)\\
#      RR^* - LL^*
#   \end{pmatrix},
# ```
# for circularly polarized measurements.
#-
# !!! note
#     In this tutorial, we stick to circular feeds but Comrade has the capabilities
#     to model linear (XX,XY, ...) and mixed basis coherencies (e.g., RX, RY, ...).
#-
#
# In reality, the measure coherencies are corrupted by both the atmosphere and the
# telescope itself. In `Comrade` we use the RIME formalism [^1] to represent these corruptions,
# namely our measured coherency matrices $V_{ij}$ are given by
# ```math
#    V_{ij} = J_iC_{ij}J_j^\dagger
# ```
# where $J$ is known as a *Jones matrix* and $ij$ denotes the baseline $ij$ with sites $i$ and $j$.
#-
# `Comrade` is highly flexible with how the Jones matrices are formed and provides several
# convenience functions that parameterize standard Jones matrices. These matrices include:
#   - [`jonesG`](@ref) which builds the set of complex gain Jones matrices
# ```math
#   G = \begin{pmatrix}
#           g_a   &0\\
#           0     &g_b\\
#       \end{pmatrix}
# ```
#   - [`jonesD`](@ref) which builds the set of complex d-terms Jones matrices
# ```math
#   D = \begin{pmatrix}
#           1   & d_a\\
#           d_b     &1\\
#       \end{pmatrix}
# ```
#   - [`jonesR`](@ref) is the basis transform matrix $T$. This transformation is special and
#      combines two things using the decomposition $T=FB$. The first, $B$, is the transformation from
#      some reference basis to the observed coherency basis (this allows for mixed basis measurements).
#      The second is the feed rotation, $F$, that transforms from some reference axis to the axis of the
#      telescope as the source moves in the sky. The feed rotation matrix `F` in terms of
#      the per station feed rotation angle $\varphi$ is
# ```math
#   F = \begin{pmatrix}
#           e^{-i\varphi}   & 0\\
#           0     & e^{i\varphi}\\
#       \end{pmatrix}
# ```
#-
#
#  In the rest of the tutorial, we are going to solve for all of these instrument model terms on
#  in addition to our image structure to reconstruct a polarized image of a synthetic dataset.

import Pkg #hide
__DIR = @__DIR__ #hide
pkg_io = open(joinpath(__DIR, "pkg.log"), "w") #hide
Pkg.activate(__DIR; io=pkg_io) #hide
Pkg.develop(; path=joinpath(__DIR, "..", ".."), io=pkg_io) #hide
Pkg.instantiate(; io=pkg_io) #hide
Pkg.precompile(; io=pkg_io) #hide
close(pkg_io) #hide


# ## Load the Data
# To get started we will load Comrade
ENV["GKSwstype"] = "nul" #hide
using Comrade

# ## Load the Data
using Pyehtim

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(14)


# Now we will load some synthetic polarized data.
obs = Pyehtim.load_uvfits_and_array(
        joinpath(__DIR, "../Data", "polarized_gaussian_all_corruptions.uvfits"),
        joinpath(__DIR, "../Data", "array.txt"), polrep="circ")


# Notice that, unlike other non-polarized tutorials, we need to include a second argument.
# This is the **array file** of the observation and is required to determine the feed rotation
# of the array.

# Now we scan average the data since the data to boost the SNR and reduce the total data volume.
obs = scan_average(obs)
#-
# Now we extract our observed/corrupted coherency matrices.
dvis = extract_table(obs, Coherencies())

# ##Building the Model/Posterior


# To build the model, we first break it down into two parts:
#    1. **The image or sky model**. In Comrade, all polarized image models are written in terms of the Stokes parameters.
#       The reason for using Stokes parameters is that it is usually what physical models consider and is
#       the often easiest to reason about since they are additive. In this tutorial, we will use a polarized image model based on Pesce (2021)[^2].
#       This model parameterizes the polarized image in terms of the [`Poincare sphere`](https://en.wikipedia.org/wiki/Unpolarized_light#Poincar%C3%A9_sphere),
#       and allows us to easily incorporate physical restrictions such as $I^2 ≥ Q^2 + U^2 + V^2$.
#    2. **The instrument model**. The instrument model specifies the model that describes the impact of instrumental and atmospheric effects.
#       We will be using the $J = GDT$ decomposition we described above. However, to parameterize the
#       R/L complex gains, we will be using a gain product and ratio decomposition. The reason for this decomposition
#       is that in realistic measurements, the gain ratios and products have different temporal characteristics.
#       Namely, many of the EHT observations tend to demonstrate constant R/L gain ratios across an
#       nights observations, compared to the gain products, which vary every scan. Additionally, the gain ratios tend to be smaller (i.e., closer to unity) than the gain products.
#       Using this apriori knowledge, we can build this into our model and reduce
#       the total number of parameters we need to model.

using StatsFuns: logistic
function sky(θ, metadata)
    (;c, σ, p, p0, pσ, angparams) = θ
    (;cache) = metadata
    cp = σ.*c.params
    rast = to_simplex(CenteredLR(), cp)
    pim = logistic.(p0 .+ pσ.*p.params)
    m = PoincareSphere2Map(rast, pim, angparams, cache)
    x0, y0 = centroid(stokes(m, :I))
    ms = shifted(m, -x0, -y0)
    return ms
end



function instrument(θ, metadata)
    (; lgR, lgL, gpR, gpL, dRx, dRy, dLx, dLy) = θ
    (; rcache, scancache, phasecache1, phasecache2, trackcache) = metadata
    ## Now construct the basis transformation cache
    jT = jonesR(rcache)

    ## Gain product parameters
    gRa = exp.(lgR .+ 0.0im)
    gLa = exp.(lgL .+ lgR .+ 0.0im) # directly model hte LR offset
    Ga = jonesG(gRa, gLa, scancache)

    gRp = exp.(1im*gpR)
    gLp = exp.(1im*(gpL))
    Gp1 = jonesG(gRp, gRp, phasecache1)
    Gp2 = jonesG(fill(complex(1.0), length(gLp)), gLp, phasecache2)

    ##D-terms
    D = jonesD(complex.(dRx, dRy), complex.(dLx, dLy), trackcache)
    ## sandwich all the jones matrices together
    J = map(Ga, Gp1, Gp2, D, jT) do ga, gp1, gp2, d, jt
        ga*gp1*gp2*d*jt
    end

    return JonesModel(J, rcache)
end

# Now, we define the model metadata required to build the model.
# We specify our image grid and cache model needed to define the polarimetric
# image model.
fovx = μas2rad(50.0)
fovy = μas2rad(50.0)
nx = 12
ny = floor(Int, fovy/fovx*nx)
grid = imagepixels(fovx, fovy, nx, ny) # image grid
pulse = BSplinePulse{3}() # pulse we will be using
cache = create_cache(NFFTAlg(dvis), grid, pulse) # cache to define the NFFT transform
skymeta = (;cache)


instrumentmeta = (;rcache, scancache, phasecache1, phasecache2, trackcache)
lklhd = RadioLikelihood(sky, instrument, dvis;
                        skymeta, instrumentmeta)




# To define the instrument models, $T$, $G$, $D$, we need to build some Jones caches (see [`JonesCache`](@ref)) that map from a flat
# vector of gain/dterms to the specific sites for each baseline.
#
# First, we will define our deterministic response cache. This define how the telescope reponds
# to a signal, including the impact of telescope field rotation due to parallactic and elevation
# changes.
rcache = ResponseCache(dvis)
#-
# Next we define our cache that maps quantities e.g., gain products, that change from scan-to-scan.
scancache = jonescache(dvis, ScanSeg())

# In addition we will assign a reference station. This is necessary for gain phases due to a trivial degeneracy being present.
# To do this we will select ALMA `AA` as the reference station as is standard in EHT analyses.
phasecache1 = jonescache(dvis, ScanSeg(); autoref = SEFDReference(complex(1.0)))
phasecache2 = jonescache(dvis, ScanSeg(); autoref = SingleReference(:AA, complex(1.0)))

#-
# Finally, we define our cache that maps quantities, e.g., gain ratios and d-terms, that are constant
# across a observation night, and we collect everything together.
trackcache = jonescache(dvis, TrackSeg())
instrumentmeta = (;rcache, scancache, phasecache1, phasecache2, trackcache)

# Moving onto our prior, we first focus on the instrument model priors.
# Each station gain requires its own prior on both the amplitudes and phases.
# For the amplitudes, we assume that the gains are apriori well calibrated around unit gains (or 0 log gain amplitudes)
# which corresponds to no instrument corruption. The gain dispersion is then set to 10% for
# all stations except LMT, representing that we expect 10% deviations from scan-to-scan. For LMT,
# we let the prior expand to 100% due to the known pointing issues LMT had in 2017.

using Distributions
using DistributionsAD
using VLBIImagePriors
distamp = station_tuple(dvis, Normal(0.0, 0.1))
distamp_ratio = station_tuple(dvis, Normal(0.0, 0.01))

#-
# For the phases, we assume that the atmosphere effectively scrambles the gains.
# Since the gain phases are periodic, we also use broad von Mises priors for all stations.
# Notice that we don't assign a prior for AA since we have already fixed it.
distphase = station_tuple(dvis, DiagonalVonMises(0.0, inv(π^2)); reference=:AA)

#-
# However, we can now also use a little additional information about the phase offsets
# where in most cases, they are much better behaved than the products
distphase_ratio = station_tuple(dvis, DiagonalVonMises(0.0, inv(0.1)); reference=:AA)


# Moving onto the d-terms, here we directly parameterize the real and complex components
# of the d-terms since they are expected to be complex numbers near the origin. To help enforce
# this smallness, a weakly informative Normal prior is used.
dist_leak = station_tuple(dvis, Normal(0.0, 0.1))


# Our image priors are:
#   - We use a Dirichlet prior, `ImageDirichlet`, with unit concentration for our stokes I image pixels, `c`.
#   - For the total polarization fraction, `p`, we assume an uncorrelated uniform prior `ImageUniform` for each pixel.
#   - To specify the orientation of the polarization, `angparams`, on the Poincare sphere,
#     we use a uniform spherical distribution, `ImageSphericalUniform`.
#-
rat = beamsize(dvis)/step(grid.X)
cmarkov = ConditionalMarkov(GMRF, grid; order=1)
dρ = truncated(InverseGamma(2.0, -log(0.1)*rat); lower=2.0, upper=max(nx, ny))
cprior = HierarchicalPrior(cmarkov, dρ)


# For all the calibration parameters, we use a helper function `CalPrior` which builds the
# prior given the named tuple of station priors and a `JonesCache`
# that specifies the segmentation scheme. For the gain products, we use the `scancache`, while
# for every other quantity, we use the `trackcache`.
fwhmfac = 2.0*sqrt(2.0*log(2.0))
prior = NamedDist(
    c = cprior,
    σ  = truncated(Normal(0.0, 0.1); lower=0.0),
    p  = cprior,
    p0 = Normal(-2.0, 1.0),
    pσ =  truncated(Normal(0.0, 1.0); lower=0.01),
    angparams = ImageSphericalUniform(nx, ny),
    lgR  = CalPrior(distamp, scancache),
    gpR  = CalPrior(distphase, phasecache1),
    lgL  = CalPrior(distamp_ratio, scancache),
    gpL  = CalPrior(distphase_ratio, phasecache2),
    dRx = CalPrior(dist_leak, trackcache),
    dRy = CalPrior(dist_leak, trackcache),
    dLx = CalPrior(dist_leak, trackcache),
    dLy = CalPrior(dist_leak, trackcache),
    )


# Putting it all together, we form our likelihood and posterior objects for optimization and
# sampling.
lklhd = RadioLikelihood(sky, instrument, dvis; skymeta, instrumentmeta)
post = Posterior(lklhd, prior)

# ## Reconstructing the Image and Instrument Effects

# To sample from this posterior, it is convenient to move from our constrained parameter space
# to an unconstrained one (i.e., the support of the transformed posterior is (-∞, ∞)). This transformation is
# done using the `asflat` function.
tpost = asflat(post)

# We can now also find the dimension of our posterior or the number of parameters we will sample.
# !!! warning
#     This can often be different from what you would expect. This difference is especially true when using
#     angular variables, where we often artificially increase the dimension
#     of the parameter space to make sampling easier.
#-

ndim = dimension(tpost)


# Now we optimize. Unlike other imaging examples, we move straight to gradient optimizers
# due to the higher dimension of the space.
using ComradeOptimization
using OptimizationOptimisers
using Zygote
Comrade.Enzyme.API.runtimeActivity!(true)
f = OptimizationFunction(tpost, Optimization.AutoZygote())
ℓ = logdensityof(tpost)
prob = Optimization.OptimizationProblem(f, rand(ndim) .- 0.5, nothing)
sol = solve(prob, OptimizationOptimisers.Adam(), maxiters=1_000, g_tol=1e-1);

# !!! warning
#     Fitting polarized images is generally much harder than Stokes I imaging. This difficulty means that
#     optimization can take a long time, and starting from a good starting location
#     is often required.
#-
# Before we analyze our solution, we need to transform it back to parameter space.
xopt = transform(tpost, sol.u)

# Now let's evaluate our fits by plotting the residuals
using Plots
residual(vlbimodel(post, xopt), dvis)

# These look reasonable, although there may be some minor overfitting.
# Let's compare our results to the ground truth values we know in this example.
# First, we will load the polarized truth
imgtrue = Comrade.load(joinpath(__DIR, "..", "Data", "polarized_gaussian.fits"), StokesIntensityMap)
# Select a reasonable zoom in of the image.
imgtruesub = regrid(imgtrue, imagepixels(fovx, fovy, nx*4, ny*4))
img = intensitymap!(copy(imgtruesub), skymodel(post, xopt))

#Plotting the results gives
import CairoMakie as CM
CM.activate!(type = "png", px_per_unit=3) #hide
fig = CM.Figure(;resolution=(450, 350));
polimage(fig[1,1], imgtruesub,
                   axis=(xreversed=true, aspect=1, title="Truth"),
                   nvec = 8,
                   length_norm=1/2, plot_total=true, pcolormap=:RdBu,
                   pcolorrange=(-0.25, 0.25),)
polimage(fig[1,2], img,
                   axis=(xreversed=true, aspect=1, title="Recon.",),
                   nvec = 8,
                   length_norm=1/2, plot_total=true, pcolormap=:RdBu,
                   pcolorrange=(-0.25, 0.25),)
CM.Colorbar(fig[2,:], colormap=:RdBu, vertical=false, colorrange=(-0.25, 0.25), label="Signed Polarization Fraction sign(V)*|p|", flipaxis=false)
CM.colgap!(fig.layout, 3)
CM.rowgap!(fig.layout, 3)
CM.hidedecorations!.(fig.content[1:2])
fig
#-

# Let's compare some image statics, like the total linear polarization fraction
ftrue = flux(imgtruesub);
@info "Linear polarization true image: $(abs(linearpol(ftrue))/ftrue.I)"
frecon = flux(img);
@info "Linear polarization recon image: $(abs(linearpol(frecon))/frecon.I)"

# And the Circular polarization fraction
@info "Circular polarization true image: $(ftrue.V/ftrue.I)"
@info "Circular polarization recon image: $(frecon.V/frecon.I)"

# Because we also fit the instrument model, we can inspect their parameters.
# To do this, `Comrade` provides a `caltable` function that converts the flattened gain parameters
# to a tabular format based on the time and its segmentation.
dR = caltable(trackcache, complex.(xopt.dRx, xopt.dRy))

# We can compare this to the ground truth d-terms
#
# |      time |            AA   |           AP   |         AZ   |        JC    |        LM    |       PV     |      SM |
# |-----------|-----------------|----------------|--------------|--------------|--------------|--------------|---------|
# | 0.0       | 0.01-0.02im      | -0.08+0.07im  |  0.09-0.10im | -0.04+0.05im |  0.03-0.02im | -0.01+0.02im | 0.08-0.07im |
#

# And same for the left-handed dterms
#
dL = caltable(trackcache, complex.(xopt.dLx, xopt.dLy))
#
# |      time |            AA   |           AP   |         AZ   |        JC    |        LM    |       PV     |      SM |
# |-----------|-----------------|----------------|--------------|--------------|--------------|--------------|---------|
# | 0.0       | 0.03-0.04im      | -0.06+0.05im  |  0.09-0.08im | -0.06+0.07im |  0.01-0.00im | -0.03+0.04im | 0.06-0.05im |
#


# Looking at the gain phase ratio
gphase_ratio = caltable(phasecache2, xopt.gpL)
#-
# we see that they are all very small. Which should be the case since this data doesn't have gain corruptions!
# Similarly our gain ratio amplitudes are also very close to unity as expected.
gamp_ratio   = caltable(scancache, exp.(xopt.lgL))
#-
# Plotting the gain phases, we see some offsets from zero. This is because the prior on the gain product
# phases is very broad, so we can't phase center the image. For realistic data
# this is always the case since the atmosphere effectively scrambles the phases.
gphase_prod = caltable(phasecache1, xopt.gpR)
plot(gphase_prod, layout=(3,3), size=(650,500))
plot!(gphase_ratio, layout=(3,3), size=(650,500))
#-
# Finally, the product gain amplitudes are all very close to unity as well, as expected since gain corruptions
# have not been added to the data.
gamp_prod = caltable(scancache, exp.(xopt.lgR))
plot(gamp_prod, layout=(3,3), size=(650,500))
plot!(gamp_ratio, layout=(3,3), size=(650,500))
#-
# At this point, you should run the sampler to recover an uncertainty estimate,
# which is identical to every other imaging example
# (see, e.g., [Stokes I Simultaneous Image and Instrument Modeling](@ref)). However,
# due to the time it takes to sample, we will skip that for this tutorial.
using ComradeAHMC
metric = DiagEuclideanMetric(ndim)
chain = sample(rng, tpost, AHMC(;metric), 700; n_adapts=500, initial_params=sol.u, progress=true)

ms = skymodel.(Ref(post), chain[500:5:end])
imgs = intensitymap.(ms, Ref(grid))
imageviz(mean(imgs), plot_total=false, nvec=10, colormap=:afmhot, adjust_length=true)

# [^1]: Hamaker J.P, Bregman J.D., Sault R.J. (1996) [https://articles.adsabs.harvard.edu/pdf/1996A%26AS..117..137H]
# [^2]: Pesce D. (2021) [https://ui.adsabs.harvard.edu/abs/2021AJ....161..178P/abstract]
#-
