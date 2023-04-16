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
#   - [`jonesT`](@ref) is the basis transform matrix $T$. This transformation is special and
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

using Comrade



# ## Load the Data
using Pkg #hide
Pkg.activate(joinpath(dirname(pathof(Comrade)), "..", "examples")) #hide

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(123)


# Now we will load some synthetic polarized data.
obs = load_ehtim_uvfits(joinpath(dirname(pathof(Comrade)), "..", "examples", "PolarizedExamples/polarized_gaussian_nogains_withdterms_withfr.uvfits"),
                        joinpath(dirname(pathof(Comrade)), "..", "examples", "PolarizedExamples/array.txt"))
# Notice that, unlike other non-polarized tutorials, we need to include a second argument.
# This is the **array file** of the observation and is required to determine the feed rotation
# of the array.

# Now we scan average the data since the data to boost the SNR and reduce the total data volume.
obs = scan_average(obs)
#-
# Now we extract our observed/corrupted coherency matrices.
dvis = extract_coherency(obs)

# ##Building the Model/Posterior


# To build the model, we first break it down into two parts:
#    1. **The image model**. In Comrade, all polarized image models are written in terms of the Stokes parameters.
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


function model(θ, metadata)
    (;c, f, p, angparams, dRx, dRy, dLx, dLy, lgp, gpp, lgr, gpr) = θ
    (; grid, cache, tcache, scancache, trackcache, phasecache) = metadata
    ## Construct the image model
    ## produce Stokes images from parameters
    imgI = f*c
    ## Converts from poincare sphere parameterization of polzarization to Stokes Parameters
    pimg = PoincareSphere2Map(imgI, p, angparams, grid)
    m = ContinuousImage(pimg, cache)

    ## Now construct the basis transformation cache
    jT = jonesT(tcache)

    ## Gain product parameters
    gPa = exp.(lgp/2 .+ 0im)
    gPp = exp.(1im.*gpp/2)
    Gpa = jonesG(gPa, gPa, scancache)
    Gpp = jonesG(gPp, gPp, phasecache)
    ## Gain ratio
    gR = exp.(lgr/2 .+ 1im.*gpr/2)
    Gr = jonesG(gR, inv.(gR), trackcache)
    ##D-terms
    D = jonesD(complex.(dRx, dRy), complex.(dLx, dLy), trackcache)
    ## sandwich all the jones matrices together
    J = Gpa*Gpp*Gr*D*jT
    ## form the complete Jones or RIME model. We use tcache here
    ## to set the reference basis of the model.
    return JonesModel(J, m, tcache)
end

# Now, we define the model metadata required to build the model.
# We specify our image grid and cache model needed to define the polarimetric
# image model.
fovx = μas2rad(50.0)
fovy = μas2rad(50.0)
nx = 5
ny = floor(Int, fovy/fovx*nx)
grid = imagepixels(fovx, fovy, nx, ny) # image grid
buffer = IntensityMap(zeros(nx, ny), grid) # buffer to store temporary image
pulse = BSplinePulse{3}() # pulse we will be using
cache = create_cache(NFFTAlg(dvis), buffer, pulse) # cache to define the NFFT transform

# To define the instrument models, $T$, $G$, $D$, we need to build some Jones caches (see [`JonesCache`](@ref)) that map from a flat
# vector of gain/dterms to the specific sites for each baseline.
#
# First, we will define our deterministic transform cache. Note that this dataset has need
# been pre-corrected for feed rotation, so we need to add those into the `tcache`.
tcache = TransformCache(dvis; add_fr=true, ehtim_fr_convention=false)
#-
# Next we define our cache that maps quantities e.g., gain products, that change from scan-to-scan.
scancache = jonescache(dvis, ScanSeg())

# In addition we will assign a reference station. This is necessary for gain phases due to a trivial degeneracy being present.
# To do this we will select ALMA `AA` as the reference station as is standard in EHT analyses.
phase_segs = station_tuple(dvis, ScanSeg(); AA=FixedSeg(1.0 + 0.0im))
phasecache = jonescache(dvis, phase_segs)

#-
# Finally, we define our cache that maps quantities, e.g., gain ratios and d-terms, that are constant
# across a observation night, and we collect everything together.
trackcache = jonescache(dvis, TrackSeg())
metadata = (;cache, grid, tcache, scancache, trackcache, phasecache)

# Moving onto our prior, we first focus on the instrument model priors.
# Each station gain requires its own prior on both the amplitudes and phases.
# For the amplitudes, we assume that the gains are apriori well calibrated around unit gains (or 0 log gain amplitudes)
# which corresponds to no instrument corruption. The gain dispersion is then set to 10% for
# all stations except LMT, representing that we expect 10% deviations from scan-to-scan. For LMT,
# we let the prior expand to 100% due to the known pointing issues LMT had in 2017.

using Distributions
using DistributionsAD
distamp = (AA = Normal(0.0, 0.1),
           AP = Normal(0.0, 0.1),
           LM = Normal(0.0, 0.1),
           AZ = Normal(0.0, 0.1),
           JC = Normal(0.0, 0.1),
           PV = Normal(0.0, 0.1),
           SM = Normal(0.0, 0.1),
           )
#-
# For the phases, we assume that the atmosphere effectively scrambles the gains.
# Since the gain phases are periodic, we also use broad von Mises priors for all stations.
# Notice that we don't assign a prior for AA since we have already fixed it.
using VLBIImagePriors
distphase = (
             AP = DiagonalVonMises(0.0, inv(π^2)),
             LM = DiagonalVonMises(0.0, inv(π^2)),
             AZ = DiagonalVonMises(0.0, inv(π^2)),
             JC = DiagonalVonMises(0.0, inv(π^2)),
             PV = DiagonalVonMises(0.0, inv(π^2)),
             SM = DiagonalVonMises(0.0, inv(π^2)),
           )
#-
# However, we can now also use a little additional information about the phase offsets
# where in most cases, they are much better behaved than the products
distphase_ratio = (
             AA = DiagonalVonMises(0.0, inv(0.01)),
             AP = DiagonalVonMises(0.0, inv(0.1^2)),
             LM = DiagonalVonMises(0.0, inv(0.1^2)),
             AZ = DiagonalVonMises(0.0, inv(0.1^2)),
             JC = DiagonalVonMises(0.0, inv(0.1^2)),
             PV = DiagonalVonMises(0.0, inv(0.1^2)),
             SM = DiagonalVonMises(0.0, inv(0.1^2)),
           )


# Moving onto the d-terms, here we directly parameterize the real and complex components
# of the d-terms since they are expected to be complex numbers near the origin. To help enforce
# this smallness, a weakly informative Normal prior is used.
distD = ( AA = Normal(0.0, 0.1),
          AP = Normal(0.0, 0.1),
          LM = Normal(0.0, 0.1),
          AZ = Normal(0.0, 0.1),
          JC = Normal(0.0, 0.1),
          PV = Normal(0.0, 0.1),
          SM = Normal(0.0, 0.1),
        )



# Our image priors are:
#   - We use a Dirichlet prior, `ImageDirichlet`, with unit concentration for our stokes I image pixels, `c`.
#   - For the total polarization fraction, `p`, we assume an uncorrelated uniform prior `ImageUniform` for each pixel.
#   - To specify the orientation of the polarization, `angparams`, on the Poincare sphere,
#     we use a uniform spherical distribution, `ImageSphericalUniform`.
#-
# For all the calibration parameters, we use a helper function `CalPrior` which builds the
# prior given the named tuple of station priors and a `JonesCache`
# that specifies the segmentation scheme. For the gain products, we use the `scancache`, while
# for every other quantity, we use the `trackcache`.
prior = (
          c = ImageDirichlet(1.0, nx, ny),
          f = Uniform(0.7, 1.2),
          p = ImageUniform(nx, ny),
          angparams = ImageSphericalUniform(nx, ny),
          dRx = CalPrior(distD, trackcache),
          dRy = CalPrior(distD, trackcache),
          dLx = CalPrior(distD, trackcache),
          dLy = CalPrior(distD, trackcache),
          lgp = CalPrior(distamp, scancache),
          gpp = CalPrior(distphase, phasecache),
          lgr = CalPrior(distamp, trackcache),
          gpr = CalPrior(distphase_ratio,trackcache),
          )


# Putting it all together, we form our likelihood and posterior objects for optimization and
# sampling.
lklhd = RadioLikelihood(model, metadata, dvis)
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
using OptimizationOptimJL
using Zygote
f = OptimizationFunction(tpost, Optimization.AutoZygote())
ℓ = logdensityof(tpost)
prob = Optimization.OptimizationProblem(f, prior_sample(rng, tpost), nothing)
sol = solve(prob, LBFGS(), maxiters=15_000, callback=((x,p)->(@info ℓ(x);false)), g_tol=1e-1)

# !!! warning
#     Fitting polarized images is generally much harder than Stokes I imaging. This difficulty means that
#     optimization can take a long time, and starting from a good starting location
#     is often required.
#-
# Before we analyze our solution, we need to transform it back to parameter space.
xopt = transform(tpost, sol)

# Now let's evaluate our fits by plotting the residuals
using Plots
residual(model(xopt, metadata), dvis)

# These look reasonable, although there may be some minor overfitting.
# Let's compare our results to the ground truth values we know in this example.
# First, we will load the polarized truth
using AxisKeys
imgtrue = Comrade.load(joinpath(dirname(pathof(Comrade)), "..", "examples", "PolarizedExamples/polarized_gaussian.fits"), StokesIntensityMap)
# Select a reasonable zoom in of the image.
imgtruesub = imgtrue(Interval(-fovx/2, fovx/2), Interval(-fovy/2, fovy/2))
# plot(imgtruesub, title="True Image", xlims=(-25.0,25.0), ylims=(-25.0,25.0))
#-
img = intensitymap!(copy(imgtruesub), model(xopt, metadata))
# plot(img, title="Reconstructed Image", xlims=(-25.0,25.0), ylims=(-25.0,25.0))

# Let's compare some image statics, like the total linear polarization fraction
using Comrade.ComradeBase: linearpol
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
gphase_ratio = caltable(trackcache, xopt.gpr)
#-
# we see that they are all very small. Which should be the case since this data doesn't have gain corruptions!
# Similarly our gain ratio amplitudes are also very close to unity as expected.
gamp_ratio   = caltable(trackcache, exp.(xopt.lgr))
#-
# Plotting the gain phases, we see some offsets from zero. This is because the prior on the gain product
# phases is very broad, so we can't phase center the image. For realistic data
# this is always the case since the atmosphere effectively scrambles the phases.
gphase_prod = caltable(phasecache, xopt.gpp)
# plot(gphase_prod, layout=(3,3), size=(650,500))
#-
# Finally, the product gain amplitudes are all very close to unity as well, as expected since gain corruptions
# have not been added to the data.
gamp_prod = caltable(scancache, exp.(xopt.lgp))
# plot(gamp_prod, layout=(3,3), size=(650,500))
#-
# At this point, you should run the sampler to recover an uncertainty estimate,
# which is identical to every other imaging example
# (see, e.g., [Stokes I Simultaneous Image and Instrument Modeling](@ref). However,
# due to the time it takes to sample, we will skip that for this tutorial. Note that on the computer environment listed
# below, 20_000 MCMC steps take 4 hours.


# [^1]: Hamaker J.P, Bregman J.D., Sault R.J. (1996) [https://articles.adsabs.harvard.edu/pdf/1996A%26AS..117..137H]
# [^2]: Pesce D. (2021) [https://ui.adsabs.harvard.edu/abs/2021AJ....161..178P/abstract]
#-

# ## Computing information
# ```
# Julia Version 1.8.5
# Commit 17cfb8e65ea (2023-01-08 06:45 UTC)
# Platform Info:
#   OS: Linux (x86_64-linux-gnu)
#   CPU: 32 × AMD Ryzen 9 7950X 16-Core Processor
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-13.0.1 (ORCJIT, znver3)
#   Threads: 1 on 32 virtual cores
# Environment:
#   JULIA_EDITOR = code
#   JULIA_NUM_THREADS = 1
# ```
