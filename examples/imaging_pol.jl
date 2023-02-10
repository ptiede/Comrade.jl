# # Polarized Image and Instrumental Modeling

# In this tutorial we will analyze an simulated simple polarized dataset to demonstrate
# Comrade's polarized imaging capabilities.

# ## Introduction to Polarized Imaging
# The EHT is a polarized interferometer. However, like all VLBI interferometers it does not
# directly measure the Stokes parameters (I, Q, U, V). Instead it measures components
# related to the electric field at the telescope along two *directions* using feeds. The EHT uses
# a circular feed telescope meaning that is measures the right and left. Most sites in the
# EHT use circular feeds meaning they measure the right (R) and left electric field (L) at each telescope.
# These products are then correlated producing **coherency matrices**,
#
# ```math
#  C_{ij} = \begin{pmatrix}
#        RR^* &  RL^*\\
#        LR^* &  LL^*
#      \end{pmatrix}.
# ```
#
# These coherency matrices are the fundamental object in interferometetry and what
# we actually observe. For a perfect interferometer these coherency matrices are related to
# usual Fourier transform of the stokes parameters by
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
#   \end{pmatrix}.
# ```
# Note however, that Comrade has the capabilities to model Linear (XX,XY, ...) and even mixed basis
# polarization (e.g., RX, RY, ...) as well.
#
# In reality the measure coherencies are corrupted by both the atmosphere and the
# telescope itself. In `Comrade` we use the RIME formalism [^1] to represent these corruptions,
# namely our measured coherency matrices $V_{ij}$ are given by
#  ```math
#     V_{ij} = J_iC_{ij}J_j^\dagger
#  ```
# where $J$ is known as a *Jones matrix* and $ij$ denotes the baseline $ij$ with sites $i$ and $j$.
# `Comrade` is highly flexible with how the Jones matrices are formed by we provide a number of
# convienience fuctions that parameterize standard Jones matrices. These include:
#   - [`jonesG`](@ref) which builds the set of complex gain Jones matrices
#        ```math
#          G = \begin{pmatrix}
#                  g_a   &0\\
#                  0     &g_b\\
#              \end{pmatrix}
#        ```
#   - [`jonesD`](@ref) which builds the set of complex d-terms Jones matrices
#        ```math
#          G = \begin{pmatrix}
#                  1   & d_a\\
#                  d_b     &1\\
#              \end{pmatrix}
#        ```
#   - [`jonesT`](@ref) is the basis transform matrix. This transform is special and
#      combines two thing. The first is the transformation from
#      some reference basis, to the observed visibility basis (this allows for mixed basis measurements).
#      The second is the feed rotation that transforms from some reference axis to the axis of the
#      telescope, as the source moves in the sky. The feed rotation matrix `F` is given by
#        ```math
#          F = \begin{pmatrix}
#                  e^{-i\varphi}   & 0\\
#                  0     & e^{i\varphi}\\
#              \end{pmatrix}
#        ```
#        where $\varphi$ is the determinstic rotation angle that can be directly obtained from the data.
#
#
#  In the rest of the tutorial we are going to solve for all of these instrument model terms on
#  in addition to our image structure to reconstruct a polarized image of a synthetic dataset.


# ## Load the Data
using Pkg; Pkg.activate(@__DIR__)

using Comrade


# To download the data visit https://doi.org/10.25739/g85n-f134
# To load the eht-imaging obsdata object we do:
obs = load_ehtim_uvfits(joinpath(@__DIR__, "PolarizedExamples/polarized_gaussian_nogains_withdterms_withfr.uvfits"),
                        joinpath(@__DIR__, "PolarizedExamples/array.txt"))
# Notive that unlike other non-polarized tutorials we need to include a second argument.
# This is the array file of the observation and is required to determine the feed rotations
# of the array.

# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      coherent.
obs = scan_average(obs)
#-
# Now we extract our coherency matrices.
dvis = extract_coherency(obs)

# ##Building the Model/Posterior

# Now we must build our intensity/visibility model and the instrument or RIME model
# mentioned above.
# The model below is broken into two parts:
#    1. The image model: In Comrade all polarized image models are written in terms of the Stokes parameters.
#       The reason for this choice is that this usually what physical models consider and it is
#       the easiest to reason about.
#       For our model in this tutorial we will be using a polarized image model based on Pesce (2021)[^2].
#       This model parameterizes the polarized image in terms of the [`Poincare sphere`](https://en.wikipedia.org/wiki/Unpolarized_light#Poincar%C3%A9_sphere),
#       and allows us to easily incorporate physical restrictions such as $I^2 ≥ Q^2 + U^2 + V^2$.
#    2. The instrument or Jones model. This is our model that describes the impact of instrumental and atmospheric effects.
#       We will be using the `J = GDF` decomposition we described above. However, to parameterize the
#       R/L complex gains we will be using a gain product and ratio decomposition. The reason for this decomposition
#       is that in realistic measurements the gain ratios and products have different temporal characteristics.
#       Namely, many of the EHT observations tend to demonstrate constant R/L gain ratios across an
#       nights observations, while the gain products tend to vary from scan to scan. Using this apriori
#       knowledge we can build this into our model and reduce the total number of parameters we need to model.


function model(θ, metadata)
    (;c, f, p, angparams, dRx, dRy, dLx, dLy, lgp, gpp, lgr, gpr) = θ
    # (;c, f, p, angparams, dRx, dRy, dLx, dLy) = θ
    # (;c, f, p, angparams) = θ
    (; grid, cache, tcache, scancache, trackcache) = metadata
    #Construct the image model
    #produce Stokes images from parameters
    imgI = f*c
    #Converts from poincare sphere parameterization of polzarization to Stokes Parameters
    pimg = PoincareSphere2Map(imgI, p, angparams, grid)
    m = ContinuousImage(pimg, cache)

    #Now construct the basis transformation cache
    jT = jonesT(tcache)

    # #Gain product parameters
    gP = exp.(lgp/2 .+ 1im.*gpp/2)
    Gp = jonesG(gP, gP, scancache)
     # Gain ratio
    gR = exp.(lgr/2 .+ 1im.*gpr/2)
    Gr = jonesG(gR, inv.(gR), trackcache)
    # # D-terms
    D = jonesD(complex.(dRx, dRy), complex.(dLx, dLy), trackcache)
    #sandwich all the jones matrices together
    J = Gp*Gr*D*jT
    # J = D*jT
    # J = jT
    #form the complete Jones or RIME model. We use tcache here
    #to set the reference basis of the model.
    return JonesModel(J, m, tcache)
end

# Now can define our model metadata to complete the specification.
# First we specify our image grid and cache model.
fovx = μas2rad(50.0)
fovy = μas2rad(50.0)
nx = 5
ny = floor(Int, fovy/fovx*nx)
grid = imagepixels(fovx, fovy, nx, ny) # image grid
buffer = IntensityMap(zeros(nx, ny), grid) # buffer to store temporary image
pulse = BSplinePulse{3}() # pulse we will be using
cache = create_cache(DFTAlg(dvis), buffer, pulse) # cache to define the NFFT transform

# To define the instrument models, T, G, D we need to define a cache (see [`JonesCache`](@ref)) that maps from a flat
# vector of gain/dterms to the specific sites for each baseline. First we will define our
# deterministic transform cache. Note that this dataset has need been pre-corrected for
# feed rotation so we need to add those into the `tcache`.
tcache = TransformCache(dvis; add_fr=true, ehtim_fr_convention=false)
#-
# Next we define our cache that maps quantities e.g., gain products, that change from scan-to-scan.
scancache = JonesCache(dvis, ScanSeg())
#-
# Finally, we define our cache that maps quantities, e.g., gain ratios and d-terms, that are constant
# across a observation night, and we collect everything together.
trackcache = JonesCache(dvis, TrackSeg())
metadata = (;cache, grid, tcache, scancache, trackcache)

# Moving onto our prior we first focus on the instrument model priors.
# Each station gain requires its own prior on both the amplitudes and phases.
# For the amplitudes we assume that the gains are aprior well calibrated around unit gains (or 0 log gain amplitudes)
# which corresponds to no instrument corruption. The gain dispersion is then set to 10% for
# all stations except LMT representing that from scan-to-scan we expect 10% deviations. For LMT
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
# For the phases we assume that the gains are effectively scrambled by the atmosphere.
# Since the gain phases are periodic we also use a von Mises priors for all stations with
# essentially a flat distribution.
using VLBIImagePriors
distphase = (AA = DiagonalVonMises(0.0, inv(1e-6)),
             AP = DiagonalVonMises(0.0, inv(π^2)),
             LM = DiagonalVonMises(0.0, inv(π^2)),
             AZ = DiagonalVonMises(0.0, inv(π^2)),
             JC = DiagonalVonMises(0.0, inv(π^2)),
             PV = DiagonalVonMises(0.0, inv(π^2)),
             SM = DiagonalVonMises(0.0, inv(π^2)),
           )
#-
# However, we can now also use a little additional information about the phase offsets
# where in most cases they are much better behaved than the products
distphase_ratio = (AA = DiagonalVonMises(0.0, inv(1e-6)),
             AP = DiagonalVonMises(0.0, inv(0.1^2)),
             LM = DiagonalVonMises(0.0, inv(0.1^2)),
             AZ = DiagonalVonMises(0.0, inv(0.1^2)),
             JC = DiagonalVonMises(0.0, inv(0.1^2)),
             PV = DiagonalVonMises(0.0, inv(0.1^2)),
             SM = DiagonalVonMises(0.0, inv(0.1^2)),
           )


# Moving onto the d-terms, here we directly parameterize the real and complex components
# of the d-terms since they are expect to be complex numbers near the origin. To help enforce
# this smallness we put a weakly informative Normal prior on the d-terms.
distD = ( AA = Normal(0.0, 0.1),
          AP = Normal(0.0, 0.1),
          LM = Normal(0.0, 0.1),
          AZ = Normal(0.0, 0.1),
          JC = Normal(0.0, 0.1),
          PV = Normal(0.0, 0.1),
          SM = Normal(0.0, 0.1),
        )



# We can now form our model parameter priors. Like our other imaging examples we use a
# Dirichlet prior for our stokes I image pixels. For the total polarization fraction we assume
# an uncorrelated uniform prior for each pixel. Finally, to specify the orientation of the
# polarization on the Poincare sphere, we use a uniform spherical distribution.
# For all the calibration parameters we use a helper function `CalPrior` which
# builds the prior given the named tuple of station priors and a `JonesCache` that specifies the
# segmentation scheme. For the gain products we use the `scancache` while for every other
# quantity we use the `trackcache`.
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
          gpp = CalPrior(distphase, scancache),
          lgr = CalPrior(distamp, trackcache),
          gpr = CalPrior(distphase_ratio,trackcache),
          )


# Putting it all together we form our likelihood and posterior objects for optimization and
# sampling.
lklhd = RadioLikelihood(model, metadata, dvis)
post = Posterior(lklhd, prior)

# ## Reconstructing the Image and Instrument Effects

# To sample from this posterior it is convienent to first move from our constrained paramter space
# to a unconstrained one (i.e., the support of the transformed posterior is (-∞, ∞)). This is
# done using the `asflat` function.
tpost = asflat(post)

# We can now also find the dimension of our posterior, or the number of parameters we are going to sample.
# !!! warning
#     This can often be different from what you would expect. This is especially true when using
#     angular variables where to make sampling easier we often artifically increase the dimension
#     of the parameter space.


ndim = dimension(tpost)


# Now we optimize. Unlike other imaging examples here we move straight to gradient optimizers
# due to the higher dimension of the space. Additionally, we will run the optimizer multiple times
# due to the overall difficult in optimization for polarized modeling.
using ComradeOptimization
using OptimizationOptimJL
using Zygote
f = OptimizationFunction(tpost, Optimization.AutoZygote())
ℓ = logdensityof(tpost)
prob = OptimizationProblem(f, prior_sample(tpost), nothing)
sol = solve(prob, LBFGS(), maxiters=15_000, callback=((x,p)->(@info ℓ(x);false)), g_tol=1e-1)

# !!! warning
#     Fitting polarized images is generally much harder than stokes I imaging. This means that
#     optimizing properly can take some effort, and starting from a reasonable staring location
#     is often required.

# Before we analyze our solution we first need to transform back to parameter space.
xopt = transform(tpost, sol)

# Now let's evaluate our fits by plotting the residuals
using Plots
residual(model(xopt, metadata), dvis)

# These look reasonable, although maybe there is some minor overfitting. This could probably be
# improved in a few ways, but that is beyond the goal of this quick tutorial.
# Let's quickly compare our results to the ground truth values which we know in this example.
# First we will load the polarized truth
using AxisKeys
imgtrue = Comrade.load(joinpath(@__DIR__, "PolarizedExamples/polarized_gaussian.fits"), StokesIntensityMap)
imgtruesub = imgtrue(Interval(-fovx/2, fovx/2), Interval(-fovy/2, fovy/2))
img = intensitymap!(copy(imgtruesub), model(xopt, metadata))


plot(img, title="Reconstructed Image", xlims=(-25.0,25.0), ylims=(-25.0,25.0))
#-
plot(imgtruesub, title="True Image", xlims=(-25.0,25.0), ylims=(-25.0,25.0))


# Now because we also fit the instrument model we can also inspect their parameters.
# To do this `Comrade` provides a `caltable` function that converts the flattened gain parameters
# to a tabular format based on the time and its segmentation.
dR = caltable(trackcache, complex.(xopt.dRx, xopt.dRy))

# We can compare this to the ground truth d-terms
# ───────────┬─────────────────────────────────────────────────────────────────────────────────────────────────────
#       time │            AA              AP             AZ           JC           LM            PV             SM
# ───────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────
#  0.0+0.0im │ 0.01-0.02im    -0.08+0.07im     0.09-0.10im  -0.04+0.05im   0.03-0.02im   -0.01+0.02im  0.08-0.07im
# ───────────┴──────────────────────────────────────────────────────────────────────────────────────────────────────

# And same for the left-handed dterms
#
dL = caltable(trackcache, complex.(xopt.dLx, xopt.dLy))
#
# ───────────┬──────────────────────────────────────────────────────────────────────────────────────────────────────────
#       time │            AA              AP             AZ             JC             LM            PV             SM
# ───────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────
#  0.0+0.0im │ 0.03-0.04im    -0.06+0.05im     0.09-0.08im    -0.06+0.07im     0.01-0.00im  -0.03+0.04im  0.06-0.05im
# ───────────┴──────────────────────────────────────────────────────────────────────────────────────────────────────────

# Looking at the gain phase ratio
gphase_ratio = caltable(trackcache, xopt.gpr)
# we see that they are all very small. Which should be the case since this data doesn't have gain corruptions!
# Similarly our gain ratio amplitudes are also very close to unity as expected.
gamp_ratio   = caltable(trackcache, exp.(xopt.lgr))
#-
# Plotting the gain phases we see some offsets from zero. This is because the prior on the gain product
# phases is very broad and so we don't have the ability to phase center the image. For realistic data
# this is always the case since the atmosphere effecticely scrambles the phases.
gphase_prod = caltable(scancache, xopt.gpp)
plot(gphase_prod, layout=(3,3), size=(650,500))
#-
# Finally the product gain amplitudes are all very close to unity as well, as expected since gain corruptions
# have not been added to the data.
gamp_prod = caltable(scancache, exp.(xopt.lgp))
plot(gamp_prod, layout=(3,3), size=(650,500))
#-
# At this point you should run the sampler to recover an uncertainty estimate,
# which is identical to every other imaging example
# (see e.g., [ Stokes I simultaneous Image and Instrument Modeling](@ref). However,
# due to the time it takes to sample we will skip that for this tutorial. Note that on computer environment listed
# below 20_000 MCMC steps takes 4 hours.

# # now we sample using hmc
# using ComradeAHMC
metric = DenseEuclideanMetric(ndim)
chain, stats = sample(post, AHMC(;metric, autodiff=Comrade.AD.ZygoteBackend()), 52_000; nadapts=42_000, init_params=xopt)

# ## References
# [^1]: Hamaker J.P, Bregman J.D., Sault R.J. (1996) [https://articles.adsabs.harvard.edu/pdf/1996A%26AS..117..137H]
# [^2]: Pesce D. (2021) [https://ui.adsabs.harvard.edu/abs/2021AJ....161..178P/abstract]


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
