# # Polarized Image and Instrumental Modeling

# In this tutorial, we will analyze a simulated simple polarized dataset to demonstrate
# Comrade's polarized imaging capabilities.

# ## Introduction to Polarized Imaging
# The EHT is a polarized interferometer. However, like all VLBI interferometers, it does not
# directly measure the Stokes parameters (I, Q, U, V). Instead, it measures components
# related to the electric field at the telescope along two *directions* using feeds.
# There are two types of feeds at telescopes: circular, which measure $R/L$ components of the
# electric field, and linear feeds, which measure $X/Y$ components of the electric field.
# Most sites in the EHT use circular feeds, meaning they measure the right (R) and left
# electric field (L) at each telescope. Although note that ALMA actually uses linear feeds.
# Currently Comrade has the ability to fit natively mixed polarization data however, the
# publically released EHT data has been converted to circular polarization.
# For a VLBI array whose feeds are purely circluar the **coherency matrices** are given by,
#
# ```math
#  C_{ij} = \begin{pmatrix}
#        RR^* &  RL^*\\
#        LR^* &  LL^*
#      \end{pmatrix}.
# ```
#
# These coherency matrices are the fundamental object in interferometry and what
# the telescope observes. For a perfect interferometer, these circular coherency matrices
# are related to the usual Fourier transform of the stokes parameters by
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
#   - [`JonesG`](@ref) which builds the set of complex gain Jones matrices
# ```math
#   G = \begin{pmatrix}
#           g_a   &0\\
#           0     &g_b\\
#       \end{pmatrix}
# ```
#   - [`JonesD`](@ref) which builds the set of complex d-terms Jones matrices
# ```math
#   D = \begin{pmatrix}
#           1   & d_a\\
#           d_b     &1\\
#       \end{pmatrix}
# ```
#   - [`JonesR`](@ref) is the basis transform matrix $T$. This transformation is special and
#      combines two things using the decomposition $T=FB$. The first, $B$, is the transformation from
#      some reference basis to the observed coherency basis (this allows for mixed basis measurements).
#      The second is the feed rotation, $F$, that transforms from some reference axis to the axis of the
#      telescope as the source moves in the sky. The feed rotation matrix `F` for circular feeds
#      in terms of the per station feed rotation angle $\varphi$ is
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
Pkg.develop(; path=joinpath(__DIR, "..", "..", ".."), io=pkg_io) #hide
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
rng = StableRNG(42)


# Now we will load some synthetic polarized data.
obs = Pyehtim.load_uvfits_and_array(
        joinpath(__DIR, "..", "..", "Data", "polarized_gaussian_all_corruptions.uvfits"),
        joinpath(__DIR, "..", "..", "Data", "array.txt"), polrep="circ")


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
#       In this tutorial, we will use a polarized image model based on Pesce (2021)[^2], and
#       parameterizes each pixel in terms of the [`Poincare sphere`](https://en.wikipedia.org/wiki/Unpolarized_light#Poincar%C3%A9_sphere).
#       This parameterization ensures that we have all the physical properties of stokes parameters.
#       Note that we also have a parameterization in terms of hyperbolic trig functions `VLBISkyModels.PolExp2Map`
#    2. **The instrument model**. The instrument model specifies the model that describes the impact of instrumental and atmospheric effects.
#       We will be using the $J = GDR$ decomposition we described above. However, to parameterize the
#       R/L complex gains, we will be using a gain product and ratio decomposition. The reason for this decomposition
#       is that in realistic measurements, the gain ratios and products have different temporal characteristics.
#       Namely, many of the EHT observations tend to demonstrate constant R/L gain ratios across an
#       nights observations, compared to the gain products, which vary every scan. Additionally,
#       the gain ratios tend to be smaller (i.e., closer to unity) than the gain products.
#       Using this apriori knowledge, we can build this into our model and reduce
#       the total number of parameters we need to model.


# First we specify our sky model. As always `Comrade` requires this to be a two argument
# function where the first argument is typically a NamedTuple of parameters we will fit
# and the second are additional metadata required to build the model.
using StatsFuns: logistic
function sky(θ, metadata)
    (;c, σ, p, p0, pσ, angparams) = θ
    (;ftot, grid) = metadata
    ## Build the stokes I model
    rast = to_simplex(CenteredLR(), c.params*σ)
    rast .= ftot.*rast
    ## The total polarization fraction is modeled in logit space so we transform it back
    pim = logistic.(p0 .+ pσ.*p.params)
    ## Build our IntensityMap
    pmap = PoincareSphere2Map(rast, pim, angparams, grid)
    ## Construct the actual image model which uses a third order B-spline pulse
    m = ContinuousImage(pmap, BSplinePulse{3}())
    ## Finally find the image centroid and shift it to be at the center
    x0, y0 = centroid(pmap)
    ms = shifted(m, -x0, -y0)
    return ms
end






# Now, we define the model metadata required to build the model.
# We specify our image grid and cache model needed to define the polarimetric
# image model. Our image will be a 10x10 raster with a 60μas FOV.
using Distributions, DistributionsAD
using VLBIImagePriors
fovx = μas2rad(60.0)
fovy = μas2rad(60.0)
nx = ny = 10
grid = imagepixels(fovx, fovy, nx, ny)

# For the image metadata we specify the grid and the total flux of the image, which is 1.0.
# Note that we specify the total flux out front since it is degenerate with an overall shift
# in the gain amplitudes.
skymeta = (;ftot=1.0, grid)


# For our image prior we will use a simpler prior than
#   - We use again use a GMRF prior. For more information see the [Imaging a Black Hole using only Closure Quantities](@ref) tutorial.
#   - For the total polarization fraction, `p`, we assume an uncorrelated uniform prior `ImageUniform` for each pixel.
#   - To specify the orientation of the polarization, `angparams`, on the Poincare sphere,
#     we use a uniform spherical distribution, `ImageSphericalUniform`.
#-
rat = beamsize(dvis)/step(grid.X)
cmarkov = ConditionalMarkov(GMRF, grid; order=1)
dρ = truncated(InverseGamma(1.0, -log(0.1)*rat); lower=2.0, upper=max(nx, ny))
cprior = HierarchicalPrior(cmarkov, dρ)


# For all the calibration parameters, we use a helper function `CalPrior` which builds the
# prior given the named tuple of station priors and a `JonesCache`
# that specifies the segmentation scheme. For the gain products, we use the `scancache`, while
# for every other quantity, we use the `trackcache`.
fwhmfac = 2.0*sqrt(2.0*log(2.0))
skyprior = (
    c = cprior,
    σ  = truncated(Normal(0.0, 0.1); lower=0.0),
    p  = cprior,
    p0 = Normal(-2.0, 2.0),
    pσ =  truncated(Normal(0.0, 1.0); lower=0.01),
    angparams = ImageSphericalUniform(nx, ny),
    )

skym = SkyModel(sky, skyprior, grid; metadata=skymeta)


# Now we build the instrument model. Due to the complexity of VLBI the instrument model is critical
# to the success of imaging and getting reliable results. For this example we will use the standard
# instrument model used in polarized EHT analyses expressed in the RIME formalism. Our Jones
# decomposition will be given by `GDR`, where `G` are the complex gains, `D` are the d-terms, and `R`
# is what we call the *ideal instrument response*, which is how an ideal interferometer using the
# feed basis we observe relative to some reference basis.
#
# Given the possible flexibility in different parameterizations of the individual Jones matrices
# each Jones matrix requires the user to specify a function that converts from parameters
# to specific parameterization f the jones matrices.

# For the complex gain matrix, we used the `JonesG` jones matrix. The first argument is now
# a function that converts from the parameters to the complex gain matrix. In this case, we
# will use a amplitude and phase decomposition of the complex gain matrix. Note that since
# the gain matrix is a diagonal 2x2 matrix the function must return a 2-element tuple.
# The first element of the tuple is the gain for the first polarization feed (R) and the
# second is the gain for the second polarization feed (L).
function fgain(x)
    gR = exp(x.lgR + 1im*x.gpR)
    gL = gR*exp(x.lgrat + 1im*x.gprat)
    return gR, gL
end
G = JonesG(fgain)
# Note that we are using the Julia `do` syntax here to define an anonymous function. This
# could've also been written as
# ```julia
# fgain(x) = (exp(x.lgR + 1im*x.gpR), exp(x.lgR + x.lgrat + 1im*(x.gpR + x.gprat)))
# G = JonesG(fgain)
# ```


# Similarly we provide a `JonesD` function for the leakage terms. Since we assume that we
# are in the small leakage limit, we will use the decomposition
# 1 d1
# d2 1
# Therefore, there are 2 free parameters for the JonesD our parameterization function
# must return a 2-element tuple. For d-terms we will use a re-im parameterization.
function fdterms(x)
    dR = complex(x.dRx, x.dRy)
    dL = complex(x.dLx, x.dLy)
    return dR, dL
end

D = JonesD(fdterms)

# Finally we define our response Jones matrix. This matrix is a basis transform matrix
# plus the feed rotation angle for each station. These are typically set by the telescope
# so there are no free parameters, so no parameterization is necessary.
R = JonesR(;add_fr=true)

# Finally, we build our total Jones matrix by using the `JonesSandwich` function. The
# first argument is a function that specifies how to combine each Jones matrix. In this case,
# we are completely standard so we just need to multiply the different jones matrices.
# Note that if no function is provided, the default is to multiply the Jones matrices,
# so we could've removed the * argument in this case.
J = JonesSandwich(splat(*), G, D, R)


intprior = (
    lgR  = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1))),
    gpR  = ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π  ^2))); refant=SEFDReference(0.0), phase=false),
    lgrat= ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1))),
    gprat= ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1)); refant = SingleReference(:AA, 0.0)),
    dRx  = ArrayPrior(IIDSitePrior(TrackSeg(), Normal(0.0, 0.2))),
    dRy  = ArrayPrior(IIDSitePrior(TrackSeg(), Normal(0.0, 0.2))),
    dLx  = ArrayPrior(IIDSitePrior(TrackSeg(), Normal(0.0, 0.2))),
    dLy  = ArrayPrior(IIDSitePrior(TrackSeg(), Normal(0.0, 0.2))),
)

intmodel = InstrumentModel(J, intprior)


# Putting it all together, we form our likelihood and posterior objects for optimization and
# sampling.
post = VLBIPosterior(skym, intmodel, dvis)

# ## Reconstructing the Image and Instrument Effects

# To sample from this posterior, it is convenient to move from our constrained parameter space
# to an unconstrained one (i.e., the support of the transformed posterior is (-∞, ∞)). This transformation is
# done using the `asflat` function.
tpost = asflat(post)

# We can also query the dimension of our posterior or the number of parameters we will sample.
# !!! warning
#     This can often be different from what you would expect. This difference is especially true when using
#     angular variables, where we often artificially increase the dimension
#     of the parameter space to make sampling easier.
#-

ndim = dimension(tpost)

using Enzyme
Enzyme.API.runtimeActivity!(true)
x = prior_sample(rng, tpost)
dx = zero(x)
autodiff(Enzyme.Reverse, logdensityof, Active, Const(tpost), Duplicated(x, dx))

using BenchmarkTools
@benchmark autodiff($Enzyme.Reverse, $logdensityof,$Active, $(Const(tpost)), Duplicated($x, fill!($dx, 0)))


# Now we optimize. Unlike other imaging examples, we move straight to gradient optimizers
# due to the higher dimension of the space. In addition the only AD package that can currently
# work with the polarized Comrade posterior is Enzyme.
using Optimization
using OptimizationOptimisers
using Enzyme
xopt, sol = comrade_opt(post, Optimisers.Adam(), AutoEnzyme(;mode=Enzyme.Reverse); initial_params=prior_sample(rng, post), maxiters=20_000)


# !!! warning
#     Fitting polarized images is generally much harder than Stokes I imaging. This difficulty means that
#     optimization can take a long time, and starting from a good starting location
#     is often required.
#-

# Now let's evaluate our fits by plotting the residuals
using Plots
residual(post, xopt)

# These look reasonable, although there may be some minor overfitting.
# Let's compare our results to the ground truth values we know in this example.
# First, we will load the polarized truth
imgtrue = load_fits(joinpath(__DIR, "..", "..", "Data", "polarized_gaussian.fits"), IntensityMap{StokesParams})
# Select a reasonable zoom in of the image.
imgtruesub = regrid(imgtrue, imagepixels(fovx, fovy, nx*4, ny*4))
img = intensitymap(Comrade.skymodel(post, xopt), axisdims(imgtruesub))

#Plotting the results gives
import CairoMakie as CM
using DisplayAs
CM.activate!(type = "png", px_per_unit=3) #hide
fig = CM.Figure(;size=(450, 350));
axs = [CM.Axis(fig[1, i], xreversed=true, aspect=1) for i in 1:2]
polimage!(axs[1], imgtruesub,
                   nvec = 8,
                   length_norm=1/2, plot_total=true, pcolormap=:RdBu,
                   pcolorrange=(-0.25, 0.25),); axs[1].title="True"
polimage!(axs[2], img,
                   nvec = 8,
                   length_norm=1/2, plot_total=true, pcolormap=:RdBu,
                   pcolorrange=(-0.25, 0.25),);axs[2].title="Recon."
CM.Colorbar(fig[2,:], colormap=:RdBu, vertical=false, colorrange=(-0.25, 0.25), label="Signed Polarization Fraction sign(V)*|p|", flipaxis=false)
CM.colgap!(fig.layout, 3)
CM.rowgap!(fig.layout, 3)
CM.hidedecorations!.(fig.content[1:2])
fig |> DisplayAs.PNG |> DisplayAs.Text
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
dR = caltable(complex.(xopt.instrument.dRx, xopt.instrument.dRy))

# We can compare this to the ground truth d-terms
#
# |      time |            AA   |           AP   |         AZ   |        JC    |        LM    |       PV     |      SM |
# |-----------|-----------------|----------------|--------------|--------------|--------------|--------------|---------|
# | 0.0       | 0.01-0.02im      | -0.08+0.07im  |  0.09-0.10im | -0.04+0.05im |  0.03-0.02im | -0.01+0.02im | 0.08-0.07im |
#

# And same for the left-handed dterms
#
dL = caltable(complex.(xopt.instrument.dLx, xopt.instrument.dLy))
#
# |      time |            AA   |           AP   |         AZ   |        JC    |        LM    |       PV     |      SM |
# |-----------|-----------------|----------------|--------------|--------------|--------------|--------------|---------|
# | 0.0       | 0.03-0.04im      | -0.06+0.05im  |  0.09-0.08im | -0.06+0.07im |  0.01-0.00im | -0.03+0.04im | 0.06-0.05im |
#


# Looking at the gain phase ratio
gphase_ratio = caltable(xopt.instrument.gprat)
#-
# we see that they are all very small. Which should be the case since this data doesn't have gain corruptions!
# Similarly our gain ratio amplitudes are also very close to unity as expected.
gamp_ratio   = caltable(exp.(xopt.instrument.lgrat))
#-
# Plotting the gain phases, we see some offsets from zero. This is because the prior on the gain product
# phases is very broad, so we can't phase center the image. For realistic data
# this is always the case since the atmosphere effectively scrambles the phases.
gphaseR = caltable(xopt.instrument.gpR)
p = Plots.plot(gphaseR, layout=(3,3), size=(650,500));
Plots.plot!(p, gphase_ratio, layout=(3,3), size=(650,500));
p |> DisplayAs.PNG |> DisplayAs.Text
#-
# Finally, the product gain amplitudes are all very close to unity as well, as expected since gain corruptions
# have not been added to the data.
gampr = caltable(exp.(xopt.instrument.lgR))
p = Plots.plot(gampr, layout=(3,3), size=(650,500))
Plots.plot!(p, gamp_ratio, layout=(3,3), size=(650,500))
p |> DisplayAs.PNG |> DisplayAs.Text
#-

# [^1]: Hamaker J.P, Bregman J.D., Sault R.J. (1996) [https://articles.adsabs.harvard.edu/pdf/1996A%26AS..117..137H]
# [^2]: Pesce D. (2021) [https://ui.adsabs.harvard.edu/abs/2021AJ....161..178P/abstract]
#-
