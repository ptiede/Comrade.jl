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
#      the per sites feed rotation angle $\varphi$ is
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
using StatsFuns
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


function f(x)
    gR = exp(x.lgR + 1im*x.gpR)
    gL = gR*exp(x.lgrat + 1im*x.gprat)
    return gR, gL
end

G = JonesG(f)


D = JonesD() do x
    return complex(x.dRx, x.dRy), complex(x.dLx, x.dLy)
end

R = JonesR()

Jones() do x
    gR = exp(x.lgR + 1im*x.gpR)
    gL = gR*exp(x.lgrat + 1im*x.gprat)
    return JonesG(gR, gL) * D(dr, dl) * R()
end

J = JonesSandwich(G, D) do g, d, r
    return (g + d0)*d*f*d2
end

sitepriors = (
    lgR   = SitePrior(Normal(0.0, 0.1), ScanSeg()); LM = SitePrior(Normal(), ScanSeg()),
    lgrat = SitePrior(Normal(0.0, 0.1), TrackSeg()),
    gpR   = SitePrior(Normal(0.0, 0.1), ScanSeg()); refant = SEFDReference(0.0),
    gprat = SitePrior(Normal(0.0, 0.1), ScanSeg()); refant = SingleReference(:AA, 0.0),
    dRx   = SitePrior(Normal(0.0, 0.1), TrackSeg()),
    dRy   = SitePrior(Normal(0.0, 0.1), TrackSeg()),
    dLx   = SitePrior(Normal(0.0, 0.1), TrackSeg()),
    dLy   = SitePrior(Normal(0.0, 0.1), TrackSeg())
)

@instrument function int(array)
    lgR   ~ SitePrior(Normal(0.0, 0.1), ScanSeg()); LM = SitePrior(Normal(), ScanSeg())
    lgrat ~ SitePrior(Normal(0.0, 0.1), TrackSeg())
    gpR   ~ SitePrior(Normal(0.0, 0.1), ScanSeg()); refant = SEFDReference(0.0)
    gprat ~ SitePrior(Normal(0.0, 0.1), ScanSeg()); refant = SingleReference(:AA, 0.0)
    dRx   ~ SitePrior(Normal(0.0, 0.1), TrackSeg())
    dRy   ~ SitePrior(Normal(0.0, 0.1), TrackSeg())
    dLx   ~ SitePrior(Normal(0.0, 0.1), TrackSeg())
    dLy   ~ SitePrior(Normal(0.0, 0.1), TrackSeg())

    J = Jones() do x
        gR = exp(x.lgR + 1im*x.gpR)
        gL = gR*exp(x.lgrat + 1im*x.gprat)
        return JonesFromG(gR, gL) * D()
    end
    return J
end

intm = InstrumentModel(J, prior, array; refbasis=CirBasis())

using Distributions
using DistributionsAD
using VLBIImagePriors


@skymodel function sky(grid, ftot, beam)
    (;c, σ, p, p0, pσ, angparams) = θ
    c ~ HierarchicalPrior(ConditionalMarkov(GMRF, grid), InverseGamma(1.0, -log(0.1)*beam/rat))
    σ ~ truncated(Normal(0.0, 0.1); lower=0.0)
    p ~ HierarchicalPrior(ConditionalMarkov(GMRF, grid), InverseGamma(1.0, -log(0.1)*beam/rat))
    p0 ~ Normal(-2.0, 1.0)
    pσ ~ truncated(Normal(0.0, 1.0); lower=0.01)
    angparams ~ ImageSphericalUniform(size(grid))

    cp = σ.*c.params
    rast = to_simplex(CenteredLR(), cp)
    pim = logistic.(p0 .+ pσ.*p.params)
    m = PoincareSphere2Map(rast, pim, angparams, grid, BSplinePulse{3}())
    x0, y0 = centroid(stokes(m, :I))
    ms = shifted(m, -x0, -y0)
    return ms
end


# We want to use a correlated image priors to model the polarized image.
# We will use a GMRF prior with a inverse gamma prior on the GMRF correlation length.
grid = imagepixels(fovx, fovy, nx, ny)
skym = skym(grid, 1.0, beam)

# Putting it all together, we form our likelihood and posterior objects for optimization and
# sampling.
post = VLBIPosterior(sky, instrument, obs...)

using ComradeOptimization
using OptimizationOptimisers
using Zygote
opt = optimal_synthetsis(post, AutoZygote(), start, alg)

# Now let's evaluate our fits by plotting the residuals
using Plots
residual(vlbimodel(post, opt), dvis)

trace = sample(post, AHMC(), 1000; n_adapts=500, progress=true)

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


# [^1]: Hamaker J.P, Bregman J.D., Sault R.J. (1996) [https://articles.adsabs.harvard.edu/pdf/1996A%26AS..117..137H]
# [^2]: Pesce D. (2021) [https://ui.adsabs.harvard.edu/abs/2021AJ....161..178P/abstract]
#-
