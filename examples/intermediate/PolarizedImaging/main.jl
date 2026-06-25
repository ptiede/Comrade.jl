import Pkg #hide
__DIR = @__DIR__ #hide
pkg_io = open(joinpath(__DIR, "pkg.log"), "w") #hide
Pkg.activate(__DIR; io = pkg_io) #hide
Pkg.develop(; path = joinpath(__DIR, "..", "..", ".."), io = pkg_io) #hide
Pkg.instantiate(; io = pkg_io) #hide
Pkg.precompile(; io = pkg_io) #hide
close(pkg_io) #hide

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
#   - [`JonesR`](@ref) is the ideal telescope response $R$. This transformation is special and
#      combines two things using the decomposition $R=FB$. The first, $B$, is the transformation from
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
#  In the rest of the tutorial, we are going to solve for all of these instrument model terms in
#  while re-creating the polarized image from the first [`EHT results on M87`](https://iopscience.iop.org/article/10.3847/2041-8213/abe71d).


# ## Load the Data
# To get started we will load Comrade
using Comrade

# ## Load the Data
using VLBIFiles

# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(48)


# Now we will load some synthetic polarized data.
fname = Base.download(
    "https://github.com/ptiede/ComradeTestData/releases/download/Data/eht_2017_data.uvfits",
    joinpath(__DIR, "m87polarized.uvfits")
)
uvd = VLBIFiles.load(VLBIFiles.UVData, fname)


# Unlike non-polarized tutorials, we also need an **array file** describing the antenna
# feed rotation parameters. We pass it via the `arrayfile` keyword on the data product.
dvis = extract_table(
    uvd, Coherencies(;
        time_average = VLBI.GapBasedScans(),
    )
)

# Often the EHT does not take into acount feed offsets. To load those we read in the standard
# ehtim style array.txt files.
Comrade.reset_mounts!(dvis, Comrade.load_array_txt(joinpath(__DIR, "../../Data/array.txt")))
# Inflate noise by 1% and drop short (uvdist < 0.1 Gλ) baselines.
add_fractional_noise!(dvis, 0.01)
dvis = flag(d -> hypot(d.baseline.U, d.baseline.V) < 0.1e9, dvis)

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


# First we specify our sky model with the `@sky` macro. Each `name ~ dist` line contributes
# an entry to the prior; everything else is the model body. The metadata fields (`mimg`,
# `ftot`, `cprior`) are passed as keyword arguments and are in scope inside both the prior
# expressions and the body.
using Distributions
using VLBIImagePriors
@sky function sky(grid; mimg, ftot, cprior)
    σs ~ ntuple(Returns(VLBITruncated(VLBIGaussian(0.0, 0.5); lower = 0.0)), 4)
    as ~ ntuple(Returns(cprior), 4)
    ## Build the stokes I model
    δs = ntuple(Val(4)) do i
        σs[i] .* as[i].params
    end

    ## Convert hyperbolic polarization parameterization to Stokes parameters.
    pmap = VLBISkyModels.PolExp2Map!(δs..., axisdims(mimg))

    ## We now add a mean image. Namely, we assume that `pmap` are multiplicative fluctuations
    ## about some mean image `mimg`. We also compute the total flux of the Stokes I image
    ## for normalization purposes below.
    bpmap = baseimage(pmap)
    bmimg = baseimage(mimg)
    ft = zero(ftot)
    bpmapI = stokes(bpmap, :I)
    bpmapQ = stokes(bpmap, :Q)
    bpmapU = stokes(bpmap, :U)
    bpmapV = stokes(bpmap, :V)

    @inbounds for i in eachindex(bpmap, bmimg)
        bpmapI[i] = bmimg[i] * bpmapI[i]
        bpmapQ[i] = bmimg[i] * bpmapQ[i]
        bpmapU[i] = bmimg[i] * bpmapU[i]
        bpmapV[i] = bmimg[i] * bpmapV[i]
        ft += bpmapI[i]
    end

    ## Renormlize to get the correct total flux.
    @inbounds for i in eachindex(bpmap)
        bpmap[i] *= ftot / ft
    end

    m = ContinuousImage(pmap, BSplinePulse{3}())
    x, y = centroid(pmap)
    return shifted(m, -x, -y)
end


# Now, we define the metadata needed to build the model.
# We specify our image grid and cache model needed to define the polarimetric
# image model. Our image will be a 10x10 raster with a 60μas FOV.
fovx = μas2rad(200.0)
fovy = μas2rad(200.0)
nx = ny = 48
grid = imagepixels(fovx, fovy, nx, ny)

fwhmfac = 2 * sqrt(2 * log(2))
mpr = modify(Gaussian(), Stretch(μas2rad(50.0) ./ fwhmfac))
mimg_raw = intensitymap(mpr, grid)
mimg = mimg_raw ./ Comrade.flux(mimg_raw)


# We use again use a GMRF prior similar to the [Imaging a Black Hole using only Closure Quantities](@ref) tutorial
# for the log-ratio transformed image. We use the same correlated image prior for the inverse-logit transformed
# total polarization. The mean total polarization fraction `p0` is centered at -2.0 with a standard deviation of 2.0
# which logit transformed puts most of the prior mass < 0.8 fractional polarization. The standard deviation of the
# total polarization fraction `pσ` again uses a Half-normal process. The angular parameters of the polarizaton are
# given by a uniform prior on the sphere.
cprior = corr_image_prior(grid, dvis; order = 2)

# The total flux is fixed to 0.6 since it is degenerate with an overall shift in the gain
# amplitudes.
skym = sky(grid; mimg, ftot = 0.6, cprior)


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

# We bundle the whole instrument model — the Jones matrices and their priors — in a single
# block with the `@instrument` macro. Each `name ~ ArrayPrior(...)` line adds an entry to
# the instrument prior; the rest of the body builds and returns the full Jones matrix. The
# `param_map`s reference the sampled-parameter names directly (no dummy `x`), and the macro
# lifts each one — and the combination function — into a named function for us (named
# functions, unlike closures, serialize reliably).
#
# The individual Jones matrices are:
#   - `JonesG` — the complex gains, in an amplitude/phase decomposition. The gain matrix is
#     diagonal, so the `param_map` returns a 2-tuple (R feed, L feed).
#   - `JonesD` — the leakage (d-)terms in the small-leakage limit `[1 d1; d2 1]`, using a
#     re-im parameterization, returning a 2-tuple.
#   - `JonesR` — the response matrix (basis transform plus feed rotation). It has no free
#     parameters.
# We combine them with `JonesSandwich` using the standard `J = adjoint(R)*G*D*R`.

# !!! note
#     For arrays with non-zero leakage, feed rotation calibration does not remove the impact
#     of feed rotations on the instrument model: the `R` and `D` matrices do not commute, so
#     feed rotation must be modelled. This is why we keep `add_fr = true` and apply
#     `adjoint(R)` in the composition. This is especially important for interferometers
#     using a mixture of linear and circular feeds.

# For the instrument prior we use IID priors for the gains and d-terms. The `IIDSitePrior`
# segments are `ScanSeg()` (independent per scan), `TrackSeg()` (constant over the track),
# and `IntegSeg()` (changes each integration time). For released EHT data the gains are
# stable over a scan, while the d-terms are stable over the track.
@instrument function instrument()
    lgR ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.2)); LM = IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 1.0)))
    lgrat ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
    gpR ~ ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))); refant = SEFDReference(0.0), phase = true)
    gprat ~ ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(0.1^2))); refant = SingleReference(:AA, 0.0), phase = true)
    dRx ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
    dRy ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
    dLx ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))
    dLy ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.2)))

    ## Complex gains: amplitude/phase decomposition, returning (gR, gL).
    G = JonesG() do
        gR = exp(complex(lgR, gpR))
        gL = gR * exp(complex(lgrat, gprat))
        return gR, gL
    end
    ## Leakage/d-terms (re-im parameterization), returning (dR, dL).
    D = JonesD((complex(dRx, dRy), complex(dLx, dLy)))
    ## The ideal response (basis transform + feed rotation); no free parameters.
    R = JonesR(; add_fr = true)
    ## Combine into the full Jones matrix J = adjoint(R)*G*D*R.
    return JonesSandwich(G, D, R) do g, d, r
        return adjoint(r) * g * d * r
    end
end

# Building the instrument model is now a single call.
intmodel = instrument()

# intmodel = InstrumentModel(JonesR(;add_fr=true))
# Putting it all together, we form our likelihood and posterior objects for optimization and
# sampling, and specifying to use Enzyme.Reverse with runtime activity for AD.
using Enzyme
post = VLBIPosterior(skym, intmodel, dvis);

# ## Reconstructing the Image and Instrument Effects

# To sample from this posterior, it is convenient to move from our constrained parameter space
# to an unconstrained one (i.e., the support of the transformed posterior is (-∞, ∞)). This transformation is
# done using the `asflat` function.
tpost = asflat(post);

# We can also query the dimension of our posterior or the number of parameters we will sample.
# !!! warning
#     This can often be different from what you would expect. This difference is especially true when using
#     angular variables, where we often artificially increase the dimension
#     of the parameter space to make sampling easier.
#-

# Now we optimize. Unlike other imaging examples, we move straight to gradient optimizers
# due to the higher dimension of the space. In addition the only AD package that can currently
# work with the polarized Comrade posterior is Enzyme.
using Optimization, OptimizationLBFGSB
xopt, sol = comrade_opt(
    post, LBFGSB();
    initial_params = prior_sample(rng, post), maxiters = 2_000
)


# Now let's evaluate our fits by plotting the residuals
using CairoMakie
using DisplayAs #hide
res = residuals(post, xopt)
fig = Figure(; size = (800, 600))
plotfields!(fig[1, 1], res[1], uvdist, x -> Comrade.measurement(x)[1, 1] / noise(x)[1, 1], axis_kwargs = (; ylabel = "RR Residual"))
plotfields!(fig[2, 1], res[1], uvdist, x -> Comrade.measurement(x)[2, 1] / noise(x)[2, 1], axis_kwargs = (; ylabel = "LR Residual"))
plotfields!(fig[1, 2], res[1], uvdist, x -> Comrade.measurement(x)[1, 2] / noise(x)[1, 2], axis_kwargs = (; ylabel = "RL Residual"))
plotfields!(fig[2, 2], res[1], uvdist, x -> Comrade.measurement(x)[2, 2] / noise(x)[2, 2], axis_kwargs = (; ylabel = "LL Residual"))
fig |> DisplayAs.PNG |> DisplayAs.Text

# These look reasonable, although there may be some minor overfitting.
# Let's plot the results
gpl = refinespatial(grid, 2)
img = intensitymap(skymodel(post, xopt), gpl)
fig = imageviz(
    img, adjust_length = true, colormap = :cmr_gothic, pcolormap = :rainbow1,
    pcolorrange = (0.0, 0.3), plot_total = false
);
fig |> DisplayAs.PNG |> DisplayAs.Text
#-

# !!! note
#     The image looks a little noisy. This is an artifact of the MAP image. To get a publication quality image
#     we recommend sampling from the posterior and averaging the samples. The results will be essentially
#     identical to the results from [EHTC VII](https://iopscience.iop.org/article/10.3847/2041-8213/abe71d).

# We can also analyze the instrument model. For example, we can look at the gain ratios and products.
# To grab the ratios and products we can use the `caltable` function which will return analyze the gprat array
# and convert it to a uniform table. We can then plot the gain phases and amplitudes.
gphase_ratio = caltable(xopt.instrument.gprat)
gamp_ratio = caltable(exp.(xopt.instrument.lgrat))

#-
# Plotting the phases first, we see large trends in the righ circular polarization phase. This is expected
# due to a lack of image centroid and the absense of absolute phase information in VLBI. However, the gain
# phase difference between the left and right circular polarization is stable and close to zero. This is
# expected since gain ratios are typically stable over the course of an observation and the constant
# offset was removed in the EHT calibration process.
gphaseR = caltable(xopt.instrument.gpR)
fig = plotcaltable(gphaseR, gphase_ratio, labels = ["R Phase", "L/R Phase"]);
fig |> DisplayAs.PNG |> DisplayAs.Text
#-
# Moving to the amplitudes we see largely stable gain amplitudes on the right circular polarization except for LMT which is
# known and due to pointing issues during the 2017 observation. Again the gain ratios are stable and close to unity. Typically
# we expect that apriori calibration should make the gain ratios close to unity.
gampr = caltable(exp.(xopt.instrument.lgR))
fig = plotcaltable(gampr, gamp_ratio, labels = ["R Amp", "L/R Amp"], axis_kwargs = (; limits = (nothing, (0.6, 1.3))));
fig |> DisplayAs.PNG |> DisplayAs.Text
#-


# To sample from the posterior, you can then just use the `sample` function from AdvancedHMC like in the
# other imaging examples. For example
# ```julia
# using AdvancedHMC
# chain = sample(rng, post, NUTS(0.8), 2_000, n_adapts = 1000, progress = true, initial_params = xopt)
# ```


# [^1]: Hamaker J.P, Bregman J.D., Sault R.J. (1996) [https://articles.adsabs.harvard.edu/pdf/1996A%26AS..117..137H]
# [^2]: Pesce D. (2021) [https://ui.adsabs.harvard.edu/abs/2021AJ....161..178P/abstract]
#-
