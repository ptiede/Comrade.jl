# ## Load the Data
# To get started we will load Comrade
using Comrade

# ## Load the Data
using Pkg #hide
Pkg.activate(joinpath(dirname(pathof(Comrade)), "..", "examples")) #hide
using Pyehtim

using JSServe: Page # hide
Page(exportable=true, offline=true) # hide


# For reproducibility we use a stable random number genreator
using StableRNGs
rng = StableRNG(125)


# Now we will load some synthetic polarized data.
obs = Pyehtim.load_uvfits_and_array(joinpath(dirname(pathof(Comrade)), "..", "examples", "PolarizedExamples/polarized_gaussian_all_corruptions.uvfits"),
                        joinpath(dirname(pathof(Comrade)), "..", "examples", "PolarizedExamples/array.txt"), polrep="circ")


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


m1 = Parameterize(MRing{2}()) do x
    return x.α, x.β
end

m2 = Parameterize(ContinuousImage(cache)) do x
    return to_simplex(CenteredLR(), x.c.params)
end

m = m1+m2

setmodel(m, prior_sample(post))

intensitymap()

@sky function polimage(ftot, cache)
    σ    ~ truncated(Normal(0.0, 0.1); lower = 0.0)
    c    ~ HierarchicalPrior(ConditionalMarkov(GMRF, cache), InverseGamma())
    p    ~ HierarchicalPrior(ConditionalMarkov(GMRF, cache), InverseGamma())
    p0   ~ Normal(-2.0, 1.0)
    pσ   ~ truncated(Normal(0.0, 0.1); lower = 0.0)
    s    ~ ImageSphericalUniform(size(K))
    img  = ftot*to_simplex(CenteredLR(), σ.*params(c))
    pimg = logistic.(p0 + pσ.*params(p))
    return PoincareMap(img, pimg, s, cache)
end

skym = polimage(1.0)


G = JonesMatrix(JonesG(prior1, prior2)) do x
    return exp.(x.lgR), exp.(x.lgR .+ x.lgrat)
end

D = JonesMatrix(JonesD(), TrackSeg(), TrackSeg(), TrackSeg()) do x
    return complex(x.dRx, x.dRy), complex(x.dLx, x.dLy)
end

R = JonesMatrix(JonesR(array))

JM = JonesModel(G, D, R) do g, d, r
    return adjoint(r)*g*d*r
end

instrumentmodel = Instrument(JM, array)







@instrument function rime(array)
    R := JonesR(array)
    F = JonesF(array)
    logGaR ~ CalPrior(pgr, array, ScanSeg())
    GaRat  ~ CalPrior(pgr, array, TrackSeg())
    Gp     ~ PhasePrior(phgr, array, ScanSeg(); autoref=SEFDReference(0.0))
    GpRat  ~ PhasePrior(phgr, array, TrackSeg(); autoref=SingleReference(:AA, 0.0))

    dRre  ~ CalPrior(dd, array, TrackSeg())
    dRim  ~ CalPrior(dd, array, TrackSeg())
    dLre  ~ CalPrior(dd, array, TrackSeg())
    dLim  ~ CalPrior(dd, array, TrackSeg())

    D =  JonesD(complex.(dRre, dRim), complex.(dLre, dLim))
    Ga = JonesG(exp.(logGaR), exp.(logGaR .+ GaRat))*JonesG(1, exp.(GaRat))
    Gp = JonesG(exp.(1im*Gp), exp.(1im*(Gp)))*JonesG(1, exp.(1im*(GpRat)))

    J = adjoint.(F).*Ga.*Gp.*D.*R

    return JonesModel(J)
end

instmument = rime(dvis)

post = VLBIPosterior(skym, instrument, dvis; grid)

out = optimize(post, LBFGS(); maxiters=15_000, g_tol=1e-1)

mopt = skymodel(post, out)

sample(post, AHMC(), nsample=10_000, n_adapts=5_000; initial_params=out)

caltabs = calibration_tables(post, samples[1])


ms = skymodel.(Ref(post), chain)

imageviz(intensitymap(ms, Ref(g)) |> mean)
