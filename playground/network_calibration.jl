using Comrade
using Enzyme
using Optimization
using OptimizationOptimisers
using AdvancedHMC
using Distributions, DistributionsAD
using CairoMakie
using Plots
using Pyehtim

function network_calibration(obs::EHTObservationTable{<:Comrade.EHTVisibilityAmplitudeDatum},
    zbl_flux::Real,
    netcal_bl::NTuple{2, Symbol}...;
    gamp_σ = 0.3)

    obsnc = Comrade.prepare_netcal_data(obs, netcal_bl...)
    skym = Comrade.NetworkCalSkyModel(zbl_flux, netcal_bl)

    netcal_prior = (
        AA = IIDSitePrior(IntegSeg(), Normal(0.0, 0.1)),
        AX = IIDSitePrior(IntegSeg(), Normal(0.0, gamp_σ)),
        SW = IIDSitePrior(IntegSeg(), Normal(0.0, gamp_σ)),
        MM = IIDSitePrior(IntegSeg(), Normal(0.0, gamp_σ)),
    )
    intprior = (
        lg = ArrayPrior(IIDSitePrior(IntegSeg(), Normal(0.0, 0.001));
                        netcal_prior...),
    )

    J = SingleStokesGain(x->@inline exp(x.lg))
    intm = InstrumentModel(J, intprior)
    # return obsnc
    post = VLBIPosterior(skym, intm, obsnc)
    return post, obsnc
end
