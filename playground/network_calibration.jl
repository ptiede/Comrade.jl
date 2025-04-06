using Comrade
using Enzyme
using Optimization
using OptimizationOptimisers
using OptimizationOptimJL
using AdvancedHMC
using Distributions, DistributionsAD
using CairoMakie
using Plots
using Pyehtim

function network_calibration(obs::EHTObservationTable{<:Comrade.EHTVisibilityAmplitudeDatum},
    zbl_flux::Real,
    netcal_bl::NTuple{2, Symbol}...;
    gamp_σ = 0.4,
    optimizer = Adam(),
    sample = false)

    obsnc = Comrade.prepare_netcal_data(obs, netcal_bl...)
    skym = Comrade.NetworkCalSkyModel(zbl_flux, netcal_bl)

    intrasites = Set(Iterators.flatten(netcal_bl))

    # intprior = (
    #     lgz = ArrayPrior(IIDSitePrior(IntegSeg(), Normal());
    #                     refant=MultiReference(setdiff(sites(obs), intrasites), 0.0)),
    #     lgμ = ArrayPrior(IIDSitePrior(TrackSeg(), Normal(0.0, 0.2));
    #                     AA = IIDSitePrior(TrackSeg(), Normal(0.0, 0.05)),
    #                     refant=MultiReference(setdiff(sites(obs), intrasites), 0.0)
    #                     ),
    #     lgσ = ArrayPrior(IIDSitePrior(TrackSeg(), Normal(log(0.3), gamp_σ));
    #                     AA = IIDSitePrior(TrackSeg(), Normal(log(0.05), 0.4)),
    #                     refant=MultiReference(setdiff(sites(obs), intrasites), 0.0)
    #                     ),
    # )
    # J = SingleStokesGain(x->@inline exp(x.lgμ + exp(x.lgσ) * x.lgz))

    intprior = (
        lg = ArrayPrior(IIDSitePrior(IntegSeg(), Normal(0.0, 0.3));
                        AA = IIDSitePrior(IntegSeg(), Normal(0.0, 0.05)),
                        refant=MultiReference(setdiff(sites(obs), intrasites), 0.0)),
    )

    J = SingleStokesGain(x->@inline(exp(x.lg)))

    intm = InstrumentModel(J, intprior)
    # return obsnc
    post = VLBIPosterior(skym, intm, obsnc)

    xopt, _ = comrade_opt(post, optimizer, AutoEnzyme(Enzyme.reverse); maxiters=10_000, g_tol=1e-1)
    return xopt, post, obsnc
end


function gain_chain(chain)
    intchain = similar(chain.instrument.lgz)
    for i in eachindex(intchain)
        S = intchain.sites[i]
        lgz = chain.instrument.lgz[i]
        lgμ = chain.instrument.lgμ[S=S][1]
        lgσ = chain.instrument.lgσ[S=S][1]
        intchain[i] = exp(lgμ + exp(lgσ) * lgz)
    end
    return intchain
end
