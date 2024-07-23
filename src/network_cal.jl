"""
    NetworkCalibrationSkyModel(zbl_flux, netcal_bl)

Constructs a SkyModel that represents what is done for network calibration.
Network calibration requires a special model since there is no actual image used. Instead
we assume that the sky has some total flux given by `zbl_flux` and the rest of the
amplitudes are the actual model parameters.

!!! note
    By default we will assume that the amplitudes are a flat
    prior from [0, zbl_flux] to be maximally permissive.

!!! note
    We need a special skymodel for network calibration since the model is not an image but
    rather we directly fit the visibility amplitudes for non-intrasite baselines.

# Arguments

 - `zbl_flux`  : The apriori measured total flux of the object.
 - `netcal_bl` : The baselines that are considered to be co-located for network calibration.
"""
Base.@kwdef struct NetworkCalSkyModel{Z<:Real, B} <: AbstractSkyModel
    zbl_flux::Z
    netcal_bl::B
end

# From LogExpFunctions.jl
@inline _logistic_bounds(::Float16) = (Float16(-16.64), Float16(7.625))
@inline _logistic_bounds(::Float32) = (-103.27893f0, 16.635532f0)
@inline _logistic_bounds(::Float64) = (-744.4400719213812, 36.7368005696771)

@inline function elogistic(x::Union{Float16, Float32, Float64})
    e = @inline exp(x)
    lower, upper = _logistic_bounds(x)
    return x < lower ? zero(x) : x > upper ? one(x) : e / (one(x) + e)
end

function set_array(m::NetworkCalSkyModel, array::AbstractArrayConfiguration)
    dtbl = datatable(array)
    sites = dtbl.sites

    netcalset = m.netcal_bl
    intrainds = findall(x->Set(x)∈Set.(netcalset), sites)
    fixvals = fill(0.0, length(intrainds))
    ampinds = setdiff(eachindex(sites), intrainds)
    dists = Distributions.MvNormal(Diagonal(fill(1.78^2, length(ampinds))))

    d = PartiallyConditionedDist(dists, ampinds, intrainds, fixvals)
    skypr = d
    f = let zblflux=m.zbl_flux, intrainds=intrainds, ampinds=ampinds
        x->(x[ampinds] .= 4 .*zblflux.*elogistic.(@view(x[ampinds])); x[intrainds] .= zblflux; y)
    end
    g = imagepixels(μas2rad(100.0), μas2rad(100.0), 256, 256)
    return ObservedSkyModel(m, FourierDualDomain(g, array, NFFTAlg()), f), skypr
end

function idealvisibilities(m::ObservedSkyModel{<:NetworkCalSkyModel}, x)
    return m.metadata(x.sky)
end

function skymodel(m::ObservedSkyModel{<:NetworkCalSkyModel}, x)
    return m.metadata(x)
end

function prepare_netcal_data(obs::EHTObservationTable{<:EHTVisibilityAmplitudeDatum}, netcal_bl...)
    S = Set(Iterators.flatten(netcal_bl))
    array = arrayconfig(obs)
    inds = findall(x->(x[1]∈S || x[2]∈S), datatable(array).sites)
    # We find all baselines that are connected to our network calibration baselines
    return obs[inds]
end
