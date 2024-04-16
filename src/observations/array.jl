
"""
    $(TYPEDEF)

This defined the abstract type for an array configuration. Namely, baseline
times, SEFD's, bandwidth, observation frequencies, etc.
"""
abstract type AbstractArrayConfiguration{F<:AbstractBaselineDatum} <: AbstractVLBITable{F} end
arrayconfig(c::AbstractArrayConfiguration) = c
build_datum(c::AbstractArrayConfiguration, i) = datatable(c)[i]

function Base.getindex(config::F, i::AbstractVector) where {F<:AbstractArrayConfiguration}
    newconftable = datatable(config)[i]
    newconfig = @set config.datatable = newconftable
    return newconfig
end


"""
    sites(d::AbstractArrayConfiguration)

Get all the sites in a observation. The result is a vector of symbols.
"""
function sites(d::AbstractArrayConfiguration)
    bl = d[:sites]
    s1 = first.(bl)
    s2 = last.(bl)
    return sort(unique(vcat(s1, s2)))
end


"""
    $(SIGNATURES)

Get the u, v, time, freq of the array as a tuple.
"""
function getuvtimefreq(ac::AbstractArrayConfiguration)
    u = ac[:U]
    v = ac[:V]
    t = ac[:T]
    ν = ac[:F]
    return (U=u, V=v, T=t, F=ν)
end


"""
    beamsize(ac::AbstractArrayConfiguration)

Calculate the approximate beam size of the array `ac` as the inverse of the longest baseline
distance.
"""
beamsize(ac::AbstractArrayConfiguration) = inv(mapreduce(hypot, max, ac[:U], ac[:V]))



"""
    $(TYPEDEF)

A single datum of an `ArrayConfiguration`
"""
struct EHTArrayBaselineDatum{T,E,V} <: AbstractBaselineDatum
    """
    u position of the data point in λ
    """
    U::T
    """
    v position of the data point in λ
    """
    V::T
    """
    time of the data point in (Hr)
    """
    T::T
    """
    frequency of the data point (Hz)
    """
    F::T
    """
    Sites codes of the baseline (u,v)
    """
    sites::Tuple{Symbol, Symbol}
    """
    Polarization basis
    """
    polbasis::E
    """
    elevation of baseline
    """
    elevation::Tuple{V,V}
    """
    parallactic angle of baseline
    """
    parallactic::Tuple{V,V}
    function EHTArrayBaselineDatum(u, v, time, freq, sites, polbasis, elevation, parallactic)
        tt, ft, ut, vt = promote(time, freq, u, v)
        T = typeof(tt)
        V = typeof(elevation[1])
        E = typeof(error)
        return new{T,E,V}(ut, vt, tt, ft, sites, polbasis, elevation, parallactic)
    end
end


"""
    $(TYPEDEF)

Table that specified pertinent details about the EHT array during an observation.
These are typically items that are known before the observation is made.

# Fields
$(FIELDS)
"""
Base.@kwdef struct EHTArrayConfiguration{A<:EHTArrayBaselineDatum,F,T,S,D<:AbstractArray{A}} <: AbstractArrayConfiguration{A}
    """
    Observing bandwith (Hz)
    """
    bandwidth::F
    """
    Telescope array file
    """
    tarr::T
    """
    Scan times
    """
    scans::S
    """
    modified julia date of the observation
    """
    mjd::Int
    """
    RA of the observation in J2000 (deg)
    """
    ra::F
    """
    DEC of the observation in J2000 (deg)
    """
    dec::F
    """
    Common source name
    """
    source::Symbol
    """
    Time zone used.
    """
    timetype::Symbol = :UTC
    """
    A struct array of `EHTArrayBaselineDatum`
    """
    datatable::D
end

function Base.show(io::IO, config::EHTArrayConfiguration)
    println(io, "EHTArrayConfiguration:")
    println(io, "  source:      ", config.source)
    println(io, "  mjd:         ", config.mjd)
    println(io, "  frequencies: ", unique(config[:F]))
    println(io, "  bandwidth:   ", config.bandwidth)
    println(io, "  sites:       ", sites(config))
    print(io, "  nsamples:    ", length(config))
end




"""
    $(TYPEDEF)
Array config file for closure quantities. This stores the design matrix `designmat`
that transforms from visibilties to closure products.

# Fields
$(FIELDS)
"""
struct ClosureConfig{F, A<:AbstractArrayConfiguration{F},D} <: AbstractArrayConfiguration{F}
    """Array configuration for visibilities"""
    ac::A
    """Closure design matrix"""
    designmat::D

    function ClosureConfig(ac::AbstractArrayConfiguration{F}, dmat) where {F}
        A = typeof(ac)
        sdmat = blockdiag(sparse.(dmat)...)
        D = typeof(sdmat)
        return new{F,A,D}(ac, sdmat)
    end
end
datatable(c::ClosureConfig) = datatable(arrayconfig(c))
arrayconfig(c::ClosureConfig) = getfield(c, :config)
build_datum(c::ClosureConfig, i) = datatable(c)[i]

function getuvtimefreq(ac::ClosureConfig)
    return getuvtimefreq(ac.ac.config)
end
