
"""
    $(TYPEDEF)

This defined the abstract type for an array configuration. Namely, baseline
times, SEFD's, bandwidth, observation frequencies, etc.
"""
abstract type AbstractArrayConfiguration{F<:AbstractBaselineDatum} <: AbstractVLBITable{F} end
arrayconfig(c::AbstractArrayConfiguration) = c
build_datum(c::AbstractArrayConfiguration, i) = datatable(c)[i]





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
struct EHTArrayBaselineDatum{T,P,V} <: AbstractBaselineDatum
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
    polbasis::P
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
        E = typeof(polbasis)
        return new{T,E,V}(ut, vt, tt, ft, sites, polbasis, elevation, parallactic)
    end
end

function flipbaseline(d::EHTArrayBaselineDatum)
    return EHTArrayBaselineDatum(-d.U, -d.V, d.T, d.F,
                                (d.sites[2], d.sites[1]),
                                (d.polbasis[2], d.polbasis[1]),
                                (d.elevation[2], d.elevation[1]),
                                (d.parallactic[2], d.parallactic[1]))
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

function VLBISkyModels.rebuild(c::EHTArrayConfiguration, table)
    newconfig    = EHTArrayConfiguration(c.bandwidth, c.tarr, c.scans,
                                         c.mjd, c.ra, c.dec, c.source,
                                         c.timetype, table)
    return newconfig
end

struct DesignMatrix{T, N, M<:AbstractSparseMatrix{T, <:Integer}, I} <: AbstractMatrix{T}
    matrix::M
    nz::I
    function DesignMatrix(m::AbstractSparseMatrix{T}) where {T}
        nz = findnz(m)
        ncl = length(findall(==(1), nz[1]))
        return new{T, ncl, typeof(m),typeof(nz)}(m, nz)
    end
end
Base.parent(m::DesignMatrix) = m.matrix
Base.getindex(m::DesignMatrix, i::Int) = getindex(m.matrix, i)
Base.size(m::DesignMatrix) = size(m.matrix)
Base.IndexStyle(::Type{<:DesignMatrix{X,N,M}}) where {X,N,M} = Base.IndexStyle(M)
Base.getindex(m::DesignMatrix, I::Vararg{Int,N}) where {N} = getindex(m.matrix, I...)
# Base.setindex!(m::DesignMatrix, v, i::Int) = setindex!(m.matrix, v, i)
# Base.setindex!(m::DesignMatrix, v, i::Vararg{Int, N}) where {N} = setindex!(m.matrix, v, i...)

Base.similar(m::DesignMatrix, ::Type{S}, dims::Dims) where {S} = similar(m.matrix, S, dims)

LinearAlgebra.mul!(out::AbstractArray, m::DesignMatrix, x::AbstractArray) = mul!(out, parent(m), x)
Base.show(io::IO, mime::MIME"text/plain",   m::DesignMatrix) = show(io, mime, m.matrix)


"""
    $(TYPEDEF)
Array config file for closure quantities. This stores the design matrix `designmat`
that transforms from visibilties to closure products.

# Fields
$(FIELDS)
"""
struct ClosureConfig{F, A<:AbstractArrayConfiguration{F}, D, V, E} <: AbstractArrayConfiguration{F}
    """Array configuration for visibilities"""
    ac::A
    """Closure design matrix"""
    designmat::D
    """visibilities to closure design matrix"""
    vis::V
    """visibility noises to closure design matrix"""
    noise::E
    function ClosureConfig(ac::AbstractArrayConfiguration{F}, dmat::Vector{<:AbstractMatrix}, vis, noise) where {F}
        A = typeof(ac)
        sdmat = DesignMatrix(blockdiag(sparse.(dmat)...))
        D = typeof(sdmat)
        return new{F,A,D, typeof(vis), typeof(noise)}(ac, sdmat, vis, noise)
    end
    function ClosureConfig(ac::AbstractArrayConfiguration{F}, dmat::DesignMatrix, vis, noise) where {F}
        A = typeof(ac)
        D = typeof(dmat)
        return new{F,A,D, typeof(vis), typeof(noise)}(ac, dmat, vis, noise)
    end
end
arrayconfig(c::ClosureConfig) = getfield(c, :ac)
Base.length(c::ClosureConfig) = size(designmat(c), 1)
Base.firstindex(c::ClosureConfig) = 1
Base.lastindex(c::ClosureConfig) = length(c)

function Base.getindex(c::ClosureConfig, i::AbstractVector)
    dmat = designmat(c)
    dmat2 = DesignMatrix(dmat[i, :])
    return ClosureConfig(arrayconfig(c), dmat2, getfield(c, :vis), getfield(c, :noise))
end

function Base.view(c::ClosureConfig, i::AbstractVector)
    dmat = designmat(c)
    dmat2 = DesignMatrix(parent(dmat)[i, :])
    return ClosureConfig(arrayconfig(c), dmat2, getfield(c, :vis), getfield(c, :noise))
end


Base.propertynames(c::ClosureConfig) = propertynames(arrayconfig(c))
function Base.getproperty(c::ClosureConfig, p::Symbol)
    getproperty(arrayconfig(c), p)
end
designmat(c::ClosureConfig) = getfield(c, :designmat)

function build_datum(arr::ClosureConfig{F, A, <:DesignMatrix{T, N}}, i::Int) where {F, A, T, N}
    arrvis = arrayconfig(arr)

    nz = designmat(arr).nz
    inds = range(firstindex(nz[2])+N*(i-1); length=N)
    J = @view(nz[2][inds])
    V = @view(nz[3][inds])
    bls = ntuple(Val(N)) do i
        return V[i] == 1.0 ? arrvis[J[i]] : flipbaseline(arrvis[J[i]])
    end
    return bls
end

function datatable(c::ClosureConfig)
    StructArray((build_datum(c, i) for i in 1:length(c)), unwrap=(T->(T<:Tuple || T<:AbstractBaselineDatum)))
end

function sites(c::ClosureConfig)
    sites(arrayconfig(c))
end



function getuvtimefreq(ac::ClosureConfig)
    return getuvtimefreq(arrayconfig(ac))
end

function Base.show(io::IO, config::ClosureConfig)
    println(io, "ClosureConfig:")
    println(io, "  source:      ", config.source)
    println(io, "  mjd:         ", config.mjd)
    println(io, "  bandwidth:   ", config.bandwidth)
    println(io, "  sites:       ", sites(arrayconfig(config)))
    print(io, "  nclosures:    ", length(config))
end
