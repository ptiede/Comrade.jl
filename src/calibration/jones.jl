export JonesCache, TrackSeg, ScanSeg, IntegSeg, jonesG, jonesD, jonesT,
       TransformCache, JonesModel

abstract type ObsSegmentation end

# Track is for quantities that remain static across an entire observation
struct TrackSeg <: ObsSegmentation end

# Scan is for quantities that are constant across a scan
struct ScanSeg <: ObsSegmentation end

# Integration is for quantities that change every integration time
struct IntegSeg <: ObsSegmentation end




function LinearAlgebra.mul!(y::AbstractArray, M::DesignMatrix, x::AbstractArray)
    LinearAlgebra.mul!(y, M.matrix, x)
end


struct JonesCache{D,S<:ObsSegmentation, ST, Ti}
    """
    Design matrix for the first station
    """
    m1::D
    """
    Design matrix for the second station
    """
    m2::D
    """
    Segmentation scheme for this cache
    """
    seg::S
    """
    station codes
    """
    stations::ST
    """
    times
    """
    times::Ti
end

stations(j::JonesCache) = j.stations
ChainRulesCore.@non_differentiable stations(j::JonesCache)

function JonesCache(obs::EHTObservation, s::TrackSeg)

    # extract relevant observation info
    times = obs[:T]
    bls = obs[:baseline]
    stats = stations(obs)

    # initialize vectors containing row and column index locations
    rowInd1 = Int[]
    colInd1 = Int[]
    rowInd2 = Int[]
    colInd2 = Int[]

    for i in 1:length(times)

        s1, s2 = bls[i]
        ind1 = findfirst(x -> x==s1, stats)
        ind2 = findfirst(x -> x==s2, stats)

        append!(colInd1,ind1)
        append!(colInd2,ind2)
        append!(rowInd1,i)
        append!(rowInd2,i)
    end

    # populate sparse design matrices
    z = fill(1.0, length(rowInd1))
    m1 = sparse(rowInd1, colInd1, z, length(times), length(stats))
    m2 = sparse(rowInd2, colInd2, z, length(times), length(stats))

    ttimes = fill(times[begin], length(stats))

    return JonesCache{typeof(m1),typeof(s), typeof(stats), typeof(ttimes)}(m1,m2,s, stats, ttimes)
end

function JonesCache(obs::EHTObservation, s::ScanSeg)

    # extract relevant observation info
    times = obs[:T]
    bls = obs[:baseline]
    gaintime, gainstat = gain_stations(scantable(obs))
    gts = collect(zip(gaintime, gainstat))

    # initialize vectors containing row and column index locations
    rowInd1 = Int[]
    colInd1 = Int[]
    rowInd2 = Int[]
    colInd2 = Int[]

    for i in 1:length(times)

        t = times[i]
        s1, s2 = bls[i]
        ind1 = findall(x -> ((x[1]==t) && (x[2]==s1)), gts)
        ind2 = findall(x -> ((x[1]==t) && (x[2]==s2)), gts)

        append!(colInd1,ind1)
        append!(colInd2,ind2)
        append!(rowInd1,fill(i, length(ind1)))
        append!(rowInd2,fill(i, length(ind2)))
    end

    # populate sparse design matrices
    z = fill(1.0, length(rowInd1))
    m1 = sparse(rowInd1, colInd1, z, length(times), length(gaintime))
    m2 = sparse(rowInd2, colInd2, z, length(times), length(gaintime))

    return JonesCache{typeof(m1),typeof(s),  typeof(gainstat), typeof(gaintime)}(m1,m2,s, gainstat, gaintime)
end

struct JonesPairs{T, M1<:AbstractVector{T}, M2<:AbstractVector{T}}
    m1::M1
    m2::M2
end

function Base.:*(x::JonesPairs, y::JonesPairs...)
    m1 = map(x->getproperty(x, :m1), (x,y...))
    m2 = map(x->getproperty(x, :m2), (x,y...))
    return _allmul(m1, m2)
end

function _allmul(m1, m2)
    m1 = reduce(.*, m1)
    m2 = reduce(.*, m2)
    return JonesPairs(m1, m2)
end


# function JonesPairs(m1::AbstractVector{T}, m2::AbstractVector{T}) where {T}
#     @assert length(m1) == length(m2) "m1 and m2 must have the same length"
#     @assert firstindex(m1) == firstindex(m2) "m1 and m2 must be indexed the same"
#     return JonesPairs{T, typeof(m1), typeof(m2)}(m1, m2)
# end

Base.length(j::JonesPairs) = length(j.m1)
Base.size(j::JonesPairs) = (length(j),)
Base.getindex(j::JonesPairs, i::Int) = (j.m1[i], j.m2[i])
function Base.setindex!(j::JonesPairs, X, i::Int)
    j.m1[i] = X[1]
    j.m2[i] = X[2]
    return j
end
Base.IndexStyle(::Type{JonesPairs}) = IndexLinear()
Base.similar(j::JonesPairs, ::Type{S}, dims::Dims{1}) where {S} = JonesPairs(similar(j.m1, S, dims), similar(j.m2, S, dims))
Base.firstindex(j::JonesPairs) = firstindex(j.m1)
Base.lastindex(j::JonesPairs) = lastindex(j.m1)

struct JonesStyle{M1,M2} <: Broadcast.AbstractArrayStyle{1} end
JonesStyle{M1,M2}(::Val{1}) where {M1,M2} = JonesStyle{M1,M2}()
Base.BroadcastStyle(::Type{<:JonesPairs{T,M1,M2}}) where {T,M1,M2} = JonesStyle{M1,M2}()

function Base.similar(bc::Broadcast.Broadcasted{JonesStyle{M1,M2}}, ::Type{ElType}) where {M1,M2,ElType}
    A = find_js(bc)
    n = length(A.m1)
    return JonesPairs(similar(M1, Eltype, n), similar(M2, Eltype, n))
end

"`A = find_aac(As)` returns the first ArrayAndChar among the arguments."
find_js(bc::Base.Broadcast.Broadcasted) = find_js(bc.args)
find_js(args::Tuple) = find_js(find_js(args[1]), Base.tail(args))
find_js(x) = x
find_js(::Tuple{}) = nothing
find_js(a::JonesPairs, rest) = a
find_js(::Any, rest) = find_js(rest)




function JonesCache(obs::EHTObservation, s::IntegSeg)

    # extract relevant observation info
    times = obs[:T]
    bls = obs[:baseline]
    stats = stations(obs)

    # organize time-station info
    tuniq = unique(times)
    tbl = Tuple{eltype(tuniq),eltype(stats)}[]
    for i in 1:length(tuniq)
        t = tuniq[i]
        ind = findall(x -> (x==t), times)
        s1 = getindex.(bls[ind],1)
        s2 = getindex.(bls[ind],2)
        statshere = unique(hcat(s1,s2))
        for j in 1:length(statshere)
            push!(tbl,(t,statshere[j]))
        end
    end

    # initialize vectors containing row and column index locations
    rowInd1 = Int[]
    colInd1 = Int[]
    rowInd2 = Int[]
    colInd2 = Int[]

    for i in 1:length(times)

        t = times[i]
        s1, s2 = bls[i]
        ind1 = findall(x -> ((x[1]==t) && (x[2]==s1)), tbl)
        ind2 = findall(x -> ((x[1]==t) && (x[2]==s2)), tbl)

        append!(colInd1,ind1)
        append!(colInd2,ind2)
        append!(rowInd1,fill(i, length(ind1)))
        append!(rowInd2,fill(i, length(ind2)))
    end

    # populate sparse design matrices
    z = fill(1.0, length(rowInd1))
    m1 = sparse(rowInd1, colInd1, z, length(times), length(tbl))
    m2 = sparse(rowInd2, colInd2, z, length(times), length(tbl))

    return JonesCache{typeof(m1),typeof(s), typeof(bls), typeof(times)}(m1,m2,s, bls, times)
end


function apply_design(gmat::T, jcache::JonesCache) where {T}
    return JonesPairs(jcache.m1*gmat, jcache.m2*gmat)
end

# GMat(g1::T, g2::T) where {T} = SMatrix{2,2,T}(g1, zero(T), zero(T), g2)
function gmat(g1, g2, m)
   S = eltype(g1)
   gs1 = m*g1
   gs2 = m*g2
   n = length(gs1)
   offdiag = fill(zero(S), n)
   StructArray{SMatrix{2,2,S,4}}((gs1, offdiag, offdiag, gs2))
end
function jonesG(g1, g2,jcache::JonesCache)
    gm1 = gmat(g1, g2, jcache.m1)
    gm2 = gmat(g1, g2, jcache.m1)
    return JonesPairs(gm1, gm2)
end

# DMat(d1::T, d2::T) where {T} = SMatrix{2,2,T}(one(T), d2, d1, one(T))
function dmat(d1, d2, m)
    S = eltype(d1)
    ds1 = m*d1
    ds2 = m*d2
    n = length(ds1)
    unit = fill(one(S), n)
    return StructArray{SMatrix{2,2,S,4}}((unit, ds2, ds1, unit))
end

function jonesD(d1::T,d2::T,jcache::JonesCache) where {T}
    dm1 = dmat(d1, d2, jcache.m1)
    dm2 = dmat(d1, d2, jcache.m2)
    return JonesPairs(dm1, dm2)
end





struct TransformCache{M, B<:PolBasis}
    """
    Transform matrices for the first stations
    """
    T1::M
    """
    Transform matrices for the second stations
    """
    T2::M
    """
    Reference polarization basis
    """
    refbasis::B
end

function TransformCache(obs::EHTObservation; ref::PolBasis=CirBasis())
    T1 = StructArray(map(x -> basis_transform(ref, x[1]), obs.data.polbasis))
    T2 = StructArray(map(x -> basis_transform(ref, x[2]), obs.data.polbasis))
    return TransformCache{typeof(T1),typeof(ref)}(T1, T2, ref)
end

jonesT(tcache::TransformCache) = JonesPairs(tcache.T1, tcache.T2)



struct JonesModel{J, M, B<:PolBasis} <: RIMEModel
    """
    Cache containing the specific Jones matrices that are to be applied to the visibilities.
    """
    jones::J
    """
    Base model that will be used to compute the uncorrupted visibilities.
    """
    model::M
    """
    Reference basis that converts from StokesParams to CoherencyMatrix
    """
    refbasis::B
end

function JonesModel(jones, model, tcache::TransformCache)
    return JonesModel(jones, model, tcache.refbasis)
end

function _visibilities(model::JonesModel{J,M,B}, u, v, time, freq) where {J,M,B}
    vis = _visibilities(model.model, u, v, time, freq)
    coh = _coherency(vis, B)
    return corrupt(coh, model.jones)
end

function corrupt(vis::AbstractArray, jones)
    # @assert length(vis) == length(jones) "visibility vector and jones pairs have mismatched dimensions!"
    vnew = jones.m1 .* vis .* adjoint.(jones.m2)
    return vnew
end

function _coherency(vis, ::Type{B}) where {B}
    return CoherencyMatrix{B,B}.(vis)
end

function ChainRulesCore.rrule(::typeof(_coherency), vis, ::Type{CirBasis})
    coh  = _coherency(vis, CirBasis)
    pd = ProjectTo(vis)
    function _coherency_pullback(Δd)
        Δvis = zero(vis)
        Δ = unthunk(Δd)
        for i in eachindex(vis)
            ΔRR = Δ[i][1,1]
            ΔLR = Δ[i][2,1]
            ΔRL = Δ[i][1,2]
            ΔLL = Δ[i][2,2]
            Δvis.I[i] = complex(ΔRR + ΔLL)
            Δvis.Q[i] = complex(ΔLR + ΔRL)
            Δvis.U[i] = 1im*(ΔLR - ΔRL)
            Δvis.V[i] = complex(ΔRR - ΔLL)
        end
        return NoTangent(), pd(Δvis), NoTangent()
    end
    return coh, _coherency_pullback
end
