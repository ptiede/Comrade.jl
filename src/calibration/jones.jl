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

    for i in eachindex(bls)
        s1, s2 = bls[i]
        ind1 = findfirst(==(s1), stats)
        ind2 = findfirst(==(s2), stats)

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
    out1 = zero(first(m1))
    out2 = zero(first(m2))
    _allmul!(out1, out2, m1, m2)
    # out1 = reduce(.*, m1)
    # out2 = reduce(.*, m2)
    return JonesPairs(out1, out2)
end

function _allmul!(out1, out2, m1, m2)
    for i in eachindex(out1, out2)
        out1[i] = mapreduce(x->getindex(x, i), Base.:*, m1)
        out2[i] = mapreduce(x->getindex(x, i), Base.:*, m2)
    end
    return nothing
end

using Enzyme
function ChainRulesCore.rrule(::typeof(_allmul), m1, m2)
    out = _allmul(m1, m2)
    pm1 = ProjectTo(m1)
    pm2 = ProjectTo(m2)
    function _allmul_pullback(Δ)
        Δm1 = zero(out.m1)
        Δm1 .= unthunk(Δ.m1)
        Δm2 = zero(out.m1)
        Δm2 .= unthunk(Δ.m2)
        dm1 = zero.(m1)
        dm2 = zero.(m2)

        out1 = zero(first(m1))
        out2 = zero(first(m2))
        autodiff(_allmul!, Const, Duplicated(out1, Δm1), Duplicated(out2, Δm2), Duplicated(m1, dm1), Duplicated(m2, dm2))
        return NoTangent(), pm1(dm1), pm2(dm2)
    end
    return out, _allmul_pullback
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

    return JonesCache{typeof(m1),typeof(s), typeof(bls), typeof(times)}(m1,m2, s, bls, times)
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
    gm2 = gmat(g1, g2, jcache.m2)
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

export jonesStokes
"""
    jonesStokes(g1::AbstractArray, gcache::JonesCache)

Construct the JonesPairs for the stokes I image only. That is, we only need to
pass a single vector corresponding to the gain for the stokes I visibility. This is
for when you only want to image Stokes I.

# Warning
In the future this functionality may be removed when stokes I fitting is replaced with the
more correct `trace(coherency)`, i.e. RR+LL for a circular basis.
"""
function jonesStokes(g, gcache)
    return JonesPairs(gcache.m1*g, gcache.m2*g)
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

function TransformCache(obs::EHTObservation; add_fr=true, ehtim_fr_convention=true, ref::PolBasis=CirBasis())
    T1 = StructArray(map(x -> basis_transform(ref, x[1]), obs.data.polbasis))
    T2 = StructArray(map(x -> basis_transform(ref, x[2]), obs.data.polbasis))
    if add_fr
        field_rotations = extract_FRs(obs; ehtim_fr_convention)
        T1 .= field_rotations.m1.*T1
        T2 .= field_rotations.m2.*T2
    end
    return TransformCache{typeof(T1), typeof(ref)}(T1, T2, ref)
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

visanalytic(::Type{<:JonesModel{J,M}}) where {J,M} = visanalytic(M)
imanalytic(::Type{<:JonesModel{J,M}}) where {J,M} = imanalytic(M)
ispolarized(::Type{<:JonesModel{J,M}}) where {J,M} = ispolarized(M)




function JonesModel(jones, model, tcache::TransformCache)
    return JonesModel(jones, model, tcache.refbasis)
end

function JonesModel(jones, model)
    return JonesModel(jones, model, CirBasis())
end

function _visibilities(model::JonesModel{J,M,B}, u, v, time, freq) where {J,M,B}
    vis = _visibilities(model.model, u, v, time, freq)
    coh = _coherency(vis, B)
    return corrupt(coh, model.jones.m1, model.jones.m2)
end

function corrupt!(vnew, vis::AbstractArray, j1, j2)
    # @assert length(vis) == length(j1) "visibility vector and jones pairs have mismatched dimensions!"
    vnew .= j1 .* vis .* adjoint.(j2)
    return nothing
end

function corrupt(vis, j1, j2)
    vnew = j1 .* vis .* adjoint.(j2)
    return vnew
end


# function ChainRulesCore.rrule(::typeof(corrupt), vis, j1, j2)
#     out = corrupt(vis, j1, j2)
#     pvis = ProjectTo(vis)
#     pj1 = ProjectTo(j1)
#     pj2 = ProjectTo(j2)
#     function _corrupt_pullback(Δ)
#         Δout = similar(out)
#         println(typeof(Δ))
#         println(typeof(out))
#         Δout .= unthunk(Δ)
#         Δvis = zero(out)
#         v2 = zero(out)
#         v2 .= vis
#         Δj1 = zero(j1)
#         Δj2 = zero(j2)
#         out = similar(out)
#         autodiff(corrupt!, Const, Duplicated(out, Δout), Duplicated(v2, Δvis), Duplicated(j1, Δj1), Duplicated(j2, Δj2))
#         return NoTangent(), pvis(Δvis), pj1(Δj1), pj2(Δj2)
#     end
#     return out, _corrupt_pullback
# end

function _coherency(vis::AbstractArray{T}, ::Type{B}) where {T<:AbstractArray,B}
    return CoherencyMatrix{B,B}.(vis)
end

function _coherency(vis::AbstractArray{<:Number}, ::Type{B}) where {B}
    return vis
end

function ChainRulesCore.rrule(::typeof(_coherency), vis::AbstractArray{<:AbstractArray}, ::Type{CirBasis})
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


function extract_FRs(obs::EHTObservation; ehtim_fr_convention=true)

    # read elevation angles for each station
    config = arrayconfig(obs)
    el1 = StructArrays.component(config.data.elevation, 1)
    el2 = StructArrays.component(config.data.elevation, 2)

    # read parallactic angles for each station
    par1 = StructArrays.component(config.data.parallactic, 1)
    par2 = StructArrays.component(config.data.parallactic, 2)


    # get ehtobservation array info
    tarr  = config.tarr
    ants  = tarr.sites
    elevs = tarr.fr_elevation
    pars  = tarr.fr_parallactic
    offs  = tarr.fr_offset

    # get station names
    bls = config.data.baseline
    ant1 = first.(bls)
    ant2 = last.(bls)

    # get multiplicative prefactors
    f_el1  = zero(el1)
    f_par1 = zero(par1)
    f_off1 = zero(el1)
    f_el2  = zero(el2)
    f_par2 = zero(par2)
    f_off2 = zero(el2)
    for i in eachindex(ant1)
        ind1 = findall(==(ant1[i]), ants) |> first
        ind2 = findall(==(ant2[i]), ants) |> first

        f_el1[i]  = elevs[ind1]
        f_el2[i]  = elevs[ind2]

        f_par1[i] = pars[ind1]
        f_par2[i] = pars[ind2]

        f_off1[i] = offs[ind1]
        f_off2[i] = offs[ind2]
    end
    # combine to get field rotations for each station
    FR1 = (f_el1 .* el1) .+ (f_par1 .* par1) .+ f_off1
    FR2 = (f_el2 .* el2) .+ (f_par2 .* par2) .+ f_off2

    if ehtim_fr_convention
        FR1 .*= 2
        FR2 .*= 2
    end
    S = Complex{eltype(FR1)}
    offdiag = fill(zero(eltype(FR1)), length(FR1))
    jF1 = StructArray{SMatrix{2,2,S,4}}((cis.(-FR1), offdiag, offdiag, cis.(FR1)))
    jF2 = StructArray{SMatrix{2,2,S,4}}((cis.(-FR2), offdiag, offdiag, cis.(FR2)))

    return JonesPairs(jF1, jF2)
end
