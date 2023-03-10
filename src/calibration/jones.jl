export JonesCache, TrackSeg, ScanSeg, FixedSeg, IntegSeg, jonesG, jonesD, jonesT,
       TransformCache, JonesModel, jonescache

"""
    $(TYPEDEF)

The data segmentation scheme to use. This is important for constructing a [`JonesCache`](@ref)
"""
abstract type ObsSegmentation end

# Track is for quantities that remain static across an entire observation
"""
    $(TYPEDEF)

Data segmentation such that the quantity is constant over a `track`, i.e., the observation "night".
"""
struct TrackSeg <: ObsSegmentation end

# Scan is for quantities that are constant across a scan
"""
    $(TYPEDEF)

Data segmentation such that the quantity is constant over a `scan`.

## Warning
Currently we do not explicity track the telescope scans. This will be fixed in a future version.
Right now `ScanSeg` and `TrackSeg` are the same
"""
struct ScanSeg{S} <: ObsSegmentation end

ScanSeg(segmented=false) = ScanSeg{false}()

# Integration is for quantities that change every integration time
"""
    $(TYPEDEF)

Data segmentation such that the quantity is constant over a correlation integration.
"""
struct IntegSeg{S} <: ObsSegmentation end

"""
    $(TYPEDEF)

Enforces that the station calibraton value will have a fixed `value`. This is most
commonly used when enforcing a reference station for gain phases.
"""
struct FixedSeg{T} <: ObsSegmentation
    value::T
end


function LinearAlgebra.mul!(y::AbstractArray, M::DesignMatrix, x::AbstractArray)
    LinearAlgebra.mul!(y, M.matrix, x)
end

struct AffineDesignMatrix{M, B}
    mat::M
    b::B
end

Base.size(d::AffineDesignMatrix) = size(d.mat)
Base.size(d::AffineDesignMatrix, i::Int) = size(d.mat, i)
Base.copy(d::AffineDesignMatrix) = AffineDesignMatrix(copy(d.mat), copy(d.b))

# function Base.:*(A::AffineDesignMatrix, v::AbstractVector)
#     return muladd(A.mat, v, A.b)
# end

function LinearAlgebra.mul!(y::AbstractArray, M::AffineDesignMatrix, x::AbstractArray)
    # println(size(M.mat))
    # println(size(x))
    mul!(y, M.mat, x)
    y .+= M.b
    return y
end

abstract type AbstractJonesCache end

"""
    $(TYPEDEF)

Holds the ancillary information for a the design matrix cache for Jones matrices. That is,
it defines the cached map that moves from model visibilities to the corrupted voltages
that are measured from the telescope.

# Fields
$(FIELDS)
"""
struct JonesCache{D1, D2, S, Sc} <: AbstractJonesCache
    """
    Design matrix for the first station
    """
    m1::D1
    """
    Design matrix for the second station
    """
    m2::D2
    """
    Segmentation schemes for this cache
    """
    seg::S
    """
    Gain Schema
    """
    schema::Sc
end

"""
    $(TYPEDEF)

Holds the ancillary information for a the design matrix cache for Jones matrices. That is,
it defines the cached map that moves from model visibilities to the corrupted voltages
that are measured from the telescope. This uses a segmented decomposition so that the
gain at a single timestamp is the sum of the previous gains. In this formulation the
gains parameters are the segmented gain offsets from timestamp to timestamp

# Fields
$(FIELDS)
"""
struct SegmentedJonesCache{D,S<:ObsSegmentation, ST, Ti} <: AbstractJonesCache
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

stations(j::AbstractJonesCache) = j.stations
ChainRulesCore.@non_differentiable stations(j::AbstractJonesCache)


struct GainSchema{S,T,G}
    sites::S
    times::T
    gts::G
    function GainSchema(sites, times)
        gts = collect(zip(times, sites))
        return new{typeof(sites), typeof(times), typeof(gts)}(sites, times, gts)
    end
end

function gain_schema(segmentation::NamedTuple, obs::EHTObservation)
    st = scantable(obs)
    times = eltype(obs[:T])[]
    sites = Symbol[]
    for i in 1:length(st)
        s = stations(st[i])
        t = st[i].time
        for j in eachindex(s)
            append_time_site!(times, sites, s[j], t, getproperty(segmentation, s[j]))
        end
    end
    return GainSchema(sites, times)
end

function gain_schema(segmentation::ObsSegmentation, obs::EHTObservation)
    sites = Tuple(stations(obs))
    return gain_schema(NamedTuple{sites}(Tuple(segmentation for _ in sites)), obs)
end

function append_time_site!(times, sites, site, t, ::ScanSeg)
    push!(sites, site)
    push!(times, t)
end

function append_time_site!(times, sites, site, t, ::IntegSeg)
    push!(sites, site)
    push!(times, t)
end

function append_time_site!(times, sites, site, t, ::TrackSeg)
    # Check is the site is already in the list, if it isn't add it to the schema
    !(site ∈ sites) && (push!(sites, site); push!(times, zero(t)))
end

function append_time_site!(times, sites, site, t, ::FixedSeg)
    nothing
end



function fill_designmat!(
        colInd::AbstractVector{Int}, rowInd::AbstractVector{Int},
        vecInd::AbstractVector{Int}, valInd::AbstractVector,
        ::TrackSeg, vis_ind, site, time, schema::GainSchema
        )
    stats = schema.sites
    ind = findfirst(==(site), stats)
    append!(colInd, [ind])
    append!(rowInd, [vis_ind])
end

function fill_designmat!(
    colInd::AbstractVector{Int}, rowInd::AbstractVector{Int},
    vecInd::AbstractVector{Int}, vals::AbstractVector,
    ::ScanSeg{false}, vis_ind, site, time, schema::GainSchema
    )

    ind = findall(x->((x[1]==time) && (x[2]==site)), schema.gts)
    append!(colInd, ind)
    append!(rowInd, fill(vis_ind, length(ind)))
end

function fill_designmat!(
    colInd::AbstractVector{Int}, rowInd::AbstractVector{Int},
    vecInd::AbstractVector{Int}, vals::AbstractVector,
    ::ScanSeg{true}, vis_ind, site, time, schema::GainSchema
    )

    ifirst = findfirst(site, schema.sites)
    t0 = scheme.times[ifirst]
    if t > t0
        ind = findall(x -> (((x[1]<=time)) && (x[2]==site)), schema.gts)
    else
        ind = findall(x -> (((x[1]==time)) && (x[2]==site)), schema.gts)
    end

    append!(colInd, ind)
    append!(rowInd, fill(i, length(ind)))
end

function fill_designmat!(
    colInd::AbstractVector{Int}, rowInd::AbstractVector{Int},
    vecInd::AbstractVector{Int}, vals::AbstractVector,
    f::FixedSeg, vis_ind, site, time, schema::GainSchema
    )
    push!(vecInd, vis_ind)
    push!(vals, f.value)
end

# Helper function that extracts the expected types for the FixedSeg segmentation
function fixed_type(segmenatation::NamedTuple)
    vs = values(segmenatation)
    inds = findall(x->(x isa FixedSeg), vs)
    return promote_type(map(x->typeof(getproperty(x, :value)), vs[inds])...)
end


function jonescache(obs::EHTObservation, segmentation::NamedTuple)
    @argcheck sort(collect(stations(obs))) == sort(collect(keys(segmentation)))

    times = obs[:T]
    bls = obs[:baseline]
    schema = gain_schema(segmentation, obs)
    rowInd1 = Int[]
    colInd1 = Int[]
    rowInd2 = Int[]
    colInd2 = Int[]

    vecInd1 = Int[]
    vecInd2 = Int[]
    T = fixed_type(segmentation)
    vals1 = T[]
    vals2 = T[]


    for i in eachindex(times)
        t = times[i]
        s1, s2 = bls[i]
        seg1 = getproperty(segmentation, s1)
        seg2 = getproperty(segmentation, s2)
        fill_designmat!(colInd1, rowInd1, vecInd1, vals1, seg1, i, s1, t, schema)
        fill_designmat!(colInd2, rowInd2, vecInd2, vals2, seg2, i, s2, t, schema)
    end

    z1 = fill(1.0, length(rowInd1))
    m1 = sparse(rowInd1, colInd1, z1, length(times), length(schema.times))
    z2 = fill(1.0, length(rowInd2))
    m2 = sparse(rowInd2, colInd2, z2, length(times), length(schema.times))

    if length(vecInd1) > 0
        v1 = sparsevec(vecInd1, vals1, length(times))
        d1 = AffineDesignMatrix(m1, v1)
    else
        d1 = m1
    end

    if length(vecInd2) > 0
        v2 = sparsevec(vecInd2, vals2, length(times))
        d2 = AffineDesignMatrix(m2, v2)
    else
        d2 = m2
    end

    return JonesCache{
                typeof(d1), typeof(d2), typeof(segmentation),
                typeof(schema)
                }(d1, d2, segmentation, schema)
end



"""
    jonescache(obs::EHTObservation, segmentation::ObsSegmentation)
    jonescache(obs::EHTObservatoin, segmentation::NamedTuple)

Constructs a `JonesCache` from a given observation `obs` using the segmentation scheme
`segmentation`. If `segmentation` is a named tuple it is assumed that each symbol in the
named tuple corresponds to a segmentation for thes sites in `obs`.

# Example
```julia-repl
# coh is a EHTObservation
julia> jonescache(coh, ScanSeg())
julia> segs = (AA = ScanSeg(), AP = TrachSeg(), AZ=FixedSegSeg())
julia> jonescache(coh, segs)
```
"""
function jonescache(obs::EHTObservation, s::ObsSegmentation)
    sites = Tuple(stations(obs))
    segs = NamedTuple{sites}(Tuple(s for _ in 1:length(sites)))
    return jonescache(obs, segs)
end

# This is an internal function that computes the set of stations from a ScanTable
function gain_stations(st::ScanTable)
    gainstat = Symbol[]
    times = eltype(st.times)[]
    for i in 1:length(st)
        s = stations(st[i])
        append!(gainstat, s)
        append!(times, fill(st[i].time, length(s)))
    end
    return times, gainstat
end


"""
    $(TYPEDEF)

Holds the pairs of Jones matrices for the first and second station of a baseline.

# Fields
$(FIELDS)
"""
struct JonesPairs{T, M1<:AbstractVector{T}, M2<:AbstractVector{T}}
    """
    Vector of jones matrices for station 1
    """
    m1::M1
    """
    Vector of jones matrices for station 2
    """
    m2::M2
end

function Base.:*(x::JonesPairs, y::JonesPairs...)
    m1 = map(x->getproperty(x, :m1), (x,y...))
    m2 = map(x->getproperty(x, :m2), (x,y...))
    o1, o2 = _allmul(m1, m2)
    JonesPairs(o1, o2)
end

function _allmul(m1, m2)
    out1 = zero(first(m1))
    out2 = zero(first(m2))
    _allmul!(out1, out2, m1, m2)
    # out1 = reduce(.*, m1)
    # out2 = reduce(.*, m2)
    return out1, out2
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
        Δm1 = zero(out[1])
        Δm1 .= unthunk(Δ[1])
        Δm2 = zero(out[2])
        Δm2 .= unthunk(Δ[2])
        dm1 = zero.(m1)
        dm2 = zero.(m2)

        out1 = zero(first(m1))
        out2 = zero(first(m2))
        autodiff(Reverse, _allmul!, Duplicated(out1, Δm1), Duplicated(out2, Δm2), Duplicated(m1, dm1), Duplicated(m2, dm2))
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


function apply_design(gmat::T, jcache::AbstractJonesCache) where {T}
    return JonesPairs(jcache.m1*gmat, jcache.m2*gmat)
end

function apply_designmats(f::F, g1, g2, m) where {F}
    return f.(m*g1), f.(m*g2)
end

function apply_designmats(::typeof(identity), g1, g2, m)
    return m*g1, m*g2
end

# GMat(g1::T, g2::T) where {T} = SMatrix{2,2,T}(g1, zero(T), zero(T), g2)
function gmat(f::F, g1, g2, m) where {F}
   S = eltype(g1)
   gs1, gs2 = apply_designmats(f, g1, g2, m)
   n = length(gs1)
   offdiag = fill(zero(S), n)
   StructArray{SMatrix{2,2,S,4}}((gs1, offdiag, offdiag, gs2))
end

"""
    jonesG(g1::AbstractVector, g2::AbstractVector, jcache::AbstractJonesCache)
    jonesG(f, g1::AbstractVector, g2::AbstractVector, jcache::AbstractJonesCache)

Constructs the pairs Jones `G` matrices for each pair of stations. The `g1` are the
gains for the first polarization basis and `g2` are the gains for the other polarization.
The first argument is optional and denotes a function that is applied to every element of
jones cache. For instance if `g1` and `g2` are the log-gains then `f=exp` will convert them
into the gains.

The layout for each matrix is as follows:
```
    g1 0
    0  g2
```
"""
function jonesG(f::F, g1::AbstractVector, g2::AbstractVector, jcache::AbstractJonesCache) where {F}
    gm1 = gmat(f, g1, g2, jcache.m1)
    gm2 = gmat(f, g1, g2, jcache.m2)
    return JonesPairs(gm1, gm2)
end
jonesG(g1::AbstractVector, g2::AbstractVector, jcache::AbstractJonesCache) = jonesG(identity, g1, g2, jcache)

@inline function jonesG_prod_ratio(gproduct, product_cache, gratio, ratio_cache)
    Gp = jonesG(gproduct, gproduct, product_cache)
    Gr = jonesG(gratio, inv.(gratio), ratio_cache)
    return Gp*Gr
end

# DMat(d1::T, d2::T) where {T} = SMatrix{2,2,T}(one(T), d2, d1, one(T))
function dmat(f::F, d1, d2, m) where {F}
    S = promote_type(eltype(d1), eltype(d2))
    ds1, ds2 = apply_designmats(f, d1, d2, m)
    n = length(ds1)
    unit = fill(one(S), n)
    return StructArray{SMatrix{2,2,S,4}}((unit, ds2, ds1, unit))
end

"""
    jonesD(d1::AbstractVector, d2::AbstractVector, jcache::AbstractJonesCache)
    jonesD(f, d1::AbstractVector, d2::AbstractVector, jcache::AbstractJonesCache)

Constructs the pairs Jones `D` matrices for each pair of stations. The `d1` are the
d-termsfor the first polarization basis and `d2` are the d-terms for the other polarization.
The first argument is optional and denotes a function that is applied to every element of
jones cache. For instance if `d1` and `d2` are the log-dterms then `f=exp` will convert them
into the dterms.

The layout for each matrix is as follows:
```
    1  d1
    d2 1
```
"""
function jonesD(f::F, d1::T,d2::T,jcache::AbstractJonesCache) where {F, T}
    dm1 = dmat(f, d1, d2, jcache.m1)
    dm2 = dmat(f, d1, d2, jcache.m2)
    return JonesPairs(dm1, dm2)
end
jonesD(d1::AbstractVector, d2::AbstractVector, jcache::AbstractJonesCache) = jonesD(identity, d1, d2, jcache)

export jonesStokes
"""
    jonesStokes(g1::AbstractArray, gcache::AbstractJonesCache)
    jonesStokes(f, g1::AbstractArray, gcache::AbstractJonesCache)

Construct the Jones Pairs for the stokes I image only. That is, we only need to
pass a single vector corresponding to the gain for the stokes I visibility. This is
for when you only want to image Stokes I.
The first argument is optional and denotes a function that is applied to every element of
jones cache. For instance if `g1` and `g2` are the log-gains then `f=exp` will convert them
into the gains.


# Warning
In the future this functionality may be removed when stokes I fitting is replaced with the
more correct `trace(coherency)`, i.e. RR+LL for a circular basis.
"""
function jonesStokes(f::F, g::AbstractVector, gcache::AbstractJonesCache) where {F}
    out1 = similar(g, size(gcache.m1, 1))
    out2 = similar(g, size(gcache.m1, 1))
    _jonesStokes!(out1, out2, f, g, gcache.m1, gcache.m2)
    return JonesPairs(out1, out2)
end
jonesStokes(g::AbstractVector, gcache::AbstractJonesCache) = jonesStokes(identity, g, gcache)

function jonesStokes(::typeof(identity), g::AbstractVector, gcache::AbstractJonesCache)
    out1 = similar(g, size(gcache.m1, 1))
    out2 = similar(g, size(gcache.m1, 1))
    _jonesStokes!(out1, out2, identity, g, gcache.m1, gcache.m2)
    return JonesPairs(out1, out2)
end

function _jonesStokes!(out1, out2, f::F, g::AbstractVector, m1, m2) where {F}
    mul!(out1, m1, g)
    mul!(out2, m2, g)
    out1 .= f.(out1)
    out2 .= f.(out2)
    return nothing
end

function _jonesStokes!(out1, out2, ::typeof(identity), g::AbstractVector, m1, m2)
    mul!(out1, m1, g)
    mul!(out2, m2, g)
    return nothing
end

function ChainRulesCore.rrule(::typeof(jonesStokes), f::F, g, gcache) where {F}
    j = jonesStokes(f, g, gcache)
    pg = ProjectTo(g)
    function _jonesStokes_pullback(Δ)
        dout1 = zero(j.m1)
        dout1 .= unthunk(Δ.m1)
        dout2 = zero(j.m2)
        dout2 .= unthunk(Δ.m2)

        dg = zero(g)

        out1 = zero(j.m1)
        out2 = zero(j.m2)

        autodiff(Reverse, _jonesStokes!, Const, Duplicated(out1, dout1), Duplicated(out2, dout2), Const(f), Duplicated(g, dg), Duplicated(gcache.m1, copy(gcache.m1)), Duplicated(gcache.m2, copy(gcache.m2)))
        return NoTangent(), NoTangent(), pg(dg), NoTangent()
    end
    return j, _jonesStokes_pullback
end

# function jonesStokes(f::F, g::AbstractVector, gcache::AbstractJonesCache) where {F}
#     return JonesPairs(f.(gcache.m1*g), f.(gcache.m2*g))
# end
# jonesStokes(g::AbstractVector, gcache::AbstractJonesCache) = jonesStokes(identity, g, gcache)
# function jonesStokes(::typeof(identity), g::AbstractVector, gcache::AbstractJonesCache)
#     return JonesPairs(gcache.m1*g, gcache.m2*g)
# end




"""
    $(TYPEDEF)

Holds various transformations that move from the measured telescope basis to the **chosen**
on sky reference basis.

# Fields
$(FIELDS)
"""
struct TransformCache{M, B<:PolBasis} <: AbstractJonesCache
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

"""
    TransformCache(obs::EHTObservation; add_fr=true, ehtim_fr_convention=true, ref::PolBasis=CirBasis())

Constructs the cache that holds the transformation from the **chosen** on-sky reference basis
to the basis that the telescope measures the electric fields given an observation `obs`.
Our convention is that the feed rotations
are included in the transformation cache with the `add_fr` toggle. The user can specify the reference basis
using the keyword argument `ref` which is the (R,L) circular basis by default.

# Notes
We use the following definition for our feed rotations

```
    exp(-iθ)  0
    0         exp(iθ)
```

# Warning
eht-imaging can sometimes pre-rotate the coherency matrices. As a result the field rotation can sometimes
be applied twice. To compensate for this we have added a `ehtim_fr_convention` which will fix this.
"""
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

"""
    jonesT(tcache::TransformCache)

Returns a `JonesPair` of matrices that transform from the model coherency matrices basis
to the on-sky coherency basis, this includes the feed rotation and choice of polarization feeds.
"""
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



"""
    JonesModel(jones::JonesPairs, model, tcache::TransformCache)
    JonesModel(jones::JonesPairs, model, refbasis::PolBasis=CirBasis())

Constructs a `JonesModel` from a `jones` pairs that describe the intrument model
and the `model` which describes the on-sky polarized visibilities. The third argument
can either be the `tcache` that converts from the model coherency basis to the instrumental
basis, or just the `refbasis` that will be used when constructing the model coherency matrices.
"""
function JonesModel(jones, model, tcache::TransformCache)
    return JonesModel(jones, model, tcache.refbasis)
end

function JonesModel(jones, model)
    return JonesModel(jones, model, CirBasis())
end

function visibilities_analytic(model::JonesModel{J,M,B}, u, v, time, freq) where {J,M,B}
    vis = visibilities_analytic(model.model, u, v, time, freq)
    coh = _coherency(vis, B)
    return corrupt(coh, model.jones.m1, model.jones.m2)
end

function visibilities_numeric(model::JonesModel{J,M,B}, u, v, time, freq) where {J,M,B}
    vis = visibilities_numeric(model.model, u, v, time, freq)
    coh = _coherency(vis, B)
    return corrupt(coh, model.jones.m1, model.jones.m2)
end


"""
    corrupt(vis, j1, j2)

Corrupts the model coherency matrices with the Jones matrices `j1` for station 1 and
`j2` for station 2.
"""
function corrupt(vis, j1, j2)
    vnew = j1 .* vis .* adjoint.(j2)
    return vnew
end


# function corrupt!(vnew, vis::AbstractArray, j1, j2)
#     # @assert length(vis) == length(j1) "visibility vector and jones pairs have mismatched dimensions!"
#     vnew .= j1 .* vis .* adjoint.(j2)
#     return nothing
# end

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

"""
    extract_FRs

Extracts the feed rotation Jones matrices (returned as a `JonesPair`) from an EHT
observation `obs`.

# Warning
eht-imaging can sometimes pre-rotate the coherency matrices. As a result the field rotation can sometimes
be applied twice. To compensate for this we have added a `ehtim_fr_convention` which will fix this.

"""
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
