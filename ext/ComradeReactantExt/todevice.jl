"""
    todevice(post::VLBIPosterior, layout::ShardLayout)
    todevice(post::VLBIPosterior, sharding::Reactant.Sharding.AbstractSharding)

Return a copy of the host posterior `post` with its arrays placed on Reactant devices.

With a `layout`, each field's dual domain is split by `layout` (see `todevice` for dual
domains): when it splits image planes such as `Ti`, every device holds the visibilities of
its own planes, padded to a common length. The data, weights and the instrument's
per-visibility station-cell indices follow the visibility points, with padding at weight 0.
Everything else is replicated: the sky function and its metadata, the instrument terms and
their per-cell values, and the image-domain likelihoods. Dim lookups stay on the host. Every
dimension named in `layout` must be a dimension of each field's image or visibility domain,
except in a `MultiSkyModel`: a component whose image lacks the split dims, such as a static
image next to a `Ti` cube, is evaluated on the visibility points of the first component that
has them, on the devices that hold those points.

Closure products and the self-calibrated station phase are replicated: every device reads
all model visibilities at the measured values, so XLA gathers them from the devices that
hold them.

With a `sharding`, every array is placed with it: `todevice(post, Reactant.Sharding.NoSharding())`
places the posterior on one device.

The prior stays on the host, as constants of the compiled program.

The model visibilities of a posterior sharded by a layout are in the blocked order of its
visibility domain; `post.data` keeps the host tables. Host utilities such as
[`residuals`](@ref), [`chi2`](@ref) and [`simulate_observation`](@ref) take the host
posterior.

`post` must be on the host: an array already on a device throws, since placing it again
would first hold it whole on one device.
"""
ComradeBase.todevice(post::VLBIPosterior, layout::ShardLayout) = _shardposterior(post, layout)
ComradeBase.todevice(post::VLBIPosterior, sharding::Reactant.Sharding.AbstractSharding) = _shardposterior(post, sharding)

function _shardposterior(post::VLBIPosterior, p)
    _ondevice(post) && throw(
        ArgumentError("the posterior holds arrays already on a device; place the posterior built on the host")
    )
    sky = _placesky(post.skymodel, p)
    pts = _points(post.skymodel, sky, p)
    return Accessors.setproperties(
        post, (;
            skymodel = sky,
            instrumentmodel = _placeinstrument(post.instrumentmodel, pts),
            lklhds = _placelikelihoods(post.lklhds, post.data, pts, p),
            lklhdsimg = _replicate(post.lklhdsimg, p),
        )
    )
end

# Whether `x` holds an array or number on a device.
_ondevice(x) = _ondevice(x, IdDict{Any, Nothing}())
function _ondevice(x, seen)
    x isa Union{Reactant.AbstractConcreteArray, Reactant.AbstractConcreteNumber} && return true
    T = typeof(x)
    (isbitstype(T) || x isa Union{Module, Type, Symbol, AbstractString}) && return false
    haskey(seen, x) && return false
    seen[x] = nothing
    if x isa AbstractArray
        isbitstype(eltype(x)) && return false
        return any(i -> isassigned(x, i) && _ondevice(x[i], seen), eachindex(x))
    end
    return any(i -> isdefined(x, i) && _ondevice(getfield(x, i), seen), 1:fieldcount(T))
end

_replicate(x, layout::ShardLayout) = todevice(x, Reactant.Sharding.Replicated(layout.mesh))
_replicate(x, sharding) = todevice(x, sharding)

# Where the per-value arrays of one table go: `host` is the table's domain of records and
# channels, `to` the sharded visibility domain of its sky, or one sharding for every array.
struct _Points{H, T}
    host::H
    to::T
end

# The per-value map `m` over `p.host`, on the devices; slots added for padding hold `pad`.
_pervalue(m::IntensityMap, p::_Points{<:Any, <:StructuredDomain}, pad) = todevice(m, p.to; pad)
_pervalue(m::IntensityMap, p::_Points, pad) = todevice(m, p.to)

## Sky

_placesky(m::Comrade.ObservedSkyModel, p) = _placesky(m, p, _placedual(m.grid, p))
_placesky(m::Comrade.ObservedSkyModel, p, grid) =
    Comrade.ObservedSkyModel(_replicate(m.f, p), grid, _replicate(m.metadata, p))

# The precomputed visibilities of a fixed model follow the visibility points.
function _placesky(m::Comrade.ObservedSkyModel{<:Comrade.FixedSkyModel}, p, grid)
    pts = _points(m, Comrade.ObservedSkyModel(m.f, grid, m.metadata), p)
    vis = _pervalue(m.metadata.vis, pts, zero(eltype(m.metadata.vis)))
    return Comrade.ObservedSkyModel(_replicate(m.f, p), grid, (; model = _replicate(m.metadata.model, p), vis))
end

_placesky(m::Comrade.MultiSkyModel, p) = Comrade.MultiSkyModel(map(sm -> _placesky(sm, p), m.skymodels))
_placesky(m::Comrade.FieldSkies, p) = Comrade.FieldSkies(map(sm -> _placesky(sm, p), m.fields))

# A component whose image lacks the dims the layout splits goes on the visibility points of the
# first component that has them.
function _placesky(m::Comrade.MultiSkyModel, layout::ShardLayout)
    names = keys(m.skymodels)
    lead = findfirst(sm -> _haslayoutdims(sm.grid, layout), values(m.skymodels))
    isnothing(lead) && return Comrade.MultiSkyModel(map(sm -> _placesky(sm, layout), m.skymodels))
    placed = _placesky(m.skymodels[lead], layout)
    vis = ComradeBase.visdomain(placed.grid)
    skies = map(names, values(m.skymodels)) do k, sm
        k === names[lead] && return placed
        _haslayoutdims(sm.grid, layout) && return _placesky(sm, layout)
        return _placesky(sm, layout, _placedual(sm.grid, layout, vis))
    end
    return Comrade.MultiSkyModel(NamedTuple{names}(skies))
end

_haslayoutdims(g, layout::ShardLayout) =
    keys(layout.axes) ⊆ map(DD.name, (DD.dims(ComradeBase.imgdomain(g))..., DD.dims(ComradeBase.visdomain(g))...))

_devicealg(alg::VLBISkyModels.ReactantNUFFTAlg, T) = alg
_devicealg(alg::Comrade.AnalyticAlg, T) = alg
_devicealg(alg, T) = VLBISkyModels.ReactantNUFFTAlg(T)

function _placedual(g::VLBISkyModels.FourierDualDomain, layout::ShardLayout)
    img = ComradeBase.imgdomain(g)
    alg = _devicealg(VLBISkyModels.algorithm(g), eltype(img))
    return todevice(VLBISkyModels.FourierDualDomain(img, ComradeBase.visdomain(g), alg), layout)
end

function _placedual(g::VLBISkyModels.FourierDualDomain, sharding)
    img = ComradeBase.imgdomain(g)
    visr = todevice(ComradeBase.visdomain(g), sharding)
    return VLBISkyModels.FourierDualDomain(img, visr, _devicealg(VLBISkyModels.algorithm(g), eltype(img)))
end

# The dual domain of `g` on the blocked visibility points `vis` of another component, split
# over the devices that hold those points.
function _placedual(g::VLBISkyModels.FourierDualDomain, ::ShardLayout, vis::StructuredDomain)
    img = ComradeBase.imgdomain(g)
    host = DD.rebuild(
        vis; coords = map(_tohost, ComradeBase.coords(vis)), executor = ComradeBase.executor(ComradeBase.visdomain(g))
    )
    alg = _devicealg(VLBISkyModels.algorithm(g), eltype(img))
    blocked = ComradeBase.executor(vis).layout
    return todevice(VLBISkyModels.FourierDualDomain(img, host, alg), ShardLayout(blocked.mesh; Node = blocked.axes.Node))
end

_tohost(c::Reactant.AbstractConcreteArray) = Array(c)
_tohost(c) = c

# The points of each field, from its host and placed sky models.
_points(host::Comrade.ObservedSkyModel, placed::Comrade.ObservedSkyModel, ::ShardLayout) =
    _Points(ComradeBase.visdomain(host.grid), ComradeBase.visdomain(placed.grid))
_points(host::Comrade.ObservedSkyModel, placed::Comrade.ObservedSkyModel, sharding) =
    _Points(ComradeBase.visdomain(host.grid), sharding)
_points(host::Comrade.FieldSkies, placed::Comrade.FieldSkies, p) = map((h, s) -> _points(h, s, p), host.fields, placed.fields)

# The components of a multi-component sky add their visibilities, so they must share points.
function _points(host::Comrade.MultiSkyModel, placed::Comrade.MultiSkyModel, p)
    pts = map((h, s) -> _points(h, s, p), values(host.skymodels), values(placed.skymodels))
    nodes = map(q -> _nodes(q.to), pts)
    allequal(nodes) || throw(
        ArgumentError("the components of a MultiSkyModel place their visibilities differently over the devices; give them image grids that the layout splits the same way")
    )
    return first(pts)
end

_nodes(d::StructuredDomain) = parent(DD.lookup(first(DD.dims(d))))
_nodes(_) = nothing

## Instrument

_placeinstrument(J, pts::_Points) = todevice(J, pts.to)
_placeinstrument(J, pts::NamedTuple) = todevice(J, map(q -> q.to, pts))

## Likelihoods

_placelikelihoods(ls::NamedTuple, data::NamedTuple, pts::NamedTuple, p) = map((l, d, q) -> _placelikelihoods(l, d, q, p), ls, data, pts)
_placelikelihoods(ls::Tuple, data::Tuple, pts, p) = map((l, d) -> _placelikelihood(l, d, pts, p), ls, data)

function _placelikelihood(l::Comrade.ConditionedLikelihood, t::VisTable, pts, p)
    meas, w = measurement(t), weight(t)
    kernel = Accessors.setproperties(l.kernel, (; w = Comrade._dense(_pervalue(w, pts, zero(eltype(w))))))
    return Comrade.ConditionedLikelihood(kernel, Comrade._dense(_pervalue(meas, pts, zero(eltype(meas)))))
end

# The linear index in the model visibilities on the devices of each value of the host domain.
_slots(pts::_Points{<:Any, <:StructuredDomain}) = ComradeBase._slots(DD.dims(pts.host), pts.to)
_slots(pts::_Points) = LinearIndices(size(pts.host))

# Closures replicated on every device, by blocks of their covariance: the closures of one time
# stamp and channel. `cols[l, i, k]` and `vals[l, i, k]` are the visibility slots and
# coefficients of the `i`-th closure of block `k`, and `W[:, :, k]` whitens the block,
# `W Σ Wᵀ = I`. Padding has zero coefficients and zero rows and columns of `W`, and reads a
# measured visibility so that its log amplitude stays finite.
struct _BlockClosures{K, I, V, M, L}
    product::K
    cols::I
    vals::V
    W::M
    L::L
end

struct _ClosureMean{B, C}
    blocks::B
    μ::C
end

function (b::_BlockClosures)(vis)
    x = Comrade._closurevalue.(Ref(b.product), vec(baseimage(vis))[b.cols])
    return _ClosureMean(b, dropdims(sum(b.vals .* x; dims = 1); dims = 1))
end

function DensityInterface.logdensityof(d::_ClosureMean, y)
    r = _closureresidual(d.blocks.product, y, d.μ)
    z = sum(d.blocks.W .* reshape(r, 1, size(r)...); dims = 2)
    return d.blocks.L - sum(abs2, z) / 2
end

_closureresidual(::Comrade.ClosurePhases, y, μ) = cis.(y) .- cis.(μ)
_closureresidual(::Comrade.LogClosureAmplitudes, y, μ) = y .- μ

function _placelikelihood(l::Comrade.ConditionedLikelihood, ct::Comrade.ClosureTable, pts, p)
    slots = _slots(pts)
    keep = Comrade._measuredvalues(getfield(ct, :vis))
    cols, vals, W, y = _closureblocks(ct, map(I -> slots[I], keep))
    blocks = _BlockClosures(ct.product, _replicate(cols, p), _replicate(vals, p), _replicate(W, p), l.kernel.L)
    return Comrade.ConditionedLikelihood(blocks, _replicate(y, p))
end

# The padded arrays of `_BlockClosures` and the measured closures by block, for the closure
# table `ct` whose `j`-th visibility is at slot `slot[j]`.
function _closureblocks(ct::Comrade.ClosureTable, slot)
    D, Σ = getfield(ct, :design), getfield(ct, :covariance)
    d = Comrade.domain(ct)
    stamp = collect(zip(d.ts, d.Fr))
    blocks = [findall(==(s), stamp) for s in unique(stamp)]
    sum(r -> nnz(Σ[r, r]), blocks) == nnz(Σ) ||
        error("the closure covariance is not block diagonal over time stamps and channels")
    Dt = sparse(D')
    T = float(eltype(D))
    nl = maximum(r -> length(nzrange(Dt, r)), axes(Dt, 2))
    nc = maximum(length, blocks)
    cols = fill(first(slot), nl, nc, length(blocks))
    vals = zeros(T, nl, nc, length(blocks))
    W = zeros(T, nc, nc, length(blocks))
    y = zeros(T, nc, length(blocks))
    x = Vector(baseimage(measurement(ct)))
    for (k, rows) in pairs(blocks)
        for (i, r) in pairs(rows), (l, q) in pairs(nzrange(Dt, r))
            cols[l, i, k] = slot[rowvals(Dt)[q]]
            vals[l, i, k] = nonzeros(Dt)[q]
        end
        n = length(rows)
        W[1:n, 1:n, k] = inv(cholesky(Symmetric(Matrix(Σ[rows, rows]))).L)
        y[1:n, k] = x[rows]
    end
    return cols, vals, W, y
end
