"""
    shard(post::VLBIPosterior, layout::ShardLayout)
    shard(post::VLBIPosterior, sharding::Reactant.Sharding.AbstractSharding)

Return a copy of the host posterior `post` with its arrays placed on Reactant devices.

With a `layout`, each field's dual domain is split by `layout` (see `shard` for dual
domains): when it splits image planes such as `Ti`, every device holds the visibilities of
its own planes, padded to a common length. The data, weights and the instrument's
per-visibility station-cell indices follow the visibility points, with padding at weight 0.
Everything else is replicated: the sky function and its metadata, the instrument terms and
their per-cell values, and the image-domain likelihoods. Dim lookups stay on the host. Every
dimension named in `layout` must be a dimension of each field's image or visibility domain.

With a `sharding`, every array is placed with it: `shard(post, Reactant.Sharding.NoSharding())`
places the posterior on one device.

The prior stays on the host, as constants of the compiled program.

The model visibilities of a posterior sharded by a layout are in the blocked order of its
visibility domain; `post.data` keeps the host tables. Host utilities such as
[`residuals`](@ref), [`chi2`](@ref) and [`simulate_observation`](@ref) take the host
posterior.

`post` must be on the host: an array already on a device throws, since placing it again
would first hold it whole on one device. Closure products and self-calibration are not
supported with a layout.
"""
ComradeBase.shard(post::VLBIPosterior, layout::ShardLayout) = _shardposterior(post, layout)
ComradeBase.shard(post::VLBIPosterior, sharding::Reactant.Sharding.AbstractSharding) = _shardposterior(post, sharding)

function _shardposterior(post::VLBIPosterior, p)
    _ondevice(post) && throw(
        ArgumentError("the posterior holds arrays already on a device; shard the posterior built on the host")
    )
    sky = _placesky(post.skymodel, p)
    pts = _points(post.skymodel, sky, p)
    return Accessors.setproperties(
        post, (;
            skymodel = sky,
            instrumentmodel = _placeinstrument(post.instrumentmodel, pts, p),
            lklhds = _placelikelihoods(post.lklhds, post.data, pts),
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

_replicate(x, layout::ShardLayout) = Reactant.to_rarray(x; sharding = Reactant.Sharding.Replicated(layout.mesh))
_replicate(x, sharding) = Reactant.to_rarray(x; sharding)

# Where the per-value arrays of one table go: `host` is the table's domain of records and
# channels, `to` the sharded visibility domain of its sky, or one sharding for every array.
struct _Points{H, T}
    host::H
    to::T
end

# The per-value map `m` over `p.host`, on the devices; slots added for padding hold `pad`.
_pervalue(m::IntensityMap, p::_Points{<:Any, <:StructuredDomain}, pad) = shard(m, p.to; pad)
_pervalue(m::IntensityMap, p::_Points, pad) = IntensityMap(Reactant.to_rarray(baseimage(m); sharding = p.to), axisdims(m))

## Sky

_placesky(m::Comrade.ObservedSkyModel, p) =
    Comrade.ObservedSkyModel(_replicate(m.f, p), _placedual(m.grid, p), _replicate(m.metadata, p))

# The precomputed visibilities of a fixed model follow the visibility points.
function _placesky(m::Comrade.ObservedSkyModel{<:Comrade.FixedSkyModel}, p)
    grid = _placedual(m.grid, p)
    pts = _points(m, Comrade.ObservedSkyModel(m.f, grid, m.metadata), p)
    vis = _pervalue(m.metadata.vis, pts, zero(eltype(m.metadata.vis)))
    return Comrade.ObservedSkyModel(_replicate(m.f, p), grid, (; model = _replicate(m.metadata.model, p), vis))
end

_placesky(m::Comrade.MultiSkyModel, p) = Comrade.MultiSkyModel(map(sm -> _placesky(sm, p), m.skymodels))
_placesky(m::Comrade.FieldSkies, p) = Comrade.FieldSkies(map(sm -> _placesky(sm, p), m.fields))

_devicealg(alg::VLBISkyModels.ReactantNUFFTAlg, T) = alg
_devicealg(alg::Comrade.AnalyticAlg, T) = alg
_devicealg(alg, T) = VLBISkyModels.ReactantNUFFTAlg(T)

function _placedual(g::VLBISkyModels.FourierDualDomain, layout::ShardLayout)
    img = ComradeBase.imgdomain(g)
    alg = _devicealg(VLBISkyModels.algorithm(g), eltype(img))
    return shard(VLBISkyModels.FourierDualDomain(img, ComradeBase.visdomain(g), alg), layout)
end

function _placedual(g::VLBISkyModels.FourierDualDomain, sharding)
    img = ComradeBase.imgdomain(g)
    visr = Reactant.to_rarray(ComradeBase.visdomain(g); sharding)
    return VLBISkyModels.FourierDualDomain(img, visr, _devicealg(VLBISkyModels.algorithm(g), eltype(img)))
end

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

_placeinstrument(J::Comrade.ObservedIdealInstrument, pts, p) = J
_placeinstrument(J::Comrade.FieldInstruments, pts::NamedTuple, p) =
    Comrade.FieldInstruments(map((j, q) -> _placeinstrument(j, q, p), J.fields, pts))

function _placeinstrument(J::Comrade.ObservedInstrumentModel, pts, p)
    return Comrade.ObservedInstrumentModel(
        _placeterms(J.instrument, p), J.refbasis, _placecells(J.layout, pts, p), J.pol, Comrade.corrected(J)
    )
end

_placeterms(J, p) = _replicate(J, p)
_placeterms(J::Comrade._StationGroups, p) = Comrade._StationGroups(_replicate(J.terms, p), J.cells)

# Each value's station cells follow the visibility points (padding reads cell 1); the per-cell
# values are replicated.
function _placecells(L::Comrade.CellLayout, pts, p)
    cell(c) = baseimage(_pervalue(IntensityMap(c, pts.host), pts, one(eltype(c))))
    coords = StructArray(map(c -> _replicate(c, p), StructArrays.components(L.coords)))
    response = isnothing(L.response) ? nothing : _replicate(L.response, p)
    return Comrade.CellLayout(coords, cell(L.cell1), cell(L.cell2), _replicate(L.entries, p), response)
end

function _placeinstrument(J::Comrade.ObservedSelfCalInstrumentModel, pts, ::ShardLayout)
    throw(ArgumentError("self-calibration does not support a ShardLayout; shard with a single sharding"))
end

function _placeinstrument(J::Comrade.ObservedSelfCalInstrumentModel, pts, sharding)
    stacks = map(st -> merge(st, (; eye = _replicate(st.eye, sharding))), J.layout.stacks)
    return Comrade.ObservedSelfCalInstrumentModel(
        _placeinstrument(J.model, pts, sharding), J.phase,
        map(w -> _replicate(w, sharding), J.wdata), map(w -> _replicate(w, sharding), J.weight),
        merge(J.layout, (; stacks)), J.squarings, J.keep, J.entry1, J.entry2
    )
end

## Likelihoods

_placelikelihoods(ls::NamedTuple, data::NamedTuple, pts::NamedTuple) = map(_placelikelihoods, ls, data, pts)
_placelikelihoods(ls::Tuple, data::Tuple, pts) = map((l, d) -> _placelikelihood(l, d, pts), ls, data)

function _placelikelihood(l::Comrade.ConditionedLikelihood, t::VisTable, pts)
    meas, w = measurement(t), weight(t)
    kernel = Accessors.setproperties(l.kernel, (; w = Comrade._dense(_pervalue(w, pts, zero(eltype(w))))))
    return Comrade.ConditionedLikelihood(kernel, Comrade._dense(_pervalue(meas, pts, zero(eltype(meas)))))
end

function _placelikelihood(l, ::Comrade.ClosureTable, pts::_Points{<:Any, <:StructuredDomain})
    throw(ArgumentError("closure products do not support a ShardLayout; shard with a single sharding"))
end
_placelikelihood(l, ::Comrade.ClosureTable, pts::_Points) = _replicate(l, pts.to)
