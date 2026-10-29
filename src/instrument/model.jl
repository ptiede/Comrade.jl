export InstrumentModel

"""
    $(TYPEDEF)

The abstract instrument model. For a concrete implementation see [`IdealInstrumentModel`](@ref) and
[`InstrumentModel`](@ref).

Any subtype must implement the following methods

 - `set_array(m::AbstractInstrumentModel, t::VisTable)`: binds `m` to the visibility table
    `t`, on whose measured values the priors are built, and returns the observed instrument
    model and its prior.
 - `apply_instrument(vis, m, x)`: applies the observed model `m` to the sky's visibilities
    `vis` at the parameters `x`, returning values over the table's `(Node, Fr)` records.
"""
abstract type AbstractInstrumentModel end

"""
    IdealInstrumentModel()

An instrument with no corruptions and no response: the sky's visibilities, in the circular
basis, are the data.
"""
struct IdealInstrumentModel <: AbstractInstrumentModel end

Base.show(io::IO, mime::MIME"text/plain", m::IdealInstrumentModel) = printstyled(io, "IdealInstrumentModel"; color = :light_cyan, bold = true)

@inline apply_instrument(vis, ::IdealInstrumentModel, x) = vis

struct InstrumentModel{J <: AbstractJonesMatrix, PI, P <: PolBasis} <: AbstractInstrumentModel
    jones::J
    prior::PI
    refbasis::P
    gaugefix::Union{ReportGauge, PinGauge}
    selfcal::Bool
    feedrotation_corrected::Bool
    function InstrumentModel{J, PI, P}(jones, prior, refbasis, gaugefix, selfcal, feedrotation_corrected) where {J, PI, P}
        gaugefix isa Union{ReportGauge, PinGauge} ||
            throw(ArgumentError("InstrumentModel gaugefix must be ReportGauge() or PinGauge(), got $(repr(gaugefix))"))
        selfcal && haskey(prior, SELFCAL_RESIDUAL) && throw(
            ArgumentError("a self-calibrated instrument model adds the parameter `$SELFCAL_RESIDUAL`; the prior already has one")
        )
        return new{J, PI, P}(jones, prior, refbasis, gaugefix, selfcal, feedrotation_corrected)
    end
end

function InstrumentModel(jones, prior, refbasis, gaugefix, selfcal, feedrotation_corrected)
    return InstrumentModel{typeof(jones), typeof(prior), typeof(refbasis)}(
        jones, prior, refbasis, gaugefix, selfcal, feedrotation_corrected
    )
end

function Base.show(io::IO, ::MIME"text/plain", m::InstrumentModel)
    printstyled(io, "InstrumentModel"; bold = true, color = :light_cyan)
    println(io)
    T = typeof(m.jones)
    ST = split(split(" $T", '{')[1], ".")[end]
    println(io, "  with Jones: ", ST)
    m.selfcal && println(io, "  with the station phase self-calibrated")
    m.feedrotation_corrected && println(io, "  for data corrected for feed rotation")
    return print(io, "  with reference basis: ", m.refbasis)
end


"""
    InstrumentModel(jones; refbasis = CirBasis(), feedrotation_corrected = false,
                    gaugefix = ReportGauge(), selfcal = false)

Builds an instrument model from the instrument term `jones`, whose terms own the priors of
their parameters (see [`Jones`](@ref)). The measured
correlations of the baseline `ab` are

    r_ab = J_a U_a X_ab U_bᴴ J_bᴴ

where `X_ab` are the sky's ideal visibilities in the global basis `refbasis`, `U_a` is the
unitary matrix from the global basis to station `a`'s local ideal feed basis, and `J_a` is
`jones` in that local basis. `U_a` is geometric and computed for every station, time and channel
from the station's `feeds`, its `feed_angles` and the feed rotation of its `mount` (see
[`ArrayInfo`](@ref)). Data already corrected for feed rotation, each station's
`U_a⁻¹` applied to its feeds, are modeled by `U_aᴴ J_a U_a` instead of `J_a U_a`: set
`feedrotation_corrected = true`. Stokes I data use no `U`.

# Arguments

 - `jones` : the instrument: a [`Jones`](@ref) term, a product of them, or a
    [`PerStation`](@ref) assignment of terms to stations.


# Optional Arguments
  - `refbasis`: the global basis of the sky's electric field, `CirBasis()` (the IAU circular
    basis, the default) or `LinBasis()`.
  - `feedrotation_corrected`: whether the data were corrected for feed rotation (default
    `false`).
  - `gaugefix`: What to do when the priors declared with `gauge = StationPhase()` (see
    [`ArrayPrior`](@ref)) add up to a station phase the data cannot determine.
    [`ReportGauge()`](@ref ReportGauge) (the default) reports the flat directions and the
    entries that would fix them; [`PinGauge()`](@ref PinGauge) pins those entries itself, adding them to
    each affected prior's referencing scheme and logging one line per pin. The analysis is
    [`Comrade.gauge_pins`](@ref) and runs when the model is bound to a visibility table.
  - `selfcal`: whether each station's common phase `φ` is written as its self-calibration
    solution `φ̂` plus a sampled residual, `φ = φ̂ + selfcal_residual`. `jones` must then omit
    the station phase. See the section below.

# Self-calibrated station phase

With `selfcal = true`, the model adds the instrument parameter `selfcal_residual`: one angle
per station and data time, with a prior uniform on the circle. In each group of stations that
baselines link at one time and frequency, the station with the largest total data weight
`Σ 1/σ²` is the reference, and its residual is fixed to 0. The priors of all other parameters
are those their terms give.

`φ̂` is the phase of the leading eigenvector of each group's Hermitian station matrix
```
A_ij = Σ_p V_ij,p conj(M_ij,p) / σ²_ij,p      (i ≠ j)
A_ii = Σ_j Σ_p |M_ij,p|² / σ²_ij,p
```
where `V` and `σ` are the data and their noise, `M` the model visibilities without the
station phase, and `p` runs over the correlation products. The off-diagonal part is the
phase-dependent term of the Gaussian log-likelihood; the diagonal makes `A` positive
semidefinite when the model matches the data. The eigenvector is approximated by the column
of `A^64`, computed by six squarings, at the reference station, so `φ̂` is zero there.

`φ̂` depends on the sky and the other instrument parameters but not on the residual, so the
change of variables from `φ` to the residual has unit Jacobian. Because the prior is uniform
on the circle, the prior of `φ` does not depend on `φ̂`, and the posterior is that of a uniform
station phase. The residual is a coordinate defined with the data, not a model parameter:
compare or reuse draws through `φ` ([`selfcal_phase`](@ref)). The model needs the data, so it
is bound to data only through a [`VLBIPosterior`](@ref) of one table of complex
visibilities or coherencies.


# Example

A Stokes I example is
```julia-repl
julia> G = Jones(
           (x, c) -> exp(complex(x.lg, x.pg));
           lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
           pg = ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))); refant = SensitiveReference(0.0)),
       )

julia> intm = InstrumentModel(G)
```

A standard polarized example is
```julia-repl
julia> G = Jones(;
           lgr = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
           gpr = ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))); refant = SensitiveReference(0.0)),
           lgrat = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
           gprat = ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))); refant = SensitiveReference(0.0)),
       ) do x, c
           gR = exp(complex(x.lgr, x.gpr))
           gL = gR * exp(complex(x.lgrat, x.gprat))
           return gains(gR, gL)
       end
julia> D = Jones(
           (x, c) -> leakage(complex(x.dRre, x.dRim), complex(x.dLre, x.dLim));
           dRre = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.1))),
           dRim = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.1))),
           dLre = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.1))),
           dLim = ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, 0.1))),
       )
julia> intm = InstrumentModel(G * D)
```
which builds the gains from the right-hand gain and the gain ratios, and `D` is the small
leakage matrix. [`@instrument`](@ref) writes the same models more compactly. An instrument
with no terms, only the response, is `InstrumentModel()`.
"""
function InstrumentModel(
        jones::AbstractJonesMatrix;
        refbasis = CirBasis(), feedrotation_corrected = false, gaugefix = ReportGauge(), selfcal = false
    )
    return InstrumentModel(jones, _priors(jones), refbasis, gaugefix, selfcal, feedrotation_corrected)
end

# The instrument with no terms of its own: only the response `U`.
struct IdentityJones <: AbstractJonesMatrix end
@inline jonesmatrix(::IdentityJones, x, c) = SMatrix{2, 2, ComplexF64, 4}(1, 0, 0, 1)

InstrumentModel(; refbasis = CirBasis(), feedrotation_corrected = false) =
    InstrumentModel(IdentityJones(), NamedTuple(), refbasis, ReportGauge(), false, feedrotation_corrected)

# `C`: whether the data were corrected for feed rotation.
struct ObservedInstrumentModel{I <: AbstractJonesMatrix, PB <: PolBasis, L <: CellLayout, P, C} <: AbstractInstrumentModel
    instrument::I
    refbasis::PB
    layout::L
    pol::P
end

function ObservedInstrumentModel(instrument, refbasis, layout, pol, corrected::Bool)
    return ObservedInstrumentModel{typeof(instrument), typeof(refbasis), typeof(layout), typeof(pol), corrected}(
        instrument, refbasis, layout, pol
    )
end
corrected(::ObservedInstrumentModel{I, PB, L, P, C}) where {I, PB, L, P, C} = C

instrument(x::ObservedInstrumentModel) = x.instrument
refbasis(x::ObservedInstrumentModel) = x.refbasis
layout(x::ObservedInstrumentModel) = x.layout
EnzymeRules.inactive(::typeof(instrument), args...) = nothing
EnzymeRules.inactive(::typeof(refbasis), args...) = nothing
EnzymeRules.inactive(::typeof(layout), args...) = nothing
EnzymeRules.inactive_type(::Type{<:ObservedInstrumentModel}) = true

function Base.show(io::IO, mime::MIME"text/plain", m::ObservedInstrumentModel)
    printstyled(io, "ObservedInstrumentModel"; bold = true, color = :light_cyan)
    println(io)
    T = typeof(m.instrument)
    ST = split(split(" $T", '{')[1], ".")[end]
    println(io, "  with Jones: ", ST)
    println(io, "  with ", length(m.layout.coords), " station cells")
    return print(io, "  with reference basis: ", m.refbasis)
end

"""
    set_array(int, t::VisTable)

The observed instrument model of `int` for the visibility table `t` and its prior. The priors
are built on the measured values of `t`.
"""
function set_array(int::InstrumentModel, t::VisTable)
    int.selfcal && int.jones isa PerStation && throw(ArgumentError("self-calibration does not support PerStation models"))
    int.selfcal && return _set_selfcal(int, t)
    prior, x = _instrumentprior(int, t)
    return _observe(int, t, x), prior
end

"""
    set_array(int, tables::NamedTuple)

The observed instrument models of `int` for the visibility tables of several fields, keyed by
field name, and their one prior. The priors are built on the measured values of all tables,
on the time axis of the earliest epoch, so a parameter is shared by the fields wherever its
segments cover the data of several.
"""
function set_array(int::InstrumentModel, tables::NamedTuple)
    int.selfcal && throw(ArgumentError("self-calibration does not support several fields"))
    _, shift = _epochshift(values(tables))
    prior, x = _instrumentprior(int, _measuredtable(values(tables)))
    observed = map((t, s) -> _observe(int, t, x, s), values(tables), shift)
    return FieldInstruments(NamedTuple{keys(tables)}(observed)), prior
end

# The prior of the instrument parameters, built on the measured values of `t`, and a draw.
function _instrumentprior(int::InstrumentModel, t::VisTable)
    jones = int.jones
    _checkstations(jones, info(t).stations.name)
    s = _measuredpoints(t).sites
    prior = _stationpriors(jones, int.prior, sort(unique(vcat(first.(s), last.(s)))))
    obs = map(Base.Fix2(ObservedArrayPrior, t), prior)
    obs = fix_phase_gauge(prior, obs, t, int.gaugefix)
    prior_obs = NamedDist(obs)
    return prior_obs, rand(prior_obs)
end

# `int` bound to the table `t`, whose times are `shift` hours after the epoch of the site arrays
# of the parameters `x`.
function _observe(int::InstrumentModel, t::VisTable, x, shift = 0.0)
    (; jones, refbasis) = int
    layout = CellLayout(t, map(siteparams, x), refbasis, _readers(jones, int.prior, info(t).stations.name); shift)
    groups, layout = _stationgroups(_withoutpriors(jones), t, layout)
    J = ObservedInstrumentModel(groups, refbasis, layout, t.pol, int.feedrotation_corrected)
    _checkgroups(J, x)
    return J
end

_checkstations(::AbstractJonesMatrix, names) = nothing
function _checkstations(J::PerStation, names)
    for s in keys(J.overrides)
        s in names || throw(ArgumentError("PerStation names the station $s, which is not in the table's stations"))
    end
    return nothing
end

# The priors with each parameter restricted to the stations whose term reads it.
# The stations, of `names`, whose term reads each parameter of `prior`.
_readers(::AbstractJonesMatrix, prior::NamedTuple, names) = map(_ -> names, prior)
function _readers(J::PerStation, prior::NamedTuple, names)
    return NamedTuple{keys(prior)}(map(p -> filter(s -> _reads(stationterm(J, s), p), names), keys(prior)))
end

_stationpriors(::AbstractJonesMatrix, prior, names) = prior
function _stationpriors(J::PerStation, prior::NamedTuple, names)
    restricted = map(keys(prior), values(prior), values(_readers(J, prior, names))) do p, d, readers
        isempty(readers) && throw(ArgumentError("no station's term reads the instrument parameter $p"))
        return length(readers) == length(names) ? d : _restrict(d, readers)
    end
    return NamedTuple{keys(prior)}(restricted)
end
function _restrict(d::ArrayPrior, stations)
    isnothing(d.stations) || return d
    return ArrayPrior(d.default_dist, d.override_dist, d.refant, d.centroid_station, d.gauge, Tuple(stations))
end
_restrict(d, stations) = d

# The terms of a `PerStation` model with the range of cells each one gives.
struct _StationGroups{T <: Tuple, C <: Tuple} <: AbstractJonesMatrix
    terms::T
    cells::C
end

# A `PerStation` model's terms, with the cells reordered so that each term's cells are contiguous:
# a traced loop over a contiguous range raises, and one over scattered cells does not.
_stationgroups(J::AbstractJonesMatrix, t, layout) = (J, layout)
function _stationgroups(J::PerStation, t::VisTable, layout::CellLayout)
    names = info(t).stations.name
    terms = (J.default, values(J.overrides)...)
    stationof = map(a -> names[a], layout.coords.ant)
    cells = (
        findall(s -> !haskey(J.overrides, s), stationof),
        map(s -> findall(==(s), stationof), keys(J.overrides))...,
    )
    used = map(!isempty, cells)
    any(used) || throw(ArgumentError("the table has no station cells"))
    T, C = terms[collect(used)], cells[collect(used)]
    ends = cumsum(map(length, C))
    ranges = map((e, c) -> (e - length(c) + 1):e, ends, C)
    return _StationGroups(Tuple(T), Tuple(ranges)), _permutecells(layout, reduce(vcat, C))
end

# `L` with its cells in the order `perm`.
function _permutecells(L::CellLayout, perm)
    newindex = invperm(perm)
    response = isnothing(L.response) ? nothing : _Columns{SMatrix{2, 2}}(map(col -> col[perm], L.response.cols))
    return CellLayout(L.coords[perm], newindex[L.cell1], newindex[L.cell2], map(e -> e[perm], L.entries), response)
end

# Every term of a `PerStation` model must give values of one type.
_checkgroups(J::ObservedInstrumentModel, x) = nothing
function _checkgroups(J::ObservedInstrumentModel{<:_StationGroups}, x)
    L = layout(J)
    params = _instrumentparams(J, (; instrument = x))
    types = map((f, c) -> typeof(_cellvalue(f, params, L, first(c))), _cellfunctions(J), instrument(J).cells)
    allequal(types) || throw(
        ArgumentError(
            "the terms of a PerStation model give values of different types, $(join(unique(types), ", ")); " *
                "write every term as a 2×2 matrix, or every term as a number"
        )
    )
    return nothing
end

"""
    ObservedIdealInstrument

An ideal instrument bound to a table: the sky's visibilities compared directly, converted to
the table's polarization in the global basis `CirBasis()`.
"""
struct ObservedIdealInstrument{P} <: AbstractInstrumentModel
    pol::P
end
EnzymeRules.inactive_type(::Type{<:ObservedIdealInstrument}) = true

set_array(::IdealInstrumentModel, t::VisTable) = (ObservedIdealInstrument(t.pol), NamedDist((;)))
set_array(::IdealInstrumentModel, tables::NamedTuple) =
    (FieldInstruments(map(t -> ObservedIdealInstrument(t.pol), tables)), NamedDist((;)))

# The instruments bound to the tables of several fields, keyed by field name, reading one set of
# parameters.
struct FieldInstruments{F <: NamedTuple} <: AbstractInstrumentModel
    fields::F
end

EnzymeRules.inactive_type(::Type{<:FieldInstruments}) = true

@inline apply_instrument(vis::NamedTuple, J::FieldInstruments, x) = map((v, j) -> apply_instrument(v, j, x), vis, J.fields)

@inline apply_instrument(vis, J::ObservedIdealInstrument, x) = _IdealValue{typeof(J.pol), CirBasis}(J.pol).(baseimage(vis))

# Instrument values are numbers or 2×2 matrices, held as one vector per component, so every
# read and write in a loop indexes a vector. Each loop is a plain loop on the host and, under
# `@trace`, a traced loop that Reactant raises to gathers, element-wise operations and stores; a
# loop raises only if it stores to each array once per iteration, hence a vector per component.
# `rgetindex` and `rsetindex!` index both.

# Values rebuilt as `T` from their component vectors `cols` (column-major for matrices).
struct _Columns{T, C <: Tuple}
    cols::C
end
_Columns{T}(cols) where {T} = _Columns{T, typeof(cols)}(cols)
Base.length(c::_Columns) = length(first(c.cols))

@inline _load(c::_Columns{<:Number}, i) = rgetindex(only(c.cols), i)
@inline _load(c::_Columns{<:StaticArraysCore.StaticMatrix}, i) = SMatrix{2, 2}(map(col -> rgetindex(col, i), c.cols))
@inline _load(c::_Columns{<:StokesParams}, i) = StokesParams(map(col -> rgetindex(col, i), c.cols)...)

@inline _store!(c::_Columns, v::Number, i) = (rsetindex!(only(c.cols), v, i); nothing)
@inline _store!(c::_Columns, v::StaticArraysCore.StaticMatrix, i) = (map((col, x) -> rsetindex!(col, x, i), c.cols, Tuple(v)); nothing)

# Storage for `n` values like `v`, in vectors like `like`.
_valuecolumns(::Number, like, n) = _Columns{Number}((similar(like, n),))
_valuecolumns(::StaticArraysCore.StaticMatrix, like, n) = _Columns{SMatrix{2, 2}}(ntuple(_ -> similar(like, n), Val(4)))
function _valuecolumns(v, like, n)
    throw(ArgumentError("an instrument term gives values of type $(typeof(v)); expected numbers or 2×2 matrices"))
end

# `n` values rebuilt as `T` from components stored back to back: component `k` of value `i` is
# `data[i + (k - 1) * n]`, the layout of a `FieldDimArray`'s storage. (A traced array cannot be
# indexed through a view at a traced index, so the components are not views.)
struct _Flat{T, D <: AbstractVector}
    data::D
    n::Int
end
_Flat{T}(data, n) where {T} = _Flat{T, typeof(data)}(data, n)
Base.length(c::_Flat) = c.n
@inline _load(c::_Flat{<:StokesParams}, i) = StokesParams(ntuple(k -> rgetindex(c.data, i + (k - 1) * c.n), Val(4))...)

# The sky's visibilities, and a vector to allocate like.
_viscolumns(v::AbstractArray{<:Number}) = _Columns{Number}((vec(v),))
_viscolumns(v::FieldDimArray{<:StokesParams}) = _Flat{StokesParams}(vec(parent(v)), length(v))
_viscolumns(v::AbstractArray{<:StokesParams}) = _Columns{StokesParams}(ntuple(k -> vec(getfield.(v, k)), Val(4)))
_template(c::_Columns) = first(c.cols)
_template(c::_Flat) = c.data

# The values over the table's `(record, channel)` shape: numbers, or a `FieldDimArray` of 2×2
# matrices.
_tablevalues(out::_Columns{<:Number}, sz) = reshape(only(out.cols), sz)
_tablevalues(out::_Columns{<:StaticArraysCore.StaticMatrix}, sz) = FieldDimArray{SMatrix{2, 2}}(reshape(vcat(out.cols...), sz..., 2, 2))

# The instrument parameters by name, as plain vectors.
@inline _instrumentparams(J::ObservedInstrumentModel, x) = map(parent ∘ siteparams, x.instrument)
@inline _instrumentparams(J::ObservedInstrumentModel{IdentityJones}, x) = (;)

# The Jones term of one cell from its parameters `xs` and coordinates `cs`, with the response
# `U` when the data need it (`R`): `J U`, or `Uᴴ J U` for data corrected for feed rotation (`C`).
struct _CellJones{I, R, C}
    instrument::I
end
function (f::_CellJones{I, R, C})(xs, cs, U) where {I, R, C}
    j = jonesmatrix(f.instrument, xs, cs)
    R || return j
    return C ? U' * j * U : j * U
end

_cellfunction(J::ObservedInstrumentModel, term) =
    _CellJones{typeof(term), !isnothing(layout(J).response), corrected(J)}(term)
_cellfunctions(J::ObservedInstrumentModel{<:_StationGroups}) = map(t -> _cellfunction(J, t), instrument(J).terms)

# The value of cell `k`: its parameters are entry `entries[p][k]` of each parameter `p`.
@inline function _cellvalue(f, params, L::CellLayout, k)
    xs = map((p, e) -> rgetindex(p, rgetindex(e, k)), params, L.entries)
    cs = map(col -> rgetindex(col, k), StructArrays.components(L.coords))
    return f(xs, cs, _cellresponse(L.response, k))
end
_cellresponse(::Nothing, k) = nothing
_cellresponse(U::_Columns, k) = _load(U, k)

# The Jones matrix `J U` of every cell, stored like `like`.
function _celljones(J::ObservedInstrumentModel, x, like)
    L = layout(J)
    f = _cellfunction(J, instrument(J))
    params = _instrumentparams(J, x)
    out = _valuecolumns(_cellvalue(f, params, L, 1), like, length(L.coords))
    _fillcells!(out, f, params, L, Base.OneTo(length(L.coords)))
    return out
end

# A `PerStation` model: each term fills its own cells.
function _celljones(J::ObservedInstrumentModel{<:_StationGroups}, x, like)
    L = layout(J)
    fs = _cellfunctions(J)
    params = _instrumentparams(J, x)
    out = _valuecolumns(_cellvalue(first(fs), params, L, 1), like, length(L.coords))
    foreach((f, cells) -> _fillcells!(out, f, params, L, cells), fs, instrument(J).cells)
    return out
end

function _fillcells!(out, f, params, L, cells)
    @trace track_numbers = false for m in 1:length(cells)
        k = rgetindex(cells, m)
        _store!(out, _cellvalue(f, params, L, k), k)
    end
    return nothing
end

@inline function apply_instrument(vis, J::ObservedInstrumentModel, x)
    L = layout(J)
    v = _viscolumns(baseimage(vis))
    like = _template(v)
    jc = _celljones(J, x, like)
    f = _CellValue{typeof(J.pol), typeof(refbasis(J))}(J.pol)
    cell1, cell2 = vec(L.cell1), vec(L.cell2)
    n = length(v)
    out = _valuecolumns(_value(f, jc, v, cell1, cell2, 1), like, n)
    @trace track_numbers = false for i in 1:n
        _store!(out, _value(f, jc, v, cell1, cell2, i), i)
    end
    return _tablevalues(out, size(baseimage(vis)))
end

@inline _value(f, jc, v, cell1, cell2, i) = f(_load(jc, rgetindex(cell1, i)), _load(v, i), _load(jc, rgetindex(cell2, i)))

# The value of the table's polarization `pol` from the Jones matrices `j1`, `j2` of a value's two
# stations and the sky's visibility `v`, in the global basis `B`.
struct _CellValue{P, B}
    pol::P
end
@inline (f::_CellValue{P, B})(j1, v, j2) where {P, B} = _projectvalue(f.pol, _sandwich(j1, _idealvalue(v, B), j2))

# The table's polarization from the sky's visibility `v` alone, in the global basis `B`.
struct _IdealValue{P, B}
    pol::P
end
@inline (f::_IdealValue{P, B})(v) where {P, B} = _projectvalue(f.pol, _idealvalue(v, B))

_idealvalue(v::Number, ::Type) = v
_idealvalue(v::StokesParams, ::Type{B}) where {B} = SMatrix{2, 2}(CoherencyMatrix{B, B}(v))

_projectvalue(::IPol, v::Number) = v
_projectvalue(::IPol, m::StaticArraysCore.StaticMatrix) = (m[1, 1] + m[2, 2]) / 2
_projectvalue(::Correlation{I, J}, m::StaticArraysCore.StaticMatrix) where {I, J} = m[I, J]
_projectvalue(::Correlations, m::StaticArraysCore.StaticMatrix) = m
function _projectvalue(pol, v)
    throw(ArgumentError("model visibilities of type $(typeof(v)) cannot give the table's polarization $(_polname(pol))"))
end

_sandwich(a, x, b) = a * x * _adjoint(b)
_adjoint(b::Number) = conj(b)
_adjoint(b::StaticArraysCore.StaticMatrix) = b'
# A number between matrix terms is that number times the identity.
_sandwich(a::StaticArraysCore.StaticMatrix, x::Number, b::StaticArraysCore.StaticMatrix) = x * (a * b')

"""
    forward_jones(J::ObservedInstrumentModel, x, t::VisTable)

The Jones matrix `J_a U_a` of every station, time stamp and channel at the parameters `x`, as an
array over `(Ant, Ts, Fr)`; stations without a record at a time hold `NaN`.
"""
function forward_jones(J::ObservedInstrumentModel, x, t::VisTable)
    L = layout(J)
    cells = _celljones(J, x, ComplexF64[])
    jc = [_load(cells, k) for k in 1:length(cells)]
    T = eltype(jc)
    out = fill(_nanof(T), length(info(t).stations), length(times(t)), size(t, 2))
    d = domain(t)
    for I in CartesianIndices(L.cell1), (a, c) in ((d.ant1[I[1]], L.cell1[I]), (d.ant2[I[1]], L.cell2[I]))
        out[a, d.ts[I[1]], I[2]] = jc[c]
    end
    return DD.DimArray(out, (Ant(info(t).stations.name), Ts(times(t)), Fr(info(t).channels.frequency)))
end

_nanof(::Type{T}) where {T <: Number} = T(NaN)
_nanof(::Type{T}) where {T <: StaticArraysCore.StaticMatrix} = T(ntuple(_ -> NaN, length(T)))
