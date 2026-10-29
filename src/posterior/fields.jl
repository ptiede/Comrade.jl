export Field

"""
    Field(name, sky::AbstractSkyModel, data...)

One source of a multi-field [`VLBIPosterior`](@ref): its `name`, its sky model `sky`, and its
data products, all formed from one visibility table of that source.
"""
struct Field{S <: AbstractSkyModel, D <: Tuple}
    name::Symbol
    sky::S
    data::D
    function Field(name, sky::AbstractSkyModel, data::Union{VisTable, ClosureTable}...)
        isempty(data) && throw(ArgumentError("the field $name has no data products"))
        return new{typeof(sky), typeof(data)}(Symbol(name), sky, data)
    end
end

"""
    VLBIPosterior(instrumentmodel::AbstractInstrumentModel, fields::Field...; admode)

The posterior of several sources, each a [`Field`](@ref) with its own sky model and data,
observed through one instrument. The parameters are
`(; sky = (; name₁ = …, name₂ = …), instrument = …)`, the sky parameters keyed by field name.

The instrument parameters are built on the measured values of all fields, on the time axis of
the earliest of their epochs, which the cell time `c.Ti` of every field reads. A parameter is
shared by the fields wherever its time segments cover the data of several: a `TrackSeg()`
leakage is shared by every field, while a `ScanSeg()` gain is shared only where the fields'
scans overlap. A term that depends on direction reads it from the cell: `c.ra` and `c.dec`
are the direction of the field's source, and `c.elevation`, `c.azimuth` and `c.parallactic`
its geometry at each station. Data products of one field must come from one table.

```julia
post = VLBIPosterior(InstrumentModel(G * D),
    Field(:SgrA, sky_sgra, dvis_sgra),
    Field(:J1924, sky_cal, dvis_cal),
)
```
"""
@noinline function VLBIPosterior(
        instrumentmodel::AbstractInstrumentModel, field::Field, fields::Field...;
        admode = EnzymeCore.set_runtime_activity(EnzymeCore.set_strong_zero(EnzymeCore.Reverse))
    )
    fs = (field, fields...)
    names = map(f -> f.name, fs)
    allunique(names) || throw(ArgumentError("the fields need distinct names, got $names"))
    tables = NamedTuple{names}(map(f -> _commonvis(f.data), fs))
    data = NamedTuple{names}(map(f -> f.data, fs))

    int, intprior = set_array(instrumentmodel, tables)
    observed = map((f, t) -> set_array(f.sky, t), fs, values(tables))
    sky = FieldSkies(NamedTuple{names}(map(first, observed)))
    skyprior = NamedDist(NamedTuple{names}(map(last, observed)))
    total_prior = combine_prior(skyprior, intprior)
    ls = map(d -> Tuple(map(makelikelihood, d)), data)

    return VLBIPosterior{
        typeof(data), Nothing, typeof(ls), Tuple{}, typeof(total_prior),
        typeof(sky), typeof(int), typeof(admode),
    }(data, nothing, ls, (), total_prior, sky, int, admode)
end

# The observed sky models of several fields, keyed by field name.
struct FieldSkies{F <: NamedTuple} <: AbstractSkyModel
    fields::F
end

idealmaps(D::VisData, m::FieldSkies, x) = _fieldmaps(D, m, x)
idealmaps(D::DualData, m::FieldSkies, x) = _fieldmaps(D, m, x)

function _fieldmaps(D, m::FieldSkies, x)
    maps = map((sm, xs) -> idealmaps(D, sm, (; sky = xs)), m.fields, x.sky)
    return map(first, maps), map(last, maps)
end

skymodel(m::FieldSkies, x) = map(skymodel, m.fields, x)
skymodel(m::FieldSkies) = map(skymodel, m.fields)
domain(m::FieldSkies) = map(domain, m.fields)

@inline logdensityofvis(lklhds::NamedTuple, vis::NamedTuple) = sum(values(map(logdensityofvis, lklhds, vis)))
@inline logdensityofimg(::Tuple{}, img::NamedTuple) = sum(values(map(i -> logdensityofimg((), i), img)))
