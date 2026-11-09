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

# The observed sky model of one source on several tables, keyed by table name; every table
# reads the same sky parameters.
struct TableSkies{F <: NamedTuple} <: AbstractSkyModel
    tables::F
end

idealmaps(D::VisData, m::TableSkies, x) = _tablemaps(D, m, x)
idealmaps(D::DualData, m::TableSkies, x) = _tablemaps(D, m, x)

function _tablemaps(D, m::TableSkies, x)
    maps = map(sm -> idealmaps(D, sm, (; sky = x.sky)), m.tables)
    return map(first, maps), map(last, maps)
end

skymodel(m::TableSkies, x) = skymodel(first(m.tables), x)
skymodel(m::TableSkies) = skymodel(first(m.tables))
domain(m::TableSkies) = map(domain, m.tables)

@inline logdensityofvis(lklhds::NamedTuple, vis::NamedTuple) = sum(values(map(logdensityofvis, lklhds, vis)))
@inline logdensityofimg(::Tuple{}, img::NamedTuple) = sum(values(map(i -> logdensityofimg((), i), img)))
@inline logdensityofimg(lklhds::NamedTuple, img::NamedTuple) = sum(values(map(logdensityofimg, lklhds, img)))
