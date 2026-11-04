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
@inline logdensityofimg(lklhds::NamedTuple, img::NamedTuple) = sum(values(map(logdensityofimg, lklhds, img)))
