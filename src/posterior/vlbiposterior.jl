struct VLBIPosterior{DV, DI, LV, LI, P, MS <: AbstractSkyModel, MI <: AbstractInstrumentModel, ADMode <: Union{Nothing, EnzymeCore.Mode}, D} <: AbstractVLBIPosterior
    data::DV
    dataimg::DI
    lklhds::LV
    lklhdsimg::LI
    prior::P
    skymodel::MS
    instrumentmodel::MI
    admode::ADMode
    # The layout or sharding `todevice` placed the posterior with; `nothing` on the host.
    placement::D
end

(post::VLBIPosterior)(θ) = logdensityof(post, θ)
admode(post::VLBIPosterior) = post.admode

@inline function datatype(::Type{<:VLBIPosterior{DV, DI}}) where {DV, DI}
    return datatype(DV, DI)
end


"""
    VLBIPosterior(dataproducts...; sky, instrument = IdealInstrumentModel(),
                  admode = set_runtime_activity(set_strong_zero(Reverse)))

The posterior of the sky model `sky` and the instrument model `instrument` given the data
products. The parameters are `(; sky = …, instrument = …)`.

The data products are visibility-domain products, [`VisTable`](@ref)s and
[`ClosureTable`](@ref)s formed from one visibility table, on whose records the models are
evaluated, and image-domain products such as [`ImgNormalData`](@ref), compared with the image of
`sky`. [`Comrade.datadomain`](@ref) sorts them.

To enable automatic differentiation, the `admode` keyword argument can be set to any `EnzymeCore.Mode` type 
of if no AD is desired then `nothing`. The default is `Enzyme.set_runtime_activity(Enzyme.Reverse)` 
which should work for essentially every problem. Note that runtime activity does have a perfomance cost, 
and as Enzyme and Comrade matures we expect this to not need runtime activity.

# Example
```julia
t = average(VisTable(uvd, IPol()), ByScan())
dlcamp, dcphase = logclosure_amplitudes(t), closure_phases(t)

function sky(θ, grid, metadata)
    (; r, a) = θ
    m = stretched(ExtendedRing(a), r, r)
    return m
end

skyprior = (r = VLBIUniform(μas2rad(10.0), μas2rad(30.0)), a = VLBIUniform(1.0, 10.0))
g  = spatialgrid(μas2rad(100.0), μas2rad(100.0), 256, 256)
skym = SkyModel(sky, skyprior, g)

G = Jones((x, c) -> exp(complex(x.lg, x.pg));
          lg = ArrayPrior(IIDStationPrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
          pg = ArrayPrior(IIDStationPrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))))
          )
intmodel = InstrumentModel(G)

post = VLBIPosterior(t, ImgNormalData(flux, 0.6, 0.1); sky = skym, instrument = intmodel)
```

# Several fields

    VLBIPosterior(data::NamedTuple; sky::NamedTuple, instrument, admode)

The posterior of several sources observed through one instrument. `sky` holds the sky model of
each field, `data` the data products of each field, with the same keys: a tuple of products, or
one product. The data products of one field must come from one table, and every field needs a
visibility-domain product. The parameters are `(; sky = (; name₁ = …, name₂ = …), instrument = …)`,
the sky parameters keyed by field name.

The instrument parameters are built on the measured values of all fields, on the time axis of
the earliest of their epochs, which the cell time `c.Ti` of every field reads. A parameter is
shared by the fields wherever its time segments cover the data of several: a `TrackSeg()`
leakage is shared by every field, while a `ScanSeg()` gain is shared only where the fields'
scans overlap. A term that depends on direction reads it from the cell: `c.ra` and `c.dec`
are the direction of the field's source, and `c.elevation`, `c.azimuth` and `c.parallactic`
its geometry at each station.

A self-calibrated instrument (`selfcal = true` in [`InstrumentModel`](@ref)) needs one table of
visibilities or coherencies per field, and fields that share no baseline at one time; the
station phase is solved against the data of all fields.

```julia
post = VLBIPosterior((; SgrA = (dvis_sgra, ImgNormalData(flux, 2.4, 0.1)), J1924 = dvis_cal);
                     sky = (; SgrA = sky_sgra, J1924 = sky_cal), instrument = InstrumentModel(G * D))
```

# Several tables of one source

    VLBIPosterior(data::NamedTuple; sky::AbstractSkyModel, instrument, admode)

The posterior of one source observed in several tables, such as different bands or epochs:
`data` holds the data products of each table, a tuple of products or one product, and the one
sky model is evaluated on the records and channels of each table, with the same parameters
`(; sky = …, instrument = …)`. A table is not padded to the channels of the others, unlike a
table joined with `vcat`. The instrument is bound to the tables as to several fields, and a
frequency-dependent sky reads each table's frequencies from its points.

```julia
post = VLBIPosterior((; Q = (dvis86, dcp86), A = dvis230); sky = skym, instrument = InstrumentModel(G))
```
"""
@noinline function VLBIPosterior(
        data...; sky, instrument = IdealInstrumentModel(),
        admode = EnzymeCore.set_runtime_activity(EnzymeCore.set_strong_zero(EnzymeCore.Reverse))
    )
    instrument isa AbstractInstrumentModel ||
        throw(ArgumentError("the instrument must be an AbstractInstrumentModel, got a $(typeof(instrument))"))
    return _posterior(sky, instrument, data, admode)
end

function _posterior(skymodel::AbstractSkyModel, instrumentmodel, data::Tuple, admode)
    vis, img = _splitdomains(data)
    t = _commonvis(instrumentmodel, vis, data)
    int, intprior = observe(instrumentmodel, t)
    sky, skyprior = observe(skymodel, t)
    total_prior = combine_prior(skyprior, intprior)

    ls = map(makelikelihood, vis)
    ils = map(makelikelihood, img)
    dataimg = _imagedata(img)
    return VLBIPosterior{
        typeof(vis), typeof(dataimg), typeof(ls), typeof(ils), typeof(total_prior),
        typeof(sky), typeof(int), typeof(admode), Nothing,
    }(vis, dataimg, ls, ils, total_prior, sky, int, admode, nothing)
end

function _posterior(skies::NamedTuple, instrumentmodel, (data,)::Tuple{NamedTuple}, admode)
    names = keys(skies)
    isempty(names) && throw(ArgumentError("a multi-field posterior needs at least one field"))
    issetequal(names, keys(data)) || throw(
        ArgumentError("the sky models and the data need the same field names, got sky models for $names and data for $(keys(data))")
    )
    for (n, s) in pairs(skies)
        s isa AbstractSkyModel || throw(ArgumentError("the sky model of the field $n must be an AbstractSkyModel, got a $(typeof(s))"))
    end
    return _severaltables(instrumentmodel, NamedTuple{names}(data), admode, "field") do tables
        observed = map(observe, skies, tables)
        return FieldSkies(map(first, observed)), NamedDist(map(last, observed))
    end
end

function _posterior(skymodel::AbstractSkyModel, instrumentmodel, (data,)::Tuple{NamedTuple}, admode)
    isempty(data) && throw(ArgumentError("a NamedTuple of data products needs at least one table"))
    return _severaltables(instrumentmodel, data, admode, "table") do tables
        observed = map(t -> observe(skymodel, t), tables)
        return TableSkies(map(first, observed)), last(first(observed))
    end
end

# The posterior of the data products of several tables, keyed like `data`; `observesky(tables)`
# returns the observed sky and its prior. `kind` names a key in error messages.
function _severaltables(observesky, instrumentmodel, data::NamedTuple, admode, kind)
    split = map(_splitdomains ∘ _products, data)
    vis, img = map(first, split), map(last, split)
    names = keys(vis)
    tables = NamedTuple{names}(map(n -> _commonvis(instrumentmodel, vis[n], _products(data[n]), "$kind $n"), names))

    int, intprior = observe(instrumentmodel, tables)
    sky, skyprior = observesky(tables)
    total_prior = combine_prior(skyprior, intprior)
    ls = map(Base.Fix1(map, makelikelihood), vis)
    dataimg = _imagedata(img)
    ils = _imagelikelihoods(dataimg)

    return VLBIPosterior{
        typeof(vis), typeof(dataimg), typeof(ls), typeof(ils), typeof(total_prior),
        typeof(sky), typeof(int), typeof(admode), Nothing,
    }(vis, dataimg, ls, ils, total_prior, sky, int, admode, nothing)
end

_posterior(::NamedTuple, _, data::Tuple, _) = throw(
    ArgumentError(
        "a NamedTuple of sky models takes one NamedTuple of data products with the same keys, " *
            "got $(length(data)) data argument(s): $(_typenames(data))"
    )
)
_posterior(sky, _, _, _) =
    throw(ArgumentError("the sky must be an AbstractSkyModel or a NamedTuple of them, got a $(typeof(sky))"))

_products(d::Tuple) = d
_products(d) = (d,)

# The visibility-domain and the image-domain products of `data`, each in the order given.
_splitdomains(data::Tuple) = (_indomain(VisDomain(), data), _indomain(ImageDomain(), data))
_indomain(_, ::Tuple{}) = ()
function _indomain(D, data::Tuple)
    d = first(data)
    rest = _indomain(D, Base.tail(data))
    return datadomain(d) === D ? (d, rest...) : rest
end

_imagedata(::Tuple{}) = nothing
_imagedata(img::Tuple) = img
_imagedata(::NamedTuple{K, <:Tuple{Vararg{Tuple{}}}}) where {K} = nothing
_imagedata(img::NamedTuple) = img
_imagelikelihoods(::Nothing) = ()
_imagelikelihoods(img::NamedTuple) = map(Base.Fix1(map, makelikelihood), img)

_typenames(::Tuple{}) = "no data products"
_typenames(data) = join(map(d -> string(nameof(typeof(d))), data), ", ")

_vis(t::VisTable) = t
_vis(ct::ClosureTable) = getfield(ct, :vis)

# The visibility table every product in `vis` is formed from; the model is evaluated on its records.
function _commonvis(vis)
    t = _vis(first(vis))
    all(d -> _samerecords(_vis(d), t), vis) || throw(
        ArgumentError("every data product of a posterior must come from one table of visibilities: same records and channels")
    )
    return t
end

# `_commonvis(vis)` after checking the products against `instrumentmodel`; `data` holds all
# products of the table, `key` names it in error messages.
function _commonvis(instrumentmodel, vis, data, key = nothing)
    isempty(vis) && throw(
        ArgumentError("$(isnothing(key) ? "a posterior" : "the $key") needs a visibility-domain data product (a VisTable or ClosureTable), got $(_typenames(data))")
    )
    _checkdata(instrumentmodel, vis, key)
    return _commonvis(vis)
end

function _samerecords(a::VisTable, b::VisTable)
    table(t) = (ComradeBase.coords(domain(t)), info(t).channels.frequency, info(t).epoch, info(t).source)
    return a === b || table(a) == table(b)
end

function combine_prior(skyprior, instrumentmodelprior)
    return NamedDist((sky = skyprior, instrument = instrumentmodelprior))
end

function combine_prior(skymodel, ::NamedDist{()})
    return NamedDist((; sky = skymodel))
end

function combine_prior(::NamedDist{()}, instrumentmodel)
    return NamedDist((; instrument = instrumentmodel))
end

function combine_prior(::NamedDist{()}, ::NamedDist{()})
    return NamedDist()
end


function Base.show(io::IO, mime::MIME"text/plain", post::VLBIPosterior)
    printstyled(io, "VLBIPosterior"; bold = true, color = :light_magenta)
    println(io)
    show(io, mime, post.skymodel)
    println()
    show(io, mime, post.instrumentmodel)
    println()
    printstyled(io, "Data Products: ", color = :light_green)
    println(io, _datadescription(post.data))
    isnothing(post.dataimg) && return nothing
    printstyled(io, "Image Data Products: ", color = :light_green)
    return println(io, _datadescription(post.dataimg))
end


# A self-calibrated model is solved against the one table of each field it is fit to.
_checkdata(::AbstractInstrumentModel, data, key = nothing) = nothing
function _checkdata(m::InstrumentModel, data, key = nothing)
    m.selfcal || return nothing
    vis = filter(d -> d isa VisTable && !(eltype(measurement(d)) <: Real), collect(data))
    of = isnothing(key) ? "" : " for the $key"
    length(data) == 1 && length(vis) == 1 || throw(
        ArgumentError(
            "a self-calibrated instrument model needs exactly one table of visibilities or coherencies$of, " *
                "got $(length(data)) table(s), $(length(vis)) of them complex"
        )
    )
    return nothing
end

"""
    simulate_observation([rng::Random.AbstractRNG], post::VLBIPosterior, θ; add_thermal_noise=true)

Create a simulated observation using the posterior and its data `post` using the parameter
values `θ`. In Bayesian terminology this is a draw from the posterior predictive distribution.

If `add_thermal_noise` is true then baseline based thermal noise is added. Otherwise, we just
return the model visibilities. A multi-field posterior gives the data products of each field,
keyed by field name.
"""
function simulate_observation(rng::Random.AbstractRNG, post::VLBIPosterior, θ; add_thermal_noise = true)
    return _simulated(rng, post.data, last(forward_model(post, θ)), add_thermal_noise)
end

# One tuple of simulated data products per field.
_simulated(rng, data::NamedTuple, vis::NamedTuple, noise) = map((d, v) -> _simulated(rng, d, v, noise), data, vis)
function _simulated(rng, data::Tuple, vis, noise)
    v = _observedvalues(rng, vis, _commonvis(data), noise)
    return map(d -> _withmodel(d, v), data)
end

simulate_observation(post::VLBIPosterior, θ; add_thermal_noise = true) = simulate_observation(Random.default_rng(), post, θ; add_thermal_noise)

"""
    selfcal_phase(post::VLBIPosterior, θ)

The station phase `φ = φ̂ + residual` at `θ` of a posterior whose instrument model
self-calibrates the station phase (see [`InstrumentModel`](@ref)), as a `StationArray`.
"""
function selfcal_phase(post::VLBIPosterior, θ)
    _, vis = idealmaps(datatype(typeof(post)), skymodel(post), θ)
    return selfcal_phase(instrumentmodel(post), vis, _instrumentvalues(θ))
end

_datadescription(data::Tuple) = join(map(_datadescription, data), ", ")
_datadescription(data::NamedTuple) = join(map((n, d) -> string(n, ": ", _datadescription(d)), keys(data), values(data)), "; ")
_datadescription(t::VisTable) = string(eltype(measurement(t)) <: Real ? "amplitudes" : "visibilities", " (", _polname(t.pol), ")")
_datadescription(ct::ClosureTable) = ct.product isa ClosurePhases ? "closure phases" : "log closure amplitudes"
_datadescription(d) = string(nameof(typeof(d)))

"""
    residuals(post::VLBIPosterior, θ)

The residuals of each data product in `post` at the parameter values `θ`: tables of the same
kind (a [`VisTable`](@ref) or [`ClosureTable`](@ref)) whose measurements are the data minus the
model. Closure phase residuals are wrapped to `(-π, π]`. A multi-field posterior gives the
residuals of each field, keyed by field name.
"""
function residuals(post::VLBIPosterior, p)
    return _residuals(post.data, last(forward_model(post, p)))
end

# One tuple of residual tables per field.
_residuals(data::NamedTuple, vis::NamedTuple) = map(_residuals, data, vis)
_residuals(data::Tuple, vis) = map(d -> _residual(d, baseimage(vis)), data)

function _residual(t::VisTable, vis)
    return VisTable(t.pol, measurement(t) .- measurement(_withmodel(t, vis)), weight(t))
end
function _residual(ct::ClosureTable, vis)
    r = measurement(ct) .- measurement(_withmodel(ct, vis))
    r = ct.product isa ClosurePhases ? atan.(sin.(r), cos.(r)) : r
    return ClosureTable(ct.product, r, getfield(ct, :vis), getfield(ct, :design), getfield(ct, :covariance))
end

"""
    chi2(post::AbstractVLBIPosterior, p; reduce=false)

Returns a tuple of the chi-squared values for each data product in the posterior `post` given the parameters `p`.
Note that the chi-square is not reduced. If you want to reduce it by dividing by the number of data
points then set `reduce=true`. 

Note that reduce doesn't take into account the number of model parameters
which for non-linear models is difficult to define globally.
"""
function chi2(post::AbstractVLBIPosterior, p; reduce = false)
    return _chi2s(residuals(post, p), reduce)
end

_chi2s(res::NamedTuple, reduce) = map(r -> _chi2s(r, reduce), res)
_chi2s(res::Tuple, reduce) = map(r -> reduce ? _chi2(r) / ndata(r) : _chi2(r), res)

# Σ |residual|² / σ² over the measured values.
_chi2(res::VisTable) = sum(_dense(weight(res)) .* abs2.(_dense(measurement(res))))
_chi2(res::ClosureTable) = sum(abs2.(baseimage(measurement(res)) ./ baseimage(noise(res))))

"""
    ndata(t)

The number of real data values in the table `t`: two per measured complex value, one per
amplitude or closure.
"""
ndata(t::VisTable) = (eltype(measurement(t)) <: Real ? 1 : 2) * count(>(0), _dense(weight(t)))
ndata(ct::ClosureTable) = length(ct)
