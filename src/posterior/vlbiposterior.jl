struct VLBIPosterior{DV, DI, LV, LI, P, MS <: AbstractSkyModel, MI <: AbstractInstrumentModel, ADMode <: Union{Nothing, EnzymeCore.Mode}} <: AbstractVLBIPosterior
    data::DV
    dataimg::DI
    lklhds::LV
    lklhdsimg::LI
    prior::P
    skymodel::MS
    instrumentmodel::MI
    admode::ADMode
end

(post::VLBIPosterior)(θ) = logdensityof(post, θ)
admode(post::VLBIPosterior) = post.admode

@inline function datatype(::Type{<:VLBIPosterior{DV, DI}}) where {DV, DI}
    return datatype(DV, DI)
end


"""
    VLBIPosterior(skymodel::SkyModel, instumentmodel::InstrumentModel, 
                  dataproducts...; 
                  admode=set_runtime_activity(set_strong_zero(Reverse)))

The posterior of the sky model `skymodel` and the instrument model `instrumentmodel` (default
[`IdealInstrumentModel`](@ref)) given the data products: [`VisTable`](@ref)s and
[`ClosureTable`](@ref)s formed from one visibility table. The models are evaluated on that
table's records.

To enable automatic differentiation, the `admode` keyword argument can be set to any `EnzymeCore.Mode` type 
of if no AD is desired then `nothing`. The default is `Enzyme.set_runtime_activity(Enzyme.Reverse)` 
which should work for essentially every problem. Note that runtime activity does have a perfomance cost, 
and as Enzyme and Comrade matures we expect this to not need runtime activity.

# Example
```julia
t = average(VisTable(uvd, IPol()), ByScan())
dlcamp, dcphase = log_closure_amplitudes(t), closure_phases(t)

function sky(θ, metadata)
    (; r, a) = θ
    m = stretched(ExtendedRing(a), r, r)
    return m
end

skyprior = (r = VLBIUniform(μas2rad(10.0), μas2rad(30.0)), a = VLBIUniform(1.0, 10.0))
g  = spatialgrid(μas2rad(100.0), μas2rad(100.0), 256, 256)
skym = SkyModel(sky, skyprior, g)

G = Jones((x, c) -> exp(complex(x.lg, x.pg));
          lg = ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1))),
          pg = ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))))
          )
intmodel = InstrumentModel(G)

post = VLBIPosterior(skym, intmodel, t)
```
"""
@noinline function VLBIPosterior(
        skymodel::AbstractSkyModel,
        instrumentmodel::AbstractInstrumentModel,
        data::Union{VisTable, ClosureTable}...;
        imgdata = nothing,
        admode = EnzymeCore.set_runtime_activity(EnzymeCore.set_strong_zero(EnzymeCore.Reverse))
    )

    t = _commonvis(data)
    _checkdata(instrumentmodel, data)
    int, intprior = set_array(instrumentmodel, t)
    sky, skyprior = set_array(skymodel, t)
    total_prior = combine_prior(skyprior, intprior)

    ls = Tuple(map(makelikelihood, data))
    if !isnothing(imgdata)
        ils = Tuple(map(makelikelihood, imgdata))
    else
        ils = ()
    end

    return VLBIPosterior{
        typeof(data), typeof(imgdata), typeof(ls), typeof(ils), typeof(total_prior),
        typeof(sky), typeof(int), typeof(admode),
    }(data, imgdata, ls, ils, total_prior, sky, int, admode)
end

VLBIPosterior(
    skymodel::AbstractSkyModel, dataproducts::Union{VisTable, ClosureTable}...;
    admode = EnzymeCore.set_runtime_activity(EnzymeCore.set_strong_zero(EnzymeCore.Reverse)), kwargs...
) =
    VLBIPosterior(skymodel, IdealInstrumentModel(), dataproducts...; admode, kwargs...)

_vis(t::VisTable) = t
_vis(ct::ClosureTable) = getfield(ct, :vis)

# The visibility table every data product is formed from; the model is evaluated on its records.
function _commonvis(data)
    t = _vis(first(data))
    for d in data
        _samerecords(_vis(d), t) || throw(
            ArgumentError("every data product of a posterior must come from one table of visibilities: same records and channels")
        )
    end
    return t
end

function _samerecords(a::VisTable, b::VisTable)
    a === b && return true
    ca, cb = ComradeBase.coords(domain(a)), ComradeBase.coords(domain(b))
    return size(a) == size(b) && ca.u == cb.u && ca.v == cb.v && ca.Ti == cb.Ti && ca.ant1 == cb.ant1 &&
        ca.ant2 == cb.ant2 && info(a).channels.frequency == info(b).channels.frequency &&
        info(a).epoch == info(b).epoch && info(a).source == info(b).source
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
    return println(io, _datadescription(post.data))
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
    v = _storage(baseimage(vis))
    t = _commonvis(data)
    if noise
        w = _dense(weight(t))
        v = rand(rng, _noiselikelihood(eltype(measurement(t)), v, w))
    end
    return map(d -> _withmodel(d, v), data)
end

_noiselikelihood(::Type{<:StaticArraysCore.StaticMatrix}, vis, w) = CoherencyLikelihood(vis, w)
_noiselikelihood(::Type, vis, w) = ComplexVisLikelihood(vis, w)
simulate_observation(post::VLBIPosterior, θ; add_thermal_noise = true) = simulate_observation(Random.default_rng(), post, θ; add_thermal_noise)

_datadescription(data::Tuple) = join(map(_datadescription, data), ", ")
_datadescription(data::NamedTuple) = join(map((n, d) -> string(n, ": ", _datadescription(d)), keys(data), values(data)), "; ")
_datadescription(t::VisTable) = string(eltype(measurement(t)) <: Real ? "amplitudes" : "visibilities", " (", _polname(t.pol), ")")
_datadescription(ct::ClosureTable) = ct.product isa ClosurePhases ? "closure phases" : "log closure amplitudes"

# `d` with its measurements replaced by those formed from the model visibilities `vis`.
function _withmodel(t::VisTable, vis)
    m = _tablevalues(eltype(measurement(t)), vis)
    w = _dense(weight(t))
    return VisTable(t.pol, IntensityMap(_elements(t.pol, _zeroed(m, w)), domain(t)), weight(t), stations(t))
end
function _withmodel(ct::ClosureTable, vis)
    x = _closurevalue.(Ref(ct.product), vis[_measuredvalues(getfield(ct, :vis))])
    m = IntensityMap(getfield(ct, :design) * x, domain(ct))
    return ClosureTable(ct.product, m, getfield(ct, :vis), getfield(ct, :design), getfield(ct, :covariance))
end

_tablevalues(::Type{<:Real}, vis) = abs.(vis)
_tablevalues(::Type, vis) = vis

# Values without a measurement are stored as 0.
_zeroed(m::AbstractArray{<:Number}, w) = ifelse.(w .> 0, m, zero(eltype(m)))
function _zeroed(m::AbstractArray{<:StaticArraysCore.StaticMatrix}, w)
    d = reshape(reinterpret(eltype(eltype(m)), collect(m)), 2, 2, size(m)...)
    return ifelse.(w .> 0, permutedims(d, (3, 4, 1, 2)), zero(eltype(d)))
end

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
    model = _withmodel(t, vis)
    r = _dense(measurement(t)) .- _dense(measurement(model))
    return VisTable(t.pol, IntensityMap(_elements(t.pol, r), domain(t)), weight(t), stations(t))
end
function _residual(ct::ClosureTable, vis)
    model = _withmodel(ct, vis)
    r = baseimage(measurement(ct)) .- baseimage(measurement(model))
    r = ct.product isa ClosurePhases ? atan.(sin.(r), cos.(r)) : r
    return ClosureTable(ct.product, IntensityMap(r, domain(ct)), getfield(ct, :vis), getfield(ct, :design), getfield(ct, :covariance))
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
