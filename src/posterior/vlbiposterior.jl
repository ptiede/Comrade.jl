struct VLBIPosterior{D, T, P, MS <: AbstractSkyModel, MI <: AbstractInstrumentModel, ADMode <: Union{Nothing, EnzymeCore.Mode}} <: AbstractVLBIPosterior
    data::D
    lklhds::T
    prior::P
    skymodel::MS
    instrumentmodel::MI
    admode::ADMode
end

(post::VLBIPosterior)(θ) = logdensityof(post, θ)
admode(post::VLBIPosterior) = post.admode


"""
    VLBIPosterior(skymodel::SkyModel, instumentmodel::InstrumentModel, 
                  dataproducts::EHTObservationTable...; 
                  admode=set_runtime_activity(Reverse))

Creates a VLBILikelihood using the `skymodel` its related metadata `skymeta`
and the `instrumentmodel` and its metadata `instumentmeta`. The `model` is a 
function that converts from parameters `θ` to a Comrade
AbstractModel which can be used to compute [`visibilitymap`](@ref) and a set of
`metadata` that is used by `model` to compute the model.

To enable automatic differentiation, the `admode` keyword argument can be set to any `EnzymeCore.Mode` type 
of if no AD is desired then `nothing`. The default is `Enzyme.set_runtime_activity(Enzyme.Reverse)` 
which should work for essentially every problem. Note that runtime activity does have a perfomance cost, 
and as Enzyme and Comrade matures we expect this to not need runtime activity.

# Warning

The `model` itself must be a two argument function where the first argument is the set
of model parameters and the second is a container that holds all the additional
information needed to construct the model. An example of this is when the model
needs some precomputed cache to define the model.

# Example
```julia
dlcamp, dcphase = extract_table(obs, LogClosureAmplitude(), ClosurePhases())
array = arrayconfiguration(dlcamp)

function sky(θ, metadata)
    (; r, a) = θ
    m = stretched(ExtendedRing(a), r, r)
    return m
end

skyprior = (r = Uniform(μas2rad(10.0), μas2rad(30.0)), a = Uniform(1.0, 10.0))
g  = imagepixels(μas2rad(100.0), μas2rad(100.0), 256, 256)
skym = SkyModel(sky, skyprior, g)

G = SingleStokesGain(x->exp(x.lg + 1im*x.pg))
intprior = (lg = ArrayPrior(IIDSitePrior(ScanSeg(), Normal(0.0, 0.1))),
            pg = ArrayPrior(IIDSitePrior(ScanSeg(), DiagVonMises(0.0, inv(π^2))))
            )
intmodel = InstrumentModel(G, intprior, array)

post = VLBIPosterior(skym, intmodel, dlcamp, dcphase)
```
"""
function VLBIPosterior(
        skymodel::AbstractSkyModel,
        instrumentmodel::AbstractInstrumentModel,
        dataproducts::EHTObservationTable...;
        admode = EnzymeCore.set_runtime_activity(EnzymeCore.Reverse)
    )


    array = arrayconfig(dataproducts[begin])
    int, intprior = set_array(instrumentmodel, array)
    sky, skyprior = set_array(skymodel, array)
    total_prior = combine_prior(skyprior, intprior)

    ls = Tuple(map(makelikelihood, dataproducts))

    return VLBIPosterior{
        typeof(dataproducts), typeof(ls), typeof(total_prior),
        typeof(sky), typeof(int), typeof(admode),
    }(dataproducts, ls, total_prior, sky, int, admode)
end

VLBIPosterior(
    skymodel::AbstractSkyModel, dataproducts::EHTObservationTable...;
    admode = EnzymeCore.set_runtime_activity(EnzymeCore.Reverse)
) =
    VLBIPosterior(skymodel, IdealInstrumentModel(), dataproducts...; admode)

function combine_prior(skyprior, instrumentmodelprior)
    return NamedDist((sky = skyprior, instrument = instrumentmodelprior))
end

function combine_prior(skymodel, ::NamedTuple{()})
    return NamedDist((;sky = skymodel,))
end

function combine_prior(::NamedTuple{()}, instrumentmodel)
    return NamedDist((; instrument = instrumentmodel))
end

function combine_prior(::NamedTuple{()}, ::NamedTuple{()})
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
    return println(io, map(x -> split(string(datumtype(x)), "{")[1], post.data)...)
    # println(io, "  Prior: ", post.prior)
end


"""
    simulate_observation([rng::Random.AbstractRNG], post::VLBIPosterior, θ; add_thermal_noise=true)

Create a simulated observation using the posterior and its data `post` using the parameter
values `θ`. In Bayesian terminology this is a draw from the posterior predictive distribution.

If `add_thermal_noise` is true then baseline based thermal noise is added. Otherwise, we just
return the model visibilities.
"""
function simulate_observation(rng::Random.AbstractRNG, post::VLBIPosterior, θ; add_thermal_noise = true)
    # ls = map(x->likelihood(x, visibilitymap(vlbimodel(post, θ), post.lklhd.ac)), post.lklhd.lklhds)
    vis = forward_model(post, θ)
    ls = map(x -> likelihood(x, vis), post.lklhds)
    if add_thermal_noise
        ms = map(x -> rand(rng, x), ls)
    else
        ms = map(x -> x.μ, ls)
    end
    configs = map(arrayconfig, post.data)
    data = post.data
    return ntuple(length(data)) do i
        di = data[i]
        return EHTObservationTable{datumtype(di)}(ms[i], noise(di), configs[i])
    end
end
simulate_observation(post::VLBIPosterior, θ; add_thermal_noise = true) = simulate_observation(Random.default_rng(), post, θ; add_thermal_noise)

"""
    residuals(post::VLBIPosterior, θ)

Compute the residuals for each data product in `post` using the parameter values `θ`.
The resturn objects are `EHTObservationTables`, where the measurements are the residuals.
"""
function residuals(post::VLBIPosterior, p)
    vis = forward_model(post, p)
    res = map(x -> residual_data(vis, x), post.data)
    return res
end

"""
    chi2(post::AbstractVLBIPosterior, p)

Returns a tuple of the chi-squared values for each data product in the posterior `post` given the parameters `p`.
Note that the chi-square is not reduced.
"""
function chi2(post::AbstractVLBIPosterior, p)
    res = residuals(post, p)
    return map(_chi2, res)
end

function _chi2(res::EHTObservationTable)
    return sum(datatable(res)) do d
        r2 = @. abs2(d.measurement / d.noise)
        # Check if residual is NaN which means that the data is missing
        isnan(r2) && return zero(r2)
        return r2
    end
end

function _chi2(res::EHTObservationTable{<:EHTCoherencyDatum})
    return sum(datatable(res)) do d
        r2 = @. abs2(d.measurement / d.noise)
        r11 = isnan(r2[1, 1]) ? zero(r2[1, 1]) : r2[1, 1]
        r12 = isnan(r2[1, 2]) ? zero(r2[1, 2]) : r2[1, 2]
        r21 = isnan(r2[2, 1]) ? zero(r2[2, 1]) : r2[2, 1]
        r22 = isnan(r2[2, 2]) ? zero(r2[2, 2]) : r2[2, 2]
        return typeof(r2)(r11, r21, r12, r22)
    end
end
