struct VLBIPosterior{D, T, P, MS<:ObservedSkyModel, MI<:AbstractInstrumentModel} <: AbstractVLBIPosterior
    data::D
    lklhds::T
    prior::P
    skymodel::MS
    instrumentmodel::MI
end

(post::VLBIPosterior)(θ) = logdensityof(post, θ)




"""
    VLBIPosterior(skymodel::SkyModel, instumentmodel::InstrumentModel, dataproducts::EHTObservationTable...)

Creates a VLBILikelihood using the `skymodel` its related metadata `skymeta`
and the `instrumentmodel` and its metadata `instumentmeta`.
. The `model`
is a function that converts from parameters `θ` to a Comrade
AbstractModel which can be used to compute [`visibilitymap`](@ref) and a set of
`metadata` that is used by `model` to compute the model.

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
        skymodel::SkyModel,
        instrumentmodel::AbstractInstrumentModel,
        dataproducts::EHTObservationTable...;
        )

    # This is needed because the prior is causing runtimeActivity
    # warnings in Enzyme
    Enzyme.API.runtimeActivity!(true)


    array = arrayconfig(dataproducts[begin])
    int, intprior = set_array(instrumentmodel, array)
    sky, skyprior = set_array(skymodel, array)

    total_prior = combine_prior(skyprior, intprior)

    ls = Tuple(map(makelikelihood, dataproducts))

    return VLBIPosterior{
                typeof(dataproducts),typeof(ls),typeof(total_prior),
                typeof(sky), typeof(int)}(dataproducts, ls, total_prior, sky, int)
end

VLBIPosterior(skymodel::SkyModel, dataproducts::EHTObservationTable...) =
    VLBIPosterior(skymodel, IdealInstrumentModel(), dataproducts...)

function combine_prior(skyprior, instrumentmodelprior)
    return NamedDist((sky=skyprior, instrument=instrumentmodelprior))
end

function combine_prior(skymodel, ::Tuple{})
    return NamedDist((sky=skymodel,))
end

function combine_prior(::Tuple{}, instrumentmodel)
    return NamedDist((instrument=skymodel.instrument,))
end


"""
    simulate_observation([rng::Random.AbstractRNG], post::VLBIPosterior, θ; add_thermal_noise=true)

Create a simulated observation using the posterior and its data `post` using the parameter
values `θ`. In Bayesian terminology this is a draw from the posterior predictive distribution.

If `add_thermal_noise` is true then baseline based thermal noise is added. Otherwise, we just
return the model visibilities.
"""
function simulate_observation(rng::Random.AbstractRNG, post::VLBIPosterior, θ; add_thermal_noise=true)
    # ls = map(x->likelihood(x, visibilitymap(vlbimodel(post, θ), post.lklhd.ac)), post.lklhd.lklhds)
    vis = forward_model(vlbimodel(post, θ), θ)
    ls = map(x->likelihood(x, vis), post.lklhd.lklhds)
    if add_thermal_noise
        ms = map(x->rand(rng, x), ls)
    else
        ms = map(x->x.μ, ls)
    end
    configs = map(arrayconfig, post.data)
    data = post.data
    return map(eachindex(ms, configs)) do i
        di = data[i]
        return EHTObservationTable{datumtype(di)}(ms[i], noise(di), configs[i])
    end
end
simulate_observation(post::VLBIPosterior, θ) = simulate_observation(Random.default_rng(), post, θ)

function residuals(post::VLBIPosterior, p)
    vis = forward_model(post, p)
    res = map(x->residuals_data(vis, x), post.data)
    return res
end
