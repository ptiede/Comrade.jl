export AbstractRadioLikelihood, logdensityof, MultiRadioLikelihood, likelihood
using LinearAlgebra
using VLBILikelihoods

export vlbimodel

abstract type AbstractRadioLikelihood end
@inline DensityInterface.DensityKind(::AbstractRadioLikelihood) = DensityInterface.IsDensity()


struct VLBILikelihood{MS,MI,D,T,A,P} <: AbstractRadioLikelihood
    skymodel::MS
    instrumentmodel::MI
    data::D
    lklhds::T
    ac::A
    positions::P
end

"""
    skymodel(post::AbstractRadioLikelihood, θ)

Returns the sky model or image of a posterior using the parameter values`θ`
"""
function skymodel(lklhd::AbstractRadioLikelihood, θ)
    lklhd.skymodel(θ)
end

"""
    skymodel(lklhd::AbstractRadioLikelihood, θ)

Returns the instrument model of a lklhderior using the parameter values`θ`
"""
function instrumentmodel(lklhd::AbstractRadioLikelihood, θ)
    lklhd.instrumentmodel(θ)
end

"""
    vlbimodel(post::Posterior, θ)

Returns the instrument model and sky model as a [`VLBIModel`](@ref) of a posterior using the parameter values `θ`

"""
function vlbimodel(d::AbstractRadioLikelihood, θ)
    skym = skymodel(d, θ)
    intm = instrumentmodel(d, θ)
    return VLBIModel(intm, skym)
end


"""
    dataproducts(d::AbstractRadioLikelihood)

Returns the data products you are fitting as a tuple. The order of the tuple corresponds
to the order of the `dataproducts` argument in [`AbstractRadioLikelihood`](@ref).
"""
function dataproducts(d::AbstractRadioLikelihood)
    return d.data
end



struct ModelMetadata{M, C}
    model::M
    metadata::C
end

function (m::ModelMetadata)(θ)
    return m.model(θ, m.metadata)
end


function _VLBILikelihood(skymodel, instrumentmodel, data::EHTObservationTable...)
    ls = Tuple(map(makelikelihood, data))
    acs = map(arrayconfig, data)
    positions = getuvtimefreq(data[1].config)
    VLBILikelihood{typeof(skymodel), typeof(instrumentmodel), typeof(data), typeof(ls), typeof(acs[1]), typeof(positions)}(skymodel, instrumentmodel, data, ls, acs[1], positions)
end

"""
    VLBILikelihood(skymodel, instumentmodel, dataproducts::EHTObservation...;
                    skymeta=nothing,
                    instrumentmeta=nothing)

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
cache = create_cache(NFFTAlg(dlcamp), IntensityMap(zeros(128,128), μas2rad(100.0), μas2rad(100.0)))

function skymodel(θ, metadata)
    (; r, a) = θ
    (; cache) = metadata
    m = stretched(ExtendedRing(a), r, r)
    return modelimage(m, metadata.cache)
end

function instrumentmodel(g, metadata)
    (;lg, gp) = g
    (;gcache) = metadata
    jonesStokes(lg.*exp.(1im.*gp), gcache)
end

prior = (
         r = Uniform(μas2rad(10.0), μas2rad(40.0)),
         a = Uniform(0.1, 5.0)
         )

VLBILikelihood(skymodel, instrumentmodel, dataproducts::EHTObservation...;
                 skymeta=(;cache,),
                 instrumentmeta=(;gcache))
```
"""
function VLBILikelihood(
        skymodel,
        instrumentmodel,
        dataproducts::EHTObservationTable...;
        skymeta = nothing,
        instrumentmeta = nothing)

    if !isnothing(skymeta)
        skym = ModelMetadata(skymodel, skymeta)
    else
        skym = skymodel
    end

    if !isnothing(instrumentmeta)
        intm = ModelMetadata(instrumentmodel, instrumentmeta)
    else
        intm = instrumentmodel
    end
    return _VLBILikelihood(skym, intm, dataproducts...)
end

"""
    AbstractRadioLikelihood(skymodel, dataproducts::EHTObservation...; skymeta=nothing)

Forms a radio likelihood from a set of data products using only a sky model.
This intrinsically assumes that the instrument model is not required since it is perfect.
This is useful when fitting closure quantities which are independent of the instrument.

If you want to form a likelihood from multiple arrays
such as when fitting different wavelengths or days, you can combine them using
[`MultiRadioLikelihood`](@ref MultiRadioLikelihood)

# Example

```julia-repl
julia> AbstractRadioLikelihood(skymodel, dcphase, dlcamp)
```
"""
function VLBILikelihood(
    skymodel,
    dataproducts::EHTObservationTable...;
    skymeta=nothing)

    if !isnothing(skymeta)
        skym = ModelMetadata(skymodel, skymeta)
    else
        skym = skymodel
    end

    return _VLBILikelihood(skym, Returns(IdealInstrument()), dataproducts...)
end



"""
    MultiRadioLikelihood(lklhd1, lklhd2, ...)
Combines multiple likelihoods into one object that is useful for fitting multiple days/frequencies.

```julia-repl
julia> lklhd1 = AbstractRadioLikelihood(dcphase1, dlcamp1)
julia> lklhd2 = AbstractRadioLikelihood(dcphase2, dlcamp2)
julia> MultiRadioLikelihood(lklhd1, lklhd2)
```

"""
struct MultiRadioLikelihood{L} <: AbstractRadioLikelihood
    lklhds::L
end

MultiRadioLikelihood(lklhds::AbstractRadioLikelihood...) = MultiRadioLikelihood(lklhds)

function Base.show(io::IO, d::MultiRadioLikelihood)
    println(io, "MultiRadioLikelihood: ")
    println(io, "  Number of obs: ", length(d.lklhds))
    for l in d.lklhds
        print(io, "\t", l)
    end
end

function DensityInterface.logdensityof(lklhds::MultiRadioLikelihood, m)
    sum(Base.Fix2(logdensityof, m), lklhds.lklhds)
end




function Base.show(io::IO, d::VLBILikelihood)
    println(io, "VLBILikelihood")
    println(io, "\tNumber of data products: ", length(d.lklhds))
end


"""
    logclosure_amplitudes(vis::AbstractArray, ac::ArrayConfiguration)

Compute the log-closure amplitudes for a set of visibilities and an array configuration

# Notes
This uses a closure design matrix for the computation.
"""
function logclosure_amplitudes(vis::AbstractArray{<:Complex}, ac::ClosureConfig)
    lva = log.(abs.(vis))
    return designmat(ac)*lva
end

"""
    closure_phases(vis::AbstractArray, ac::ArrayConfiguration)

Compute the closure phases for a set of visibilities and an array configuration

# Notes
This uses a closure design matrix for the computation.
"""
function closure_phases(vis::AbstractArray{<:Complex}, ac::ClosureConfig)
    ph = angle.(vis)
    return designmat(ac)*ph
end

amplitudes(vis::AbstractArray{<:Complex}, ac::AbstractArrayConfiguration) = abs.(vis)
phase(vis::AbstractArray{<:Complex}) = angle.(vis)



function DensityInterface.logdensityof(d::AbstractRadioLikelihood, θ)
    ac = d.positions
    m = vlbimodel(d, θ)
    # Convert because of conventions
    vis = visibilitymap(m, ac)
    return _logdensityofvis(d, vis)
end

function _logdensityofvis(d::AbstractRadioLikelihood, vis::AbstractArray)
    # This sum is fast and more generic than a loop
    #f = Base.Fix2(logdensityof, vis)
    #return sum(f, d.lklhds)
        # We use a for loop here since Zygote plays nice with this
        acc = logdensityof(first(d.lklhds), vis)
        @inbounds for l in d.lklhds[begin+1:end]
            acc += logdensityof(l, vis)
        end
        return acc
end

struct ConditionedLikelihood{F, O} <: AbstractRadioLikelihood
    kernel::F
    obs::O
end

DensityInterface.logdensityof(d::ConditionedLikelihood, μ) = logdensityof(d.kernel(μ), d.obs)

function AbstractRadioLikelihood(model, data::ConditionedLikelihood...)
    return AbstractRadioLikelihood{typeof(model), typeof(data), Nothing}(model, data, nothing)
end


"""
    likelihood(d::ConditionedLikelihood, μ)

Returns the likelihood of the model, with parameters μ. That is, we return the
distribution of the data given the model parameters μ. This is an actual probability
distribution.
"""
likelihood(d::ConditionedLikelihood, μ) = d.kernel(μ)



# internal function that creates the likelihood for a set of complex visibilities
function makelikelihood(data::Comrade.EHTObservationTable{<:Real, <:Comrade.EHTComplexVisibilityDatum})
    Σ = data[:error].^2
    vis = data[:measurement]
    lnorm = VLBILikelihoods.lognorm(ComplexVisLikelihood(vis, Σ))
    ℓ = ConditionedLikelihood(vis) do μ
        ComplexVisLikelihood(μ, Σ, lnorm)
    end
    return ℓ
end

function makelikelihood(data::Comrade.EHTObservationTable{T, <:Comrade.EHTCoherencyDatum}) where {T}
    Σ = deepcopy(map(x->x.^2, data[:noise]))
    vis = data[:measurement]
    # lnorm = VLBILikelihoods.lognorm(CoherencyLikelihood(vis, Σ))
    ℓ = ConditionedLikelihood(vis) do μ
        CoherencyLikelihood(μ, Σ, 0.0)
    end
    return ℓ
end

# internal function that creates the likelihood for a set of visibility amplitudes
function makelikelihood(data::Comrade.EHTObservationTable{<:Real, <:Comrade.EHTVisibilityAmplitudeDatum})
    Σ = data[:noise].^2
    amp = data[:measurement]
    ℓ = ConditionedLikelihood(amp) do μ
        RiceAmplitudeLikelihood(abs.(μ), Σ)
    end
    return ℓ
end

# internal function that creates the likelihood for a set of log closure amplitudes
function makelikelihood(data::Comrade.EHTObservationTable{<:Real, <:Comrade.EHTLogClosureAmplitudeDatum})
    dmat = data.config.designmat
    amp2 = abs2.(data.config.ac.data.measurement)
    Σlamp = data.config.ac.data.error.^2 ./ amp2

    # Form the closure covariance matrix
    Σlca = VLBILikelihoods.CholeskyFactor(sparse(dmat*Diagonal(Σlamp)*transpose(dmat)))
    f = Base.Fix2(logclosure_amplitudes, data.config)
    amp = data[:measurement]
    lnorm = VLBILikelihoods.lognorm(AmplitudeLikelihood(amp, Σlca))
    ℓ = ConditionedLikelihood(amp) do μ
        AmplitudeLikelihood(f(μ), Σlca, lnorm)
    end
    return ℓ
end


# internal function that creates the likelihood for a set of closure phase datum
function makelikelihood(data::Comrade.EHTObservationTable{<:Real, <:Comrade.EHTClosurePhaseDatum})
    dmat = data.config.designmat
    amp2 = abs2.(data.config.ac.data.measurement)
    Σphase = data.config.ac.data.error.^2 ./ amp2

    # Form the closure covariance matrix
    Σcp =  VLBILikelihoods.CholeskyFactor(sparse(dmat*Diagonal(Σphase)*transpose(dmat)))
    f = Base.Fix2(closure_phases, data.config)
    phase = data[:measurement]
    lnorm = VLBILikelihoods.lognorm(ClosurePhaseLikelihood(phase, Σcp))
    ℓ = ConditionedLikelihood(phase) do μ
        ClosurePhaseLikelihood(f(μ), Σcp, lnorm)
    end

    return ℓ
end
