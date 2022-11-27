using MeasureBase: logdensityof, Likelihood
export RadioLikelihood, logdensityof, MultiRadioLikelihood
using LinearAlgebra

"""
    RadioLikelihood(model, data1, data2, ...)

Forms a radio likelihood from a set of data products. These data products must share
the same array data/configuration. If you want to form a likelihood from multiple arrays
such as when fitting different wavelengths or days, you can combine them using
[`MultiRadioLikelihood`](@ref MultiRadioLikelihood)

# Example

```julia-repl
julia> RadioLikelihood(model, dcphase1, dlcamp1)
```
"""
struct RadioLikelihood{M,T,A} <: MB.AbstractMeasure
    model::M
    lklhds::T
    ac::A
end

struct ModelMetadata{M, C}
    model::M
    metadata::C
end

function (m::ModelMetadata)(θ)
    return m.model(θ, m.metadata)
end




function RadioLikelihood(model, data::EHTObservation...)
    ls = Tuple(map(makelikelihood, data))
    acs = arrayconfig.(data)
    #@argcheck acs[1] == acs[2]
    RadioLikelihood{typeof(model), typeof(ls), typeof(acs[1])}(model, ls, acs[1])
end

"""
    RadioLikelihood(model, metadata, obs::EHTObservation...)

Creates a RadioLikelihood using the `model` and its related `metadata`. The `model`
is a function that converts from parameters `θ` to a Comrade
AbstractModel which can be used to compute [`visibilities`](@ref) and a set of
`metadata` that is used by `model` to compute the model.

# Warning

The `model` itself must be a two argument function where the first argument is the set
of model parameters and the second is a container that holds all the additional
information needed to construct the model. An example of this is when the model
needs some precomputed cache to define the model.

# Example
```julia

cache = create_cache(FFTAlg(), IntensityMap(zeros(128,128), μas2rad(100.0), μas2rad(100.0)))

function model(θ, metadata)
    (; r, a) = θ
    m = stretched(ExtendedRing(a), r, r)
    return modelimage(m, metadata.cache)
end

prior = (
         r = Uniform(μas2rad(10.0), μas2rad(40.0)),
         a = Uniform(0.1, 5.0)
         )

RadioLikelihood(model, (cache = cache), obs)
```
"""
function RadioLikelihood(model, metadata::NamedTuple, data::EHTObservation...)
    ls = Tuple(map(makelikelihood, data))
    acs = arrayconfig.(data)
    #@argcheck acs[1] == acs[2]
    mms = ModelMetadata(model, metadata)
    RadioLikelihood{typeof(mms), typeof(ls), typeof(acs[1])}(mms, ls, acs[1])
end

"""
    MultiRadioLikelihood(lklhd1, lklhd2, ...)
Combines multiple likelihoods into one object that is useful for fitting multiple days/frequencies.

```julia-repl
julia> lklhd1 = RadioLikelihood(dcphase1, dlcamp1)
julia> lklhd2 = RadioLikelihood(dcphase2, dlcamp2)
julia> MultiRadioLikelihood(lklhd1, lklhd2)
```

"""
struct MultiRadioLikelihood{L} <: MB.AbstractMeasure
    lklhds::L
end

MultiRadioLikelihood(lklhds::RadioLikelihood...) = MultiRadioLikelihood(lklhds)

function Base.show(io::IO, d::MultiRadioLikelihood)
    println(io, "MultiRadioLikelihood: ")
    println(io, "  Number of obs: ", length(d.lklhds))
    for l in d.lklhds
        print(io, "\t", l)
    end
end

function MB.logdensityof(lklhds::MultiRadioLikelihood, m)
    sum(Base.Fix2(logdensityof, m), lklhds.lklhds)
end


function RadioLikelihood(model, data::MT.Likelihood...)
    return RadioLikelihood{typeof(model), typeof(data), Nothing}(model, data, nothing)
end


function Base.show(io::IO, d::RadioLikelihood{T}) where {T}
    println(io, "RadioLikelihood")
    println(io, "\tNumber of data products: ", length(d.lklhds))
end


"""
    logclosure_amplitudes(vis::AbstractArray, ac::ArrayConfiguration)

Compute the log-closure amplitudes for a set of visibilities and an array configuration

# Notes
This uses a closure design matrix for the computation.
"""
function logclosure_amplitudes(vis::AbstractArray{<:Complex}, ac::ArrayConfiguration)
    lva = log.(abs.(vis))
    return ac.designmat*lva
end

"""
    closure_phases(vis::AbstractArray, ac::ArrayConfiguration)

Compute the closure phases for a set of visibilities and an array configuration

# Notes
This uses a closure design matrix for the computation.
"""
function closure_phases(vis::AbstractArray{<:Complex}, ac::ArrayConfiguration)
    ph = angle.(vis)
    return ac.designmat*ph
end

amplitudes(vis::AbstractArray{<:Complex}, ac::ArrayConfiguration) = abs.(vis)
phase(vis::AbstractArray{<:Complex}) = angle.(vis)



function MB.logdensityof(d::RadioLikelihood, θ::NamedTuple)
    ac = d.ac
    m = d.model(θ)
    vis = visibilities(m, ac)
    return _logdensityofvis(d, vis)
end

function _logdensityofvis(d::RadioLikelihood, vis::AbstractArray)
    # We use a for loop here since Zygote plays nice with this
    acc = logdensityof(first(d.lklhds), vis)
    @inbounds for l in d.lklhds[begin+1:end]
        acc += logdensityof(l, vis)
    end
    return acc
end

#function MB.logdensityof(d::RadioLikelihood, m::ComradeBase.AbstractModel)
#    acc = logdensityof(first(d.lklhds), m)
#    @inbounds for i in 2:length(d.lklhds)
#        acc += logdensityof(d.lklhds[i], m)
#    end
#    return acc
#end

# internal function that creates the likelihood for a set of complex visibilities
function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTVisibilityDatum})
    errinv = inv.(data[:error])
    vis = StructArray{Complex{eltype(data[:visr])}}((data[:visr],data[:visi]))
    ℓ = Likelihood(vis) do (μ, )
        ComplexNormal{(:μ, :τ)}(μ, errinv)
    end
    return ℓ
end


# internal function that creates the likelihood for a set of visibility amplitudes
function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTVisibilityAmplitudeDatum})
    τ = inv.(data[:error])
    amp = getdata(data, :amp)
    ℓ = Likelihood(amp) do (μ, )
        AmpNormal{(:μ, :τ)}(abs.(μ), τ)
    end
    return ℓ
end


# internal function that creates the likelihood for a set of log closure amplitudes
function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTLogClosureAmplitudeDatum})
    dmat = data.config.designmat
    τ  = inv.(data[:error])
    #covm = σvis*dmat*transpose(σvis)

    #C = Cholesky(Symmetric(covm))
    #σ = C.L

    f = Base.Fix2(logclosure_amplitudes, data.config)
    amp = data[:amp]
    ℓ = Likelihood(amp) do (μ,)
        AmpNormal{(:μ, :τ)}(f(μ), τ)
    end
    #return AmpNormal{(:μ, :τ)}(amp, τ)
    return ℓ
end


# internal function that creates the likelihood for a set of closure phase datum
function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTClosurePhaseDatum})
    τ = inv.(getdata(data, :error)).^2
    f = Base.Fix2(closure_phases, data.config)
    phase = data[:phase]
    ℓ = Likelihood(phase) do (μ,)
        CPVonMises{(:μ, :κ)}(f(μ), τ)
    end

    return ℓ
end
