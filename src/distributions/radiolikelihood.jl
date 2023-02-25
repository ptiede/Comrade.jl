export RadioLikelihood, logdensityof, MultiRadioLikelihood, likelihood
using LinearAlgebra
using VLBILikelihoods

abstract type VLBILikelihood end
@inline DensityInterface.DensityKind(::VLBILikelihood) = DensityInterface.IsDensity()


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
struct RadioLikelihood{M,T,A,P,V} <: VLBILikelihood
    model::M
    lklhds::T
    ac::A
    positions::P
    vis::V
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
    positions = getuvtimefreq(data[1].config)
    vis = similar(positions.U, Complex{eltype(positions.U)})
    vis .= 0
    #@argcheck acs[1] == acs[2]
    RadioLikelihood{typeof(model), typeof(ls), typeof(acs[1]), typeof(positions), typeof(vis)}(model, ls, acs[1], positions, vis)
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
    positions = getuvtimefreq(data[1].config)
    #@argcheck acs[1] == acs[2]
    mms = ModelMetadata(model, metadata)
    RadioLikelihood{typeof(mms), typeof(ls), typeof(acs[1]), typeof(positions)}(mms, ls, acs[1], positions)
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
struct MultiRadioLikelihood{L} <: VLBILikelihood
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

function DensityInterface.logdensityof(lklhds::MultiRadioLikelihood, m)
    sum(Base.Fix2(logdensityof, m), lklhds.lklhds)
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




function DensityInterface.logdensityof(d::RadioLikelihood, θ::NamedTuple)
    ac = d.positions
    m = d.model(θ)
    # Convert because of conventions
    vis = visibilitymap!(d.vis, m, ac)
    return _logdensityofvis(d, vis)
end

function _logdensityofvis(d::RadioLikelihood, vis::AbstractArray)
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

struct ConditionedLikelihood{F, O} <: VLBILikelihood
    kernel::F
    obs::O
end

DensityInterface.logdensityof(d::ConditionedLikelihood, μ) = logdensityof(d.kernel(μ), d.obs)

function RadioLikelihood(model, data::ConditionedLikelihood...)
    return RadioLikelihood{typeof(model), typeof(data), Nothing}(model, data, nothing)
end


"""
    likelihood(d::ConditionedLikelihood, μ)

Returns the likelihood of the model, with parameters μ. That is, we return the
distribution of the data given the model parameters μ. This is an actual probability
distribution.
"""
likelihood(d::ConditionedLikelihood, μ) = d.kernel(μ)



# internal function that creates the likelihood for a set of complex visibilities
function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTVisibilityDatum})
    Σ = data[:error].^2
    vis = data[:measurement]
    ℓ = ConditionedLikelihood(vis) do μ
        ComplexVisLikelihood(μ, Σ, 0.0)
    end
    return ℓ
end

function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTCoherencyDatum})
    Σ = map(x->x.^2, data[:error])
    vis = data[:measurement]
    ℓ = ConditionedLikelihood(vis) do μ
        CoherencyLikelihood(μ, Σ, 0.0)
    end
    return ℓ
end

# internal function that creates the likelihood for a set of visibility amplitudes
function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTVisibilityAmplitudeDatum})
    Σ = data[:error].^2
    amp = data[:measurement]
    ℓ = ConditionedLikelihood(amp) do μ
        AmplitudeLikelihood(abs.(μ), Σ, 0.0)
    end
    return ℓ
end

# internal function that creates the likelihood for a set of log closure amplitudes
function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTLogClosureAmplitudeDatum})
    dmat = data.config.designmat
    #amp2 = abs.(data.config.ac.data.measurement)
    #Σlamp = data.config.ac.data.error.^2 ./ amp2

    # Form the closure covariance matrix
    #Σlca = PDMat(Matrix(dmat*Diagonal(Σlamp)*transpose(dmat)))
    Σlca = data[:error].^2
    f = Base.Fix2(logclosure_amplitudes, data.config)
    amp = data[:measurement]
    ℓ = ConditionedLikelihood(amp) do μ
        AmplitudeLikelihood(f(μ), Σlca, 0.0)
    end
    return ℓ
end


# internal function that creates the likelihood for a set of closure phase datum
function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTClosurePhaseDatum})
    dmat = data.config.designmat
    #amp2 = abs2.(data.config.ac.data.measurement)
    #Σphase = data.config.ac.data.error.^2 ./ amp2

    # Form the closure covariance matrix
    #Σcp = PDMat(Matrix(dmat*Diagonal(Σphase)*transpose(dmat)))
    Σcp = data[:error].^2
    f = Base.Fix2(closure_phases, data.config)
    phase = data[:measurement]
    ℓ = ConditionedLikelihood(phase) do μ
        ClosurePhaseLikelihood(f(μ), Σcp, 0.0)
    end

    return ℓ
end
