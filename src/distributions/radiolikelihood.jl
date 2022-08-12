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
julia> RadioLikelihood(dcphase1, dlcamp1)
```
"""
struct RadioLikelihood{M,T,A} <: MB.AbstractMeasure
    model::M
    lklhds::T
    ac::A
end


function RadioLikelihood(model, data::EHTObservation...)
    ls = Tuple(map(makelikelihood, data))
    acs = arrayconfig.(data)
    #@argcheck acs[1] == acs[2]
    RadioLikelihood{typeof(model), typeof(ls), typeof(acs[1])}(model, ls, acs[1])
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

function MB.logdensityof(lklhds::MultiRadioLikelihood, m)
    sum(Base.Fix2(logdensityof, m), lklhds.lklhds)
end


function RadioLikelihood(model, data::MT.Likelihood...)
    return RadioLikelihood{typeof(model), typeof(data), Nothing}(model, data, nothing)
end


function Base.show(io::IO, d::RadioLikelihood{T}) where {T}
    println(io, "RadioLikelihood{$T}")
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
