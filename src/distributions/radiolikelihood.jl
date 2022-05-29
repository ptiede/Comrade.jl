using MeasureTheory
export RadioLikelihood, logdensityof, MultiRadioLikelihood
using LinearAlgebra


struct RadioLikelihood{T,A} <: MeasureBase.AbstractMeasure
    lklhds::T
    ac::A
end


struct MultiRadioLikelihood{L} <: MeasureBase.AbstractMeasure
    lklhds::L
end

"""
    `MultiRadioLikelihood(lklhd1, lklhd2, ...)`
Combines multiple likelihoods into one object that is useful for fitting multiple days/frequencies.
"""
MultiRadioLikelihood(lklhds::RadioLikelihood...) = MultiRadioLikelihood(lklhds)

function MeasureBase.logdensityof(lklhds::MultiRadioLikelihood, m)
    sum(x->logdensityof(x, m), lklhds.lklhds)
end

"""
    `RadioLikelihood(data1, data2, ...)`
Forms a radio likelihood from a set of data products. These data products must share
the same array data/configuration. If you want to form a likelihood from multiple arrays
such as when fitting different wavelengths or days, you can combine them using
`MultiRadioLikelihood`

```julia
lklhd1 = RadioLikelihood(dcphase1, dlcamp1)
lklhd2 = RadioLikelihood(dcphase2, dlcamp2)

lklhd = MultiRadioLikelihood(lklhd1, lklhd2)
```
"""
function RadioLikelihood(data::EHTObservation...)
    ls = Tuple(map(makelikelihood, data))
    acs = arrayconfig.(data)
    #@argcheck acs[1] == acs[2]
    RadioLikelihood{typeof(ls), typeof(acs[1])}(ls, acs[1])
end

function RadioLikelihood(data::Likelihood...)
    return RadioLikelihood{typeof(data), Nothing}(data, nothing)
end


function Base.show(io::IO, d::RadioLikelihood{T}) where {T}
    println(io, "RadioLikelihood{$T}")
    println(io, "\tNumber of data products: ", length(d.lklhds))
end


amplitudes(vis::AbstractArray{<:Complex}) = abs.(vis)
phase(vis::AbstractArray{<:Complex}) = angle.(vis)



function MeasureBase.logdensityof(d::RadioLikelihood, m::ComradeBase.AbstractModel)
    ac = d.ac
    vis = visibilities(m, ac)
    return logdensityof(d, vis)
end

function MeasureBase.logdensityof(d::RadioLikelihood, vis::AbstractArray)
    # We use a for loop here since Zygote plays nice with this
    acc = logdensityof(first(d.lklhds), vis)
    @inbounds for l in d.lklhds[begin+1:end]
        acc += logdensityof(l, vis)
    end
    return acc
end

#function MeasureBase.logdensityof(d::RadioLikelihood, m::ComradeBase.AbstractModel)
#    acc = logdensityof(first(d.lklhds), m)
#    @inbounds for i in 2:length(d.lklhds)
#        acc += logdensityof(d.lklhds[i], m)
#    end
#    return acc
#end

function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTVisibilityDatum})
    errinv = inv.(data[:error])
    #k = kernel(ComplexNormal{(:μ, :τ)},
    #            τ = τ,
    #            μ = μ
    #        )
#
    vis = StructArray{Complex{eltype(data[:visr])}}((data[:visr],data[:visi]))
    ℓ = Likelihood(vis) do (μ, )
        ComplexNormal{(:μ, :τ)}(μ, errinv)
    end
    return ℓ
end


function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTVisibilityAmplitudeDatum})
    τ = inv.(data[:error])
    #τf(_) = τ
    #f(v) = abs.(v)

    #k = kernel(AmpNormal{(:μ, :τ)},
    #               τ = τf,
    #               μ = f
    #         )
    amp = getdata(data, :amp)
    ℓ = Likelihood(amp) do (μ, )
        AmpNormal{(:μ, :τ)}(abs.(μ), τ)
    end
    return ℓ
end


function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTLogClosureAmplitudeDatum})
    dmat = data.config.designmat
    τ  = inv.(data[:error])
    #covm = σvis*dmat*transpose(σvis)

    #C = Cholesky(Symmetric(covm))
    #σ = C.L
    f = let dmat=dmat
        vis->dmat*log.(abs.(vis))
    end
    #k = kernel(AmpNormal{(:μ, :τ)},
    #            vis->(τ = τ,
    #                  μ = f(vis))
    #          )

    amp = data[:amp]
    ℓ = Likelihood(amp) do (μ,)
        AmpNormal{(:μ, :τ)}(f(μ), τ)
    end
    #return AmpNormal{(:μ, :τ)}(amp, τ)
    return ℓ
end

function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTClosurePhaseDatum})
    τ = inv.(getdata(data, :error)).^2
    dmat = data.config.designmat
    f(vis) = dmat*angle.(vis)
    #d = CPVonMises{(:μ, :κ)}
    #k = kernel(d,
    #            vis->(κ = τ,
    #                  μ = f(vis))
    #          )
    phase = data[:phase]
    ℓ = Likelihood(phase) do (μ,)
        CPVonMises{(:μ, :κ)}(f(μ), τ)
    end

    return ℓ
end
