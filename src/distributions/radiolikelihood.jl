using MeasureTheory
export RadioLikelihood, logdensity
using LinearAlgebra

struct RadioLikelihood{T} <: MeasureBase.AbstractMeasure
    lklhds::T
    function RadioLikelihood(data...)
        ls = map(makelikelihood, data)
        new{typeof(ls)}(ls)
    end
end

function Base.show(io::IO, d::RadioLikelihood{T}) where {T}
    println(io, "RadioLikelihood{$T}")
    println(io, "\tNumber of data products: ", length(d.lklhds))
end


@inline MeasureBase.logdensity(d::RadioLikelihood, m) = sum(x->logdensity(x, m), d.lklhds)


function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTVisibilityAmplitudeDatum})
    u,v = getdata(data, :u), getdata(data, :v)
    τ = inv.(getdata(data, :error))
    f(m) = amplitude.(Ref(m), u, v)

    k = kernel(AmpNormal{(:μ, :τ)},
                   τ = x->τ,
                   μ = x->f(x)
             )
    amp = getdata(data, :amp)
    return Likelihood(k, amp)
end


function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTLogClosureAmplitudeDatum})
    u1, v1 = getdata(data, :u1), getdata(data, :v1)
    u2, v2 = getdata(data, :u2), getdata(data, :v2)
    u3, v3 = getdata(data, :u3), getdata(data, :v3)
    u4, v4 = getdata(data, :u4), getdata(data, :v4)
    τ = inv.(getdata(data, :error))
    f(m) = logclosure_amplitudes(m, u1, v1, u2, v2, u3, v3, u4, v4)
    k = kernel(AmpNormal{(:μ, :τ)},
                τ = x->τ,
                μ = x->f(x)
              )

    amp = getdata(data, :amp)
    return Likelihood(k, amp)
end

function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTClosurePhaseDatum})
    u1, v1 = getdata(data, :u1), getdata(data, :v1)
    u2, v2 = getdata(data, :u2), getdata(data, :v2)
    u3, v3 = getdata(data, :u3), getdata(data, :v3)
    τ = inv.(getdata(data, :error)).^2
    f(m) = closure_phases(m, u1, v1, u2, v2, u3, v3)
    d = CPVonMises{(:μ, :κ)}
    k = kernel(d,
                κ = x->τ,
                μ = x->f(x)
              )

    amp = getdata(data, :phase)
    return Likelihood(k, amp)
end
