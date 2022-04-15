using MeasureTheory
export RadioLikelihood, logdensity
using LinearAlgebra

struct RadioLikelihood{T,A} <: MeasureBase.AbstractMeasure
    lklhds::T
    ac::A
    function RadioLikelihood(data...)
        ls = Tuple(map(makelikelihood, data))
        ac = arrayconfig(first(data))
        new{typeof(ls), typeof(ac)}(ls, ac)
    end
end

function Base.show(io::IO, d::RadioLikelihood{T}) where {T}
    println(io, "RadioLikelihood{$T}")
    println(io, "\tNumber of data products: ", length(d.lklhds))
end


amplitudes(vis::AbstractArray{<:Complex}) = abs.(vis)
phase(vis::AbstractArray{<:Complex}) = angle.(vis)



function MeasureBase.logdensity(d::RadioLikelihood, m)
    ac = d.ac
    vis = visibilities(m, ac.ac.data.u, ac.ac.data.v)
    return MeasureBase.logdensity(d, vis)
end

function MeasureBase.logdensity(d::RadioLikelihood, vis::AbstractArray)
    # We use a for loop here since Zygote plays nice with this
    acc = logdensity(first(d.lklhds), vis)
    @inbounds for i in 2:length(d.lklhds)
        acc += logdensity(d.lklhds[i], vis)
    end
    return acc
end

function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTVisibilityDatum})
    errinv = inv.(data[:error])
    τ(_) = errinv
    μ(_) = x
    k = kernel(ComplexNormal{(:μ, :τ)},
                τ = τ,
                μ = μ
            )
    vis = StructArray{Complex{eltype(data[:visr])}}((data[:visr],data[:visi]))
    return Likelihood(k, vis)
end


function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTVisibilityAmplitudeDatum})
    τ = inv.(getdata(data, :error))
    τf(_) = τ
    f(v) = abs.(v)

    k = kernel(AmpNormal{(:μ, :τ)},
                   τ = τf,
                   μ = f
             )
    amp = getdata(data, :amp)
    return Likelihood(k, amp)
end


function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTLogClosureAmplitudeDatum})
    dmat = data.config.designmat
    #σ   = diag(data[:error])
    #covm = σ*dmat*transpose(σ)
    τ = inv.(data[:error])
    τf(x) = τ
    f(vis) = dmat*log.(abs.(vis))
    k = kernel(AmpNormal{(:μ, :τ)},
                vis->(τ = τ,
                      μ = f(vis))
              )

    amp = data[:amp]
    #return AmpNormal{(:μ, :τ)}(amp, τ)
    return Likelihood(k, amp)
end

function makelikelihood(data::Comrade.EHTObservation{<:Real, <:Comrade.EHTClosurePhaseDatum})
    τ = inv.(getdata(data, :error)).^2
    τf(_) = τ
    dmat = data.config.designmat
    f(vis) = dmat*angle.(vis)
    d = CPVonMises{(:μ, :κ)}
    k = kernel(d,
                vis->(κ = τ,
                      μ = f(vis))
              )

    phase = data[:phase]
    #return CPVonMises{(:μ, :κ)}(phase, τ)
    return Likelihood(k, phase)
end
