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



function MeasureBase.logdensityof(d::RadioLikelihood, m)
    ac = d.ac
    vis = visibilities(m, ac)
    return MeasureBase.logdensityof(d, vis)
end

function MeasureBase.logdensityof(d::RadioLikelihood, vis::AbstractArray)
    # We use a for loop here since Zygote plays nice with this
    acc = logdensityof(first(d.lklhds), vis)
    @inbounds for i in 2:length(d.lklhds)
        acc += logdensityof(d.lklhds[i], vis)
    end
    return acc
end

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
