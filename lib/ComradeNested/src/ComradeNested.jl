module ComradeNested

using Comrade
using AbstractMCMC
using TupleVectors
using Reexport

@reexport using NestedSamplers

Comrade.samplertype(::Type{<:Nested}) = Comrade.IsCube()

function AbstractMCMC.sample(post::Comrade.TransformedPosterior, sampler::Nested, args...; kwargs...)
    ℓ = logdensityof(post)
    model = NestedModel(ℓ, identity)

    samples, stats = sample(model, sampler, args...; chain_type=Array, kwargs...)
    weights = samples[:, end]
    chain = transform.(Ref(post), eachrow(samples[:,1:end-1]))
    return TupleVector(chain), merge((;weights,), stats)
end

end
