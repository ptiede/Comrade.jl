module ComradeNested

using Comrade
using NestedSamplers
using TupleVectors

Comrade.samplertype(::Type{<:Nested}) = IsCube()

function AbstractMCMC.sample(post::Comrade.TransformedPosterior, sampler::Nested, args...; kwargs...)
    ℓ(x) = logdensityof(post, x)
    model = NestedModel(ℓ, identity)

    samples, stats = sample(model, sampler, args...; chain_type=Array, kwargs...)
    weights = samples[:, end]
    chain = transform.(Ref(post), eachrow(samples[:,1:end-1]))
    return TupleVector(chain), merge((;weights,), stats)
end

end
