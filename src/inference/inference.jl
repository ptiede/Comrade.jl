using AbstractMCMC, PrettyTables

export chain, samplerstats, samplerinfo

"""
    $(TYPEDEF)

Specifies that the sampling algorithm usually expects a hypercube transform
"""
struct IsCube end

"""
    $(TYPEDEF)

Specifies that the sampling algorithm usually expects a uncontrained transform
"""
struct IsFlat end

struct HasDeriv end
struct NoDeriv end

export sample


struct PosteriorSamples{T, C<:AbstractVector{T}, S, M} <: AbstractVector{T}
    chain::C
    stats::S
    metadata::M
    function PosteriorSamples(chain::C, stats::S, metadata::M) where {C,S,M}
        new{eltype(chain), C, S, M}(_convert2structarr(chain), stats, metadata)
    end
end

function PosteriorSamples(chain::StructArray, stats::StructArray, metadata)
    length(chain) != length(stats) && throw(ArgumentError("chain and stats must have the same length"))
    PosteriorSamples(chain, stats, metadata)
end

_convert2structarr(x::Union{<:NamedTuple, <:AbstractVector}) = StructArray(x, unwrap=(T->!(T<:Real)))
_convert2structarr(x::StructArray) = x

function PosteriorSamples(chain::Union{<:NamedTuple, <:AbstractVector}, stats::Union{<:NamedTuple, <:AbstractVector}, metadata)
    return PosteriorSamples(_convert2structarr(chain), StructArray(stats), metadata)
end

function PosteriorSamples(chain::Union{<:NamedTuple, <:AbstractVector}, metadata)
    return PosteriorSamples(_convert2structarr(chain), metadata)
end



function Base.parent(s::PosteriorSamples)
    return chain(s)
end

function Base.size(s::PosteriorSamples)
    return size(parent(s))
end

function Base.getindex(s::PosteriorSamples, i::Int)
    return getindex(chain(s), i)
end

function Base.setindex!(s::PosteriorSamples, v, i::Int)
    return setindex!(chain(s), v, i)
end

function Base.similar(s::PosteriorSamples, ::Type{S}, dims::Dims) where {S}
    isnothing(samplerstats(s)) && return PosteriorSamples(similar(chain(s), S, dims), nothing, samplerinfo(s))
    return PosteriorSamples(similar(chain(s), S, dims), similar(samplerstats(s), dims), samplerinfo(s))
end

Base.IndexStyle(::Type{<:PosteriorSamples{T,C}}) where {T,C} = Base.IndexStyle(C)


function Base.show(io::IO, ::MIME"text/plain", s::PosteriorSamples)
    println(io, "PosteriorSamples")
    println(io, "  $(length(s[1])) parameters")
    println(io, "  sampler used: ", get(samplerinfo(s), :sampler, "unknown"))
    println(io, "Samples: ")
    ct = chain(s)
    pretty_table(io, Tables.columns(ct);
                    #  formatters = (v,i,j)->round(v, digits=3)
                )
end


function Base.propertynames(s::PosteriorSamples)
    return propertynames(chain(s))
end

function Base.getproperty(s::PosteriorSamples, p::Symbol)
    return getproperty(chain(s), p)
end

function chain(s::PosteriorSamples)
    return getfield(s, :chain)
end

function samplerstats(s::PosteriorSamples)
    return getfield(s, :stats)
end

function samplerinfo(s::PosteriorSamples)
    return getfield(s, :metadata)
end

"""
    samplertype(::Type)

Sampler type specifies whether to use a unit hypercube or unconstrained transformation.
"""
samplertype(::Type) = ArgumentError("samplertype not specified")


"""
    sample(post::Posterior, sampler::S, args...; initial_params=nothing, kwargs...)

Sample a posterior `post` using the `sampler`. You can optionally pass the starting location
of the sampler using `initial_params`, otherwise a random draw from the prior will be used.
"""
function AbstractMCMC.sample(rng::Random.AbstractRNG, post::Posterior, sampler::S, args...; initial_params=nothing, kwargs...) where {S}
    θ0 = initial_params
    if isnothing(initial_params)
        θ0 = prior_sample(post)
    end
    return _sample(samplertype(S), rng, post, sampler, args...; initial_params=θ0, kwargs...)
end

function AbstractMCMC.sample(post::Posterior, sampler, args...; initial_params=nothing, kwargs...)
    sample(Random.default_rng(), post, sampler, args...; initial_params, kwargs...)
end

function _sample(::IsFlat, rng, post, sampler, args...; initial_params, kwargs...)
    tpost = asflat(post)
    θ0 = HypercubeTransform.inverse(tpost, initial_params)
    return sample(rng, tpost, sampler, args...; initial_params=θ0, kwargs...)
end

function _sample(::IsCube, rng, post, sampler, args...; initial_params, kwargs...)
    tpost = ascube(post)
    θ0 = HypercubeTransform.inverse(tpost, initial_params)
    return sample(rng, tpost, sampler, args...; initial_params=θ0, kwargs...)
end


# include(joinpath(@__DIR__, "fishermatrix.jl"))
