"""
    AbstractSkyModel

The abstract type for Comrade Sky Models. For a concrete implementation see [`SkyModel`](@ref).

Any subtype must implement the following methods
 - `ObservedSkyModel(m::AbstractSkyModel, array::AbstractArrayConfiguration)`: Constructs an observed sky model
    given the sky model `m` and the array configuration `array`. This method is used to compute the visibilities
    and image of the sky model.

The following methods have default implementations:
 - `idealvisibilities(m::AbstractSkyModel, x)`: Computes the ideal visibilities of the sky model `m`
    given the model parameters `x`.
 - `skymodel(m::AbstractSkyModel, x)`: Returns the sky model image given the model parameters `x`.
 - `domain(m::AbstractSkyModel)`: Returns the domain of the sky model `m`.
 - `set_array(m::AbstractSkyModel, array::AbstractArrayConfiguration)`: Sets the array configuration
    for the sky model `m` and returns the observed sky model and prior.
 - `set_prior(m::AbstractSkyModel, array::AbstractArrayConfiguration)`: Sets the prior for the sky model
    `m` given the array configuration `array`. This is used to set the priors for the model parameters.

"""
abstract type AbstractSkyModel end

function Base.show(io::IO, mime::MIME"text/plain", m::AbstractSkyModel)
    T = typeof(m)
    ST = split(split(" $T", '{')[1], ".")[end]
    printstyled(io, ST; bold = true, color = :blue)
    println(io)
    println(io, "  with map: $(skymodel(m))")
    # GT = typeof(domain(m))
    # SGT = split("$GT", '{')[1]
    print(io, "   on grid: \n")
    show(io, mime, domain(m))
    return print(io, "\n   )\n")
end

skymodel(m::AbstractSkyModel) = getfield(m, :f)


function set_array(m::AbstractSkyModel, array::AbstractArrayConfiguration)
    return ObservedSkyModel(m, array), NamedDist(set_prior(m, array))
end

function domain(m::AbstractSkyModel; kwargs...)
    return getfield(m, :grid)
end

"""
    idealvisibilities(m::AbstractSkyModel, x)

Computes the ideal non-corrupted visibilities of the sky model `m` given the model parameters `x`.
"""
function idealvisibilities(m::AbstractSkyModel, x)
    skym = skymodel(m, x.sky)
    return visibilitymap(skym, domain(m))
end

function skymodel(m::AbstractSkyModel, x)
    return m.f(x, m.metadata)
end

function set_prior(m::AbstractSkyModel, array::AbstractArrayConfiguration)
    return getfield(m, :prior)
end

struct ObservedSkyModel{F, G <: VLBISkyModels.AbstractDomain, M} <: AbstractSkyModel
    f::F
    grid::G
    metadata::M
end


include("models.jl")
include("fixed.jl")
include("multi.jl")
