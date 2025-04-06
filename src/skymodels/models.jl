export SkyModel, FixedSkyModel


"""
    AbstractSkyModel

The abstract type for Comrade Sky Models. For a concrete implementation see [`SkyModel`](@ref).

Any subtype must implement the following methods

 - `set_array(m::AbstractSkyModel, array::AbstractArrayConfiguration)`: Sets the array configuration
    for the sky model `m` and returns the observed sky model and prior.

The following methods have default implementations:
 - `idealvisibilities(m::AbstractSkyModel, x)`: Computes the ideal visibilities of the sky model `m`
    given the model parameters `x`.
 - `skymodel(m::AbstractSkyModel, x)`: Returns the sky model image given the model parameters `x`.
 - `domain(m::AbstractSkyModel)`: Returns the domain of the sky model `m`.
"""
abstract type AbstractSkyModel end


struct SkyModel{F, P, G <: AbstractDomain, A <: FourierTransform, M} <: AbstractSkyModel
    f::F
    prior::P
    grid::G
    algorithm::A
    metadata::M
end

function Base.show(io::IO, mime::MIME"text/plain", m::AbstractSkyModel)
    T = typeof(m)
    ST = split(split(" $T", '{')[1], ".")[end]
    printstyled(io, ST; bold = true, color = :blue)
    println(io)
    println(io, "  with map: $(m.f)")
    GT = typeof(m.grid)
    SGT = split("$GT", '{')[1]
    return print(io, "   on grid: $SGT")
end

"""
    SkyModel(f, prior, grid::AbstractRectiGrid; algorithm = NFFTAlg(), metadata=nothing)

Construct a sky model using the function map `f` with parameter priors `prior`, where the image
model is defined on the domain `grid`. If the underlying model produced by `f` is non-analytic,
then `algorithm` is used to numerically Fourier transform the model image. The `metadata` option
contains additional information needed by the model `f`.

# Arguments

 - `f(x, p)` : A function must be two arguments, where `x` are the model parameters and `p` is the metadata.
    typically `x` and `p` are named tuples.
- `prior` : A named tuple of priors for the model parameters defined in `x`. Each model parameter
   used in `x` must have a defined element in the prior.
- `grid` : The domain on which the model is defined. This defines the field of view and resolution
   of the model. Note that if `f` produces a analytic model then this field of view isn't used
   directly in the computation of the visibilities.

# Optional Arguments
- `algorithm` : The Fourier transform algorithm used to compute the visibilities. The default is
   `NFFTAlg()` which uses a non-uniform fast Fourier transform. Other options can be found by using
   `subtypes(VLBISkyModels.FourierTransform)`
- `metadata` : Additional information needed by the model `f`. These are the addtional arguments
   passed to the model function `f`.
"""
function SkyModel(f, prior, grid::AbstractRectiGrid; algorithm = NFFTAlg(), metadata = nothing)
    return SkyModel(f, prior, grid, algorithm, metadata)
end

function VLBISkyModels.FourierDualDomain(grid::AbstractRectiGrid, array::AbstractArrayConfiguration, alg::FourierTransform; executor = Serial())
    return FourierDualDomain(grid, domain(array; executor), alg)
end

struct ObservedSkyModel{F, G <: VLBISkyModels.AbstractDomain, M} <: AbstractSkyModel
    f::F
    grid::G
    metadata::M
end

function domain(m::AbstractSkyModel; kwargs...)
    return getfield(m, :grid)
end

# If we are using a analytic model then we don't need to plan the FT and we
# can save some memory by not storing the plans.
struct AnalyticAlg <: FourierTransform end
struct AnalyticPlan <: VLBISkyModels.AbstractPlan end
VLBISkyModels.getplan(::AnalyticPlan) = nothing
VLBISkyModels.getphases(::AnalyticPlan) = nothing
VLBISkyModels.create_plans(::AnalyticAlg, imgdomain, visdomain) = (AnalyticPlan(), AnalyticPlan())

"""
    ObservedSkyModel(sky::AbstractSkyModel, array::AbstractArrayConfiguration)

Constructs a sky model that has been theoretically observed by an array with configuration `array`.
Note that this is not a public facing method and is used internally to construct the observed sky
model for use in [`VLBIPosterior`](@ref). Users should construct a [`SkyModel`](@ref) and
pass that to a [`VLBIPosterior`](@ref) object instead.

"""
function ObservedSkyModel(m::SkyModel, arr::AbstractArrayConfiguration)
    x = rand(NamedDist(m.prior))
    ms = m.f(x, m.metadata)
    # if analytic don't bother planning the FT
    if ComradeBase.visanalytic(typeof(ms)) === ComradeBase.IsAnalytic()
        g = FourierDualDomain(m.grid, arr, AnalyticAlg())
    else
        g = FourierDualDomain(m.grid, arr, m.algorithm)
    end
    return ObservedSkyModel(m.f, g, m.metadata)
end


function set_array(m::AbstractSkyModel, array::AbstractArrayConfiguration)
    return ObservedSkyModel(m, array), m.prior
end

"""
    idealvisibilities(m::AbstractSkyModel, x)

Computes the ideal non-corrupted visibilities of the sky model `m` given the model parameters `x`.
"""
function idealvisibilities(m::AbstractSkyModel, x)
    @info typeof(x)
    skym = skymodel(m, x.sky)
    return visibilitymap(skym, domain(m))
end

function skymodel(m::AbstractSkyModel, x)
    return m.f(x, m.metadata)
end

"""
    FixedSkyModel(m::AbstractModel, grid::AbstractRectiGrid; algorithm = NFFTAlg())

Construct a sky model that has no free parameters. This is useful for models where the
image structure is known apriori but the instrument model is unknown.

# Arguments

 - `m` : The fixed sky model.
- `grid` : The domain on which the model is defined. This defines the field of view and resolution
   of the model. Note that if `f` produces a analytic model then this field of view isn't used
   directly in the computation of the visibilities.

# Optional Arguments
- `algorithm` : The Fourier transform algorithm used to compute the visibilities. The default is
   `NFFTAlg()` which uses a non-uniform fast Fourier transform. Other options can be found by using
   `subtypes(VLBISkyModels.FourierTransform)`
"""
Base.@kwdef struct FixedSkyModel{M <: AbstractModel, G, A <: FourierTransform} <: AbstractSkyModel
    model::M
    grid::G
    algorithm::A = NFFTAlg()
end

function FixedSkyModel(m::AbstractModel, grid::AbstractRectiGrid; algorithm = NFFTAlg())
    return FixedSkyModel(m, grid, algorithm)
end

function ObservedSkyModel(m::FixedSkyModel, arr::AbstractArrayConfiguration)
    gfour = FourierDualDomain(m.grid, arr, m.algorithm)
    img = intensitymap(m.model, gfour)
    vis = visibilitymap(m.model, gfour)
    return ObservedSkyModel(m, gfour, (; img, vis))
end

function set_array(m::FixedSkyModel, array::AbstractArrayConfiguration)
    return ObservedSkyModel(m, array), NamedTuple()
end

function idealvisibilities(m::ObservedSkyModel{<:FixedSkyModel}, x)
    return m.metadata.vis
end

function skymodel(m::ObservedSkyModel{<:FixedSkyModel}, x)
    return m.f.model
end
