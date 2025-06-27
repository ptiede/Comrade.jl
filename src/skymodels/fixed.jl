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

skymodel(m::FixedSkyModel) = m.model

function FixedSkyModel(m::AbstractModel, grid::AbstractRectiGrid; algorithm = NFFTAlg())
    return FixedSkyModel(m, grid, algorithm)
end

function ObservedSkyModel(m::FixedSkyModel, arr::AbstractArrayConfiguration)
    gfour = FourierDualDomain(m.grid, arr, m.algorithm)
    img = intensitymap(m.model, gfour)
    vis = visibilitymap(m.model, gfour)
    return ObservedSkyModel(m, gfour, (; img, vis))
end

function set_prior(::FixedSkyModel, ::AbstractArrayConfiguration)
    return (;)
end

function idealvisibilities(m::ObservedSkyModel{<:FixedSkyModel}, x)
    return m.metadata.vis
end

function skymodel(m::ObservedSkyModel{<:FixedSkyModel}, x)
    return m.f.model
end
