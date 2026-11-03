export SkyModel, FixedSkyModel


struct SkyModel{F, P, G <: AbstractDomain, A <: FourierTransform, M} <: AbstractSkyModel
    f::F
    prior::P
    grid::G
    algorithm::A
    metadata::M
end


"""
    SkyModel(f, prior, grid::AbstractRectiGrid; algorithm = NFFTAlg(), metadata=nothing)

Construct a sky model using the function map `f` with parameter priors `prior`, where the image
model is defined on the domain `grid`. If the underlying model produced by `f` is non-analytic,
then `algorithm` is used to numerically Fourier transform the model image. The `metadata` option
contains additional information needed by the model `f`.

# Arguments

 - `f(x, grid, metadata)` : The sky function. `x` are the model parameters, `grid` is the image
   grid, and `metadata` is the `metadata` option. A model that needs no metadata is written
   `f(x, grid, _)`.
- `prior` : A named tuple of priors for the model parameters defined in `x`. Each model parameter
   used in `x` must have a defined element in the prior.
- `grid` : The domain on which the model is defined. This defines the field of view and resolution
   of the model. Note that if `f` produces a analytic model then this field of view isn't used
   directly in the computation of the visibilities.

# Optional Arguments
- `algorithm` : The Fourier transform algorithm used to compute the visibilities. The default is
   `NFFTAlg()` which uses a non-uniform fast Fourier transform. Other options can be found by using
   `subtypes(VLBISkyModels.FourierTransform)`
- `metadata` : Additional information needed by the model `f`, passed as its third argument.
   It can be any object. The grid is not part of it: `f` receives the grid as its second
   argument.
"""
function SkyModel(f, prior, grid::AbstractRectiGrid; algorithm = NFFTAlg(), metadata = nothing)
    return SkyModel(f, prior, grid, algorithm, metadata)
end

# If we are using a analytic model then we don't need to plan the FT and we
# can save some memory by not storing the plans.
struct AnalyticAlg <: FourierTransform end
struct AnalyticPlan <: VLBISkyModels.AbstractPlan end
VLBISkyModels.getplan(::AnalyticPlan) = nothing
VLBISkyModels.getphases(::AnalyticPlan) = nothing
VLBISkyModels.create_plans(::AnalyticAlg, imgdomain, visdomain) = (AnalyticPlan(), AnalyticPlan())

"""
    ObservedSkyModel(sky::AbstractSkyModel, t::VisTable)

Constructs a sky model that has been theoretically observed by an array with configuration `array`.
Note that this is not a public facing method and is used internally to construct the observed sky
model for use in [`VLBIPosterior`](@ref). Users should construct a [`SkyModel`](@ref) and
pass that to a [`VLBIPosterior`](@ref) object instead.

"""
function ObservedSkyModel(m::SkyModel, arr::VisTable)
    x = rand(NamedDist(m.prior))
    ms = m.f(x, m.grid, m.metadata)
    # if analytic don't bother planning the FT
    if ComradeBase.visanalytic(typeof(ms)) === ComradeBase.IsAnalytic()
        g = FourierDualDomain(m.grid, arr, AnalyticAlg())
    else
        g = FourierDualDomain(m.grid, arr, m.algorithm)
    end
    return ObservedSkyModel(m.f, g, m.metadata)
end
