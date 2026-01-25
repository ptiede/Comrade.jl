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

function VLBISkyModels.FourierDualDomain(grid::AbstractRectiGrid, array::AbstractArrayConfiguration, alg::FourierTransform; executor = executor(grid))
    return FourierDualDomain(grid, domain(array; executor), alg)
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
