import ComradeBase: AbstractModel, IsPrimitive, NotPrimitive, IsAnalytic, NotAnalytic,
                    IsPolarized, NotPolarized,
                    visanalytic, imanalytic, isprimitive, ispolarized
import ComradeBase: visibility_point,
                    intensitymap, intensitymap!, intensity_point,
                    flux

export visibility, amplitude, closure_phase, logclosure_amplitude, bispectrum,
       visibilities, amplitudes, closure_phases, logclosure_amplitudes, bispectra,
       flux, intensitymap, intensitymap!, PolarizedModel, convolve!


abstract type AbstractModelImage{M} <: ComradeBase.AbstractModel end

# ChainRulesCore.@non_differentiable visanalytic(M)
# ChainRulesCore.@non_differentiable imanalytic(M)
# ChainRulesCore.@non_differentiable isprimitive(M)
#


include(joinpath(@__DIR__, "methods.jl"))
include(joinpath(@__DIR__, "pulse.jl"))
include(joinpath(@__DIR__, "geometric_models.jl"))
include(joinpath(@__DIR__, "modelimage/modelimage.jl"))
include(joinpath(@__DIR__, "modifiers.jl"))
include(joinpath(@__DIR__, "combinators.jl"))
include(joinpath(@__DIR__, "polarized.jl"))
include(joinpath(@__DIR__, "continuous_image.jl"))
include(joinpath(@__DIR__, "test.jl"))
include(joinpath(@__DIR__, "misc.jl"))
include(joinpath(@__DIR__, "threaded.jl"))
