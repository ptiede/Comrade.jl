"""
    Comrade
Radio Observation Sampling Exploration
"""
module Comrade

using AbstractFFTs
using AbstractMCMC
using Accessors: @set
using ArgCheck: @argcheck
using BasicInterpolators
using BlockDiagonals
using DocStringExtensions
using ChainRulesCore
using ComradeBase
using ForwardDiff
using FFTW: fft, fftfreq, fftshift, ifft, ifft!, ifftshift, plan_fft
using FLoops
using MappedArrays: mappedarray
using MeasureBase, MeasureTheory
using NFFT: nfft, plan_nfft
using PaddedViews
using PyCall: pyimport, PyNULL, PyObject
using SpecialFunctions
using Reexport
using Requires: @require
using StructArrays: StructVector, StructArray, append!!
using Tables
# Write your package code here.

@reexport using ComradeBase
using ComradeBase: visanalytic, imanalytic

export SequentialEx, ThreadedEx

const ehtim = PyNULL()

"""
    $(SIGNATURES)
Loads the [eht-imaging](https://github.com/achael/eht-imaging) library and stores it in the
`ehtim` variable.

# Notes
This will fail if ehtim isn't installed in the python installation that PyCall references.
"""
function load_ehtim()
    #try
    copy!(ehtim, pyimport("ehtim"))
    #catch
    #    @warn "No ehtim installation found in python path. Some data functionality will not work"
    #end
end


export rad2μas, μas2rad, ehtim, load_ehtim

"""
    rad2μas(x)
Converts a number from radians to μas
"""
@inline rad2μas(x) = 180/π*3600*1e6*x

"""
    μas2rad(x)
Converts a number from μas to rad
"""
@inline μas2rad(x) = x/(180/π*3600*1e6)


#include("interface.jl")
#include("images/images.jl")
import ComradeBase: flux, radialextent, intensitymap, intensitymap!
export create_cache
include("observations/observations.jl")
include("models/models.jl")
include("distributions/distributions.jl")
include("distributions/radiolikelihood.jl")
include("visualizations/visualizations.jl")
include("bayes/bayes.jl")
include("inference/inference.jl")
include("calibration/gains.jl")

function __init__()
    @require AdvancedHMC="0bf59076-c3b1-5ca4-86bd-e02cd72cde3d" include("inference/advancedhmc.jl")
    @require AdaptiveMCMC="717c3277-546c-407d-8270-09a39a0919a0" include("inference/adaptivemcmc.jl")
    @require NestedSamplers="41ceaf6f-1696-4a54-9b49-2e7a9ec3782e" include("inference/nested.jl")
    @require Dynesty="eb527566-0f3e-4aab-bb5f-9d2e403dba70" include("inference/dynesty.jl")
    @require Pathfinder="b1d3bc72-d0e7-4279-b92f-7fa5d6d2d454" include("inference/pathfinder.jl")
    @require GalacticOptim="a75be94c-b780-496d-a8a9-0878b188d577" include("inference/galacticoptim.jl")
end

end
