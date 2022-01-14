"""
    Comrade
Radio Observation Sampling Exploration
"""
module Comrade

using AbstractFFTs
using Accessors: @set
using BasicInterpolators
using DocStringExtensions
using ChainRulesCore
using ComradeBase
using ForwardDiff
using FFTW: fft, fftfreq, fftshift, ifft, ifft!, ifftshift, plan_fft
using Interpolations: interpolate, scale, extrapolate, BSpline, Cubic, Line, OnGrid
using LoopVectorization: @turbo
using MappedArrays: mappedarray
using MeasureBase
using NFFT: nfft, plan_nfft
using PaddedViews
using PyCall: pyimport, PyNULL, PyObject
using SpecialFunctions
using Reexport
using StaticArrays: FieldVector, FieldMatrix
using StructArrays: StructArray
using ZygoteRules
# Write your package code here.

@reexport using ComradeBase


const ehtim = PyNULL()

"""
    $(SIGNATURES)
Loads the [eht-imaging](https://github.com/achael/eht-imaging) library and stores it in the
`ehtim` variable.

# Notes
This will fail if ehtim isn't installed in the python installation that PyCall references.
"""
function load_ehtim()
    try
       copy!(ehtim, pyimport("ehtim"))
    catch
        @warn "No ehtim installation found in python path. Some data functionality will not work"
    end
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
import ComradeBase: flux
include("observations/observations.jl")
include("models/models.jl")
include("distributions/distributions.jl")

function __init__()
    # try
    #     copy!(ehtim, pyimport("ehtim"))
    # catch
    #     @warn "No ehtim installation found in python path. Some data functionality will not work"
    # end
end

end
