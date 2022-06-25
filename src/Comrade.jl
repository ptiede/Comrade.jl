"""
    Comrade
Composable Modeling of Radio Emission
"""
module Comrade

using AbstractFFTs
using AbstractMCMC
using Accessors: @set
using ArgCheck: @argcheck
using BasicInterpolators
using DocStringExtensions
using ChainRulesCore
using ComradeBase
using FITSIO
using FileIO
using ForwardDiff
using FFTW: fft, fftfreq, fftshift, ifft, ifft!, ifftshift, plan_fft
using FLoops
#using MappedArrays: mappedarray
import MeasureBase as MB
import MeasureTheory as MT
using NFFT: nfft, plan_nfft
using PaddedViews
using PyCall: pyimport, PyNULL, PyObject
using SpecialFunctions
using Reexport
using Requires: @require
using StructArrays: StructVector, StructArray, append!!
using Tables
using UUIDs
# Write your package code here.

@reexport using ComradeBase
export AD
using ComradeBase: visanalytic, imanalytic

export SequentialEx, ThreadedEx

const ehtim = PyNULL()

"""
    load_ehtim()

Loads the [eht-imaging](https://github.com/achael/eht-imaging) library and stores it in the
exported `ehtim` variable.

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


export rad2μas, μas2rad, ehtim, load_ehtim, logdensity_def, logdensityof

"""
    rad2μas(x)
Converts a number from radians to micro-arcseconds (μas)
"""
@inline rad2μas(x) = 180/π*3600*1e6*x

"""
    μas2rad(x)
Converts a number from micro-arcseconds (μas) to rad
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
    # FIX THIS
    del_format(format"FITS")
    add_format(format"FITS",
        # See https://www.loc.gov/preservation/digital/formats/fdd/fdd000317.shtml#sign
        [0x53,0x49,0x4d,0x50,0x4c,0x45,0x20,0x20,0x3d,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0x54],
        [".fit", ".fits", ".fts", ".FIT", ".FITS", ".FTS", ".fit",],
        [:FITSIO => UUID("525bcba6-941b-5504-bd06-fd0dc1a4d2eb")],
        [:AstroImages => UUID("fe3fc30c-9b16-11e9-1c73-17dabf39f4ad")],
        [:Comrade => UUID("99d987ce-9a1e-4df8-bc0b-1ea019aa547b")]
    )

end

end
