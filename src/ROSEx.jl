"""
    ROSEx
Radio Observation Sampling Exploration
"""
module ROSEx

using Accessors: @set
using DocStringExtensions
using FFTW: fft, fftfreq, fftshift, ifft, ifft!, ifftshift, plan_fft
using PaddedViews
using MappedArrays: mappedarray
using MeasureBase
using LoopVectorization: @turbo
using NFFT: nfft, plan_nfft
using SpecialFunctions
using StructArrays
using Interpolations: interpolate, scale, extrapolate, BSpline, Cubic, Line, OnGrid
#using ImageFiltering: imfilter, imfilter!, Kernel.gaussian, Fill
using PyCall: pyimport, PyNULL, PyObject
# Write your package code here.

const ehtim = PyNULL()

function load_ehtim()
    try
        copy!(ehtim, pyimport("ehtim"))
    catch
        @warn "No ehtim installation found in python path. Some data functionality will not work"
    end
end


export rad2μas, μas2rad, ehtim, load_ehtim

@inline rad2μas(x) = 180/π*3600*1e6*x
@inline μas2rad(x) = x/(180/π*3600*1e6)

include("stokesimage.jl")
include("observations/observations.jl")
include("models/models.jl")
include("likelihoods/likelihoods.jl")

function __init__()
    # try
    #     copy!(ehtim, pyimport("ehtim"))
    # catch
    #     @warn "No ehtim installation found in python path. Some data functionality will not work"
    # end
end

end
