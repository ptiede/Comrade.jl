"""
    R.O.S.E.
Radio Observation Sampling Extravagance
"""
module ROSE

using DocStringExtensions
using FFTW
using LoopVectorization
using Memoization
using SpecialFunctions
using StructArrays
using Interpolations
using ImageFiltering: imfilter, Kernel.gaussian, Fill, Algorithm.FFT
# Write your package code here.

export Disk, Gaussian, ConcordanceCrescent,
       intensity, visibility, flux,
       RImage, SqExpKernel, BSplineKernel,
       stretched, shifted, rotated, smoothed, renormed,
       load_tpy,
       getdata, renorm, pixel_iterator,
       stokesmatrix, stokesmatrix!

include("models/models.jl")
include("observations/observations.jl")
end
