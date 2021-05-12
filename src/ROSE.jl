"""
    R.O.S.E.
Radio Observation Sampling Extravagance
"""
module ROSE

using DocStringExtensions
using FFTW: fft, fft
using PaddedViews
using Memoization
using LoopVectorization: @avx
using SpecialFunctions
using StructArrays
using Interpolations: interpolate, scale, extrapolate, BSpline, Cubic, Line, OnGrid
using ImageFiltering: imfilter, Kernel.gaussian, Fill, Algorithm.FFT
# Write your package code here.

export Disk, Gaussian, ConcordanceCrescent, MRing,
       intensity, visibility, flux,
       RImage, SqExpKernel, BSplineKernel,
       stretched, shifted, rotated, smoothed, renormed,
       load_tpy,
       getdata, renorm, imagepixels,
       intensitymap, intensitymap!,
       CPNormal

include("stokesimage.jl")
include("observations/observations.jl")
include("models/models.jl")
include("likelihoods/likelihoods.jl")
end
