"""
    R.O.B.S
Radio Observation Bayesian
"""
module ROSE

using DocStringExtensions
using FFTW
using LoopVectorization
using Memoize
using SpecialFunctions
using StructArrays
# Write your package code here.

export Disk, Gaussian,
       intensity, visibility, flux,
       RImage, SqExpKernel, BSplineKernel,
       scale, shift, rotate,
       stokesmatrix, stokesmatrix!

include("models/models.jl")
end
