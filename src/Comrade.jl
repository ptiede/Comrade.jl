"""
    Comrade
Composable Modeling of Radio Emission
"""
module Comrade

using AbstractMCMC
using Accessors: @set
using ArgCheck: @argcheck
using DensityInterface
import Distributions as Dists
using DocStringExtensions
using ChainRulesCore
using EnzymeCore
using EnzymeCore: EnzymeRules
using FillArrays: Fill
using IntervalSets
using LogDensityProblems
using LinearAlgebra
import HypercubeTransform: ascube, asflat, NamedDist, NamedDist, transform, inverse
using HypercubeTransform
#using MappedArrays: mappedarray
using Measurements
using NamedTupleTools
using Printf
using Random
using RecipesBase
using Reexport
using SparseArrays
using StaticArraysCore
using StructArrays: StructVector, StructArray, append!!
import StructArrays
using Tables
import TransformVariables as TV
using ComradeBase: AbstractDomain, AbstractSingleDomain, AbstractRectiGrid
using VLBISkyModels: FourierTransform, FourierDualDomain

# Reexport the core libraries for Comrade
@reexport using VLBISkyModels
@reexport using ComradeBase
@reexport using PolarizedTypes
@reexport using VLBIImagePriors

export linearpol, mbreve, evpa
using ComradeBase: AbstractRectiGrid, AbstractDomain, UnstructuredDomain,
    AbstractModel, AbstractPolarizedModel, AbstractHeader


export rad2μas, μas2rad, logdensity_def, logdensityof


import ComradeBase: flux, radialextent, intensitymap, intensitymap!,
    intensitymap_analytic, intensitymap_analytic!,
    intensitymap_numeric, intensitymap_numeric!,
    visibilitymap, visibilitymap!,
    _visibilitymap, _visibilitymap!,
    visibilitymap_analytic, visibilitymap_analytic!,
    visibilitymap_numeric, visibilitymap_numeric!,
    visanalytic, imanalytic, ispolarized,
    NotAnalytic, IsAnalytic, NotPolarized, IsPolarized,
    visibility_point, intensity_point,
    closure_phase, closure_phasemap,
    logclosure_amplitude, logclosure_amplitudemap,
    visibility, amplitude,
    amplitudemap
include("observations/observations.jl")
include("instrument/instrument.jl")
include("skymodels/models.jl")
include("posterior/abstract.jl")
include("inference/inference.jl")
include("visualizations/visualizations.jl")
include("dirty_image.jl")
include("mrf_image.jl")
include("rules.jl")

end
