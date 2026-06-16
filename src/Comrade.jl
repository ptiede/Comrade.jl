module Comrade

using Adapt
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
import ProbabilityTransports as PT
import ProbabilityTransports: transport_node, pfwd_step, pback_step!,
    pback_eltype, latent_pfwd, latent_pback, latent_pfwd_and_logdensity,
    transport_to, TVFlat, StdUniform, TransportedDistribution
# `asflat`/`ascube`/`transform`/`inverse` come from VLBIImagePriors' compat shim
# (which delegates to ProbabilityTransports). Import them so Comrade *extends* the
# same bindings on its posterior types (rather than creating shadowing copies that
# would clash with VLBIImagePriors under `using Comrade, VLBIImagePriors`).
import TransformVariables: transform, inverse, dimension
import VLBIImagePriors: asflat, ascube
#using MappedArrays: mappedarray
using Measurements
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
using ReactantCore
using ComradeBase: rgetindex, rsetindex!

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
include("macro_utils.jl")
include("observations/observations.jl")
include("instrument/instrument.jl")
include("skymodels/abstract.jl")
include("posterior/abstract.jl")
include("inference/inference.jl")
include("visualizations/visualizations.jl")
include("dirty_image.jl")
include("mrf_image.jl")
include("rules.jl")

end
