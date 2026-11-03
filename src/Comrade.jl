module Comrade

using Adapt
using AbstractMCMC
using DensityInterface
import Distributions as Dists
using DocStringExtensions
using ChainRulesCore
using EnzymeCore
using EnzymeCore: EnzymeRules
using LogDensityProblems
using LinearAlgebra
import ProbabilityTransports as PT
# Only the names used unqualified below; node-protocol methods (`pfwd_step` etc.) are
# always extended `PT.`-qualified.
import ProbabilityTransports: transport_node, latent_pfwd, latent_pback,
    latent_pfwd_and_logdensity, transport_to, TVFlat, StdUniform, TransportedDistribution
# `asflat`/`ascube`/`transform`/`inverse` come from VLBIImagePriors' compat shim
# (which delegates to ProbabilityTransports). Import them so Comrade *extends* the
# same bindings on its posterior types (rather than creating shadowing copies that
# would clash with VLBIImagePriors under `using Comrade, VLBIImagePriors`).
import TransformVariables: transform, inverse, dimension
import TransformVariables as TV
import VLBIImagePriors: asflat, ascube
#using MappedArrays: mappedarray
using Accessors: setproperties
using Printf
using Random
using Reexport
using StaticArraysCore
using Statistics
using StructArrays: StructArray
import StructArrays
using ComradeBase: AbstractDomain, AbstractRectiGrid
using VLBISkyModels: FourierTransform, FourierDualDomain
import DimensionalData as DD

# Reexport the core libraries for Comrade
@reexport using VLBISkyModels
@reexport using ComradeBase
@reexport using PolarizedTypes
@reexport using VLBIImagePriors
@reexport using VLBIInstrumentModels

import VLBIInstrumentModels: observe, domain, simulate_observation, selfcal_phase
using VLBIInstrumentModels: AbstractInstrumentModel, IdealInstrumentModel, FieldInstruments,
    ObservedInstrumentModel, ObservedSelfCalInstrumentModel, instrument,
    _closurevalue, _dense, _elements, _measuredvalues, _polname,
    _observedvalues, _storage, _withmodel,
    _kwname, _scan_used_names, _split_tildes,
    WhitenedHierarchicalTransform, StdWhitenedHierarchicalTransform,
    MarkovChains, chain_term, _flat_dim, _fused_flat, _std_dim

using ComradeBase: AbstractRectiGrid, AbstractDomain, AbstractModel


export rad2μas, μas2rad


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
include("datatypes.jl")
include("skymodels/abstract.jl")
include("posterior/abstract.jl")
include("inference/inference.jl")
include("visualizations/visualizations.jl")
include("mrf_image.jl")
include("rules.jl")

@static if VERSION >= v"1.11"
    eval(
        Expr(
            :public,
            :AbstractVLBIPosterior, :TransformedVLBIPosterior,
            :AbstractSkyModel, :ObservedSkyModel,
            :logprior, :likelihood, :idealmaps, :forward_model_map, :rmap,
        )
    )
end

end
