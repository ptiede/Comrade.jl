module ComradeReactantExt
using Comrade
using Comrade: Adapt
using ComradeBase: ComradeBase, DD, ShardLayout, StructuredDomain
using Reactant
using StructArrays
using StaticArraysCore
using Accessors
using PolarizedTypes
using DensityInterface
using LinearAlgebra
using SparseArrays

include("todevice.jl")
include("misc.jl")
include("sample.jl")


end
