module ComradeAHMC

using AbstractMCMC
using Reexport
@reexport using AdvancedHMC
using Comrade
using DocStringExtensions
using LogDensityProblems, LogDensityProblemsAD
using ArgCheck: @argcheck
using Random
using Accessors
using Printf
using StatsBase
using AbstractMCMC: Sample
using Serialization


function __init__()
    return @warn "ComradeDynesty is deprecated, AdvancedHMC.jl is now an extension."
end


end
