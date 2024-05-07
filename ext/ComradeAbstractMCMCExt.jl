module ComradeAbstractMCMCExt

using Comrade
if isdefined(Base, :get_extension)
    using AbstractMCMC
    using ADTypes
    using JLD2
    using LogDensityProblems
    using LogDensityProblemsAD
    using Serialization
    using StructArrays
else
    using ..AbstractMCMC
    using ..ADTypes
    using ..JLD2
    using ..LogDensityProblems
    using ..LogDensityProblemsAD
    using ..Serialization
    using ..StructArrays
end



end
