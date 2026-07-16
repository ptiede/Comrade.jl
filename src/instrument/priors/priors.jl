export ArrayPrior

abstract type AbstractInstrumentPrior <: Distributions.ContinuousMultivariateDistribution end

include("segmentation.jl")
include("independent.jl")
include("array_priors.jl")
include("refant.jl")
include("processes.jl")
include("markov_chain.jl")
include("chained.jl")
include("gauss_markov.jl")
