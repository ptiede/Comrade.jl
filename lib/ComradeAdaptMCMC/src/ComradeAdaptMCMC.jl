module ComradeAdaptMCMC

using AdaptiveMCMC
using Comrade
using TypedTables
using AbstractMCMC

export AdaptMCMC

"""
    AdaptMCMC(;ntemp,
              swap=:nonrev,
              algorithm = :ram,
              fulladapt = true,
              acc_sw = 0.234,
              all_levels = false
              )

Create an `AdaptMCMC.jl` sampler. This sampler uses the [`AdaptiveMCMC.jl`](https://github.com/mvihola/AdaptiveMCMC.jl)
package to sample from the posterior. Namely, this is a parallel tempering algorithm with an adaptive
exploration and tempering sampler. For more information please see [https://github.com/mvihola/AdaptiveMCMC.jl].

The arguments of the function are:
 - `ntemp`: Number of temperature to run in parallel tempering
 - `swap`: Which temperature swapping strategy to use, options are:
   - `:norev` (default) uses a non-reversible tempering scheme (still ergodic)
   - `:single` single randomly picked swap
   - `:randperm` swap in random order
   - `:sweep` upward or downward sweeps picked at random
 - `algorithm`: exploration MCMC algorithm (default is :ram which uses robust adaptive metropolis-hastings) options are:
   - `:ram` (default) Robust adaptive metropolis
   - `:am` Adaptive metropolis
   - `:asm` Adaptive scaling metropolis
   - `:aswam` Adaptive scaling within adaptive metropolis
 - `fulladapt`: whether we adapt both the tempering ladder and the exploration kernel (default is true, i.e. adapt everything)
 - `acc_sw`: The target acceptance rate for temperature swaps
 - `all_levels`: Store all tempering levels to memory (warning this can use a lot of memory)
"""
Base.@kwdef struct AdaptMCMC
    ntemp::Int
    swap::Symbol = :nonrev
    algorithm::Symbol = :ram
    fulladapt::Bool = true
    acc_sw::Float64 = 0.234
    all_levels::Bool = false
end

Comrade.samplertype(::Type{<:AdaptMCMC}) = Comrade.IsCube()

"""
    sample(post::Posterior, sampler::AdaptMCMC, nsamples, burnin=nsamples÷2, args...; init_params=nothing, kwargs...)

Sample the posterior `post` using the `AdaptMCMC` sampler. This will produce `nsamples`
with the first `burnin` steps removed. The `init_params` indicate where to start the sampler from
and it is expected to be a `NamedTuple` of parameters.

Possible additional kwargs are:

- `thin::Int = 1`: which says to save only every `thin` sample to memory
- `rng`: Specify a random number generator (default uses GLOBAL_RNG)

This return a tuple where:
 - First element are the chains from the sampler. If `all_levels=false` the only the unit temperature (posterior) chain is returned
 - Second element is the additional ancilliary information about the samples including the
   loglikelihood `logl`, sampler state `state`, average exploration kernel acceptance rate
   `accexp` for each tempering level, and average temperate swap acceptance rates `accswp`
    for each tempering level.
"""
function AbstractMCMC.sample(post::Comrade.TransformedPosterior, sampler::AdaptMCMC, nsamples, burnin=nsamples÷2, args...; init_params=nothing, kwargs...)
    ℓ = logdensityof(post)
    function lpr(xx)
        for x in xx
            (x > 1.0 || x < 0.0) && return -Inf
        end
        return 0.0
    end

    θ0 = init_params
    if isnothing(init_params)
        @warn "No starting location chosen, picking start from random"
        θ0 = prior_sample(post)
    end

    apt = adaptive_rwm(θ0, ℓ, nsamples;
                       algorithm = sampler.algorithm,
                       b = burnin,
                       fulladapt=sampler.fulladapt,
                       L = sampler.ntemp,
                       acc_sw = sampler.acc_sw,
                       log_pr = lpr,
                       all_levels=sampler.all_levels,
                       swaps = sampler.swap,
                       kwargs...
                       )

    stats = (logl = apt.D, state = apt.R, accexp = apt.accRWM, accswp=apt.accSW)
    if sampler.all_levels
        chains = Tuple(Table(transform.(Ref(post), eachcol(apt.allX[i]))) for i in eachindex(apt.allX))
    else
        chains = transform.(Ref(post), eachcol(apt.X)) |> Table
    end
    return chains, stats
end


end
