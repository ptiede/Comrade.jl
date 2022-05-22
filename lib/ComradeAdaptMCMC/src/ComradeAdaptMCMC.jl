module ComradeAdaptMCMC

using AdaptiveMCMC
using Comrade
using TupleVectors
using AbstractMCMC

export AdaptMCMC

Base.@kwdef struct AdaptMCMC
    ntemp::Int
    swap::Symbol = :nonrev
    algorithm::Symbol = :ram
    fulladapt::Bool = true
    acc_sw::Float64 = 0.234
    all_levels::Bool = false
end

Comrade.samplertype(::Type{<:AdaptMCMC}) = Comrade.IsCube()

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
        θ0 = transform(tpost, rand(tpost.prior))
    end

    apt = adaptive_rwm(Comrade.HypercubeTransform.inverse(post, θ0), ℓ, nsamples;
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

    chain = transform.(Ref(post), eachcol(apt.X)) |> TupleVector
    stats = (logl = apt.D, state = apt.R, accexp = apt.accRWM, accswp=apt.accSW)
    return chain, stats
end


end
