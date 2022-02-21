using .Pathfinder

"""
    `pathfinder(post::Posterior, ndraws::Int; init_params=nothing, kwargs...)`

Estimates the posterior using the path finder algorithm from [Zhang 2021](https://arxiv.org/abs/2108.03782).

# Arguments

- post::Posterior: is the posterior you want to fit
- ndraws::Int: is the number of draws from the fit to the posterior
- init_params: Initial location to start pathfinder. If nothing the starting location
               will be a random draw from the prior.

# Returns

- q::Distributions::MvNormal: ELBO maximizing variational approximation to the posterior.
                              **Note** this is the approximation in the transformed space,
                              not the original model parameter space.
- ϕ::TupleVector: `ndraws` from th variational estimated draws from the target distribution
- logqϕ::Vector{Int}: logdensity of draws `ϕ`

# Notes
For a list of potential keyword arguments please the the [Pathfinder.jl](https://sethaxen.github.io/Pathfinder.jl/dev/)
documentation.
"""
function Pathfinder.pathfinder(post::Posterior,
                               ndraws::Int;
                               init_params=nothing,
                               kwargs...)
    # transform to unit space
    tpost = asflat(post)

    # set initial location, if nothing provided grab a draw from the prior
    p0 = init_params
    if isnothing(p0)
        p0 = rand(post.prior)
    end
    θ0 = inverse(tpost.transform, p0)
    ℓ(x) = logdensity(tpost, x)

    q, ϕ, logqϕ = pathfinder(ℓ, θ0, ndraws; kwargs...)
    out = transform_and_logjac.(Ref(tpost.transform), eachcol(ϕ))
    tϕ = first.(out)
    tlogqϕ = logqϕ .+ last.(out)
    return q, TupleVector(tϕ), tlogqϕ
end


"""
    `mulitpathfinder(post::Posterior, ndraws::Int; init_params=50, kwargs...)`

Estimates the posterior using the path finder algorithm from [Zhang 2021](https://arxiv.org/abs/2108.03782).

# Arguments

- post::Posterior: is the posterior you want to fit
- ndraws::Int: is the number of draws from the fit to the posterior
- init_params: If an integers this is the number of pathfinder algorithms you want to run.
               Otherwise, it must be a Vector where each element is a starting location of
               the model.

# Returns

- q::Distributions::MixtureModel: Uniformly weights mixture of ELBO-maximizing MvNormals.
                                  **Note** this normal approximation is done in the transformed
                                  space, not the model parameter space!
- ϕ::TupleVector: `ndraws` from th variational estimated draws from the target distribution
- component_inds::Vector{Int}: Indices `k` of components in `q` from which each element in
                               ϕ is drawn.

# Notes
For a list of potential keyword arguments please the the [Pathfinder.jl](https://sethaxen.github.io/Pathfinder.jl/dev/)
documentation.
"""
function Pathfinder.multipathfinder(post::Posterior, ndraws::Int; init_params=50, kwargs...)
    tpost = asflat(post)
    ℓ(x) = logdensity(tpost, x)

    # if init_params is an integer then we draw init_params from the prior
    if init_params isa Integer
        θ0 = inverse.(Ref(tpost.transform), rand(post.prior, init_params))
    else
        θ0 = inverse.(Ref(tpost.transform), init_params)
    end


    q, ϕ, inds = multipathfinder(ℓ, θ0, ndraws; kwargs...)
    out = transform.(Ref(tpost), eachcol(ϕ))
    return q, TupleVector(out), inds
end


"""
    `mulitpathfinder(post::TransformedPosterior, ndraws::Int; init_params=50, kwargs...)`

Estimates the posterior using the path finder algorithm from [Zhang 2021](https://arxiv.org/abs/2108.03782).

# Arguments

- post::Posterior: is the posterior you want to fit
- ndraws::Int: is the number of draws from the fit to the posterior
- init_params: If an integers this is the number of pathfinder algorithms you want to run.
               Otherwise, it must be a Vector where each element is a starting location of
               the model.

# Returns

- q::Distributions::MixtureModel: Uniformly weights mixture of ELBO-maximizing MvNormals.
                                  **Note** this normal approximation is done in the transformed
                                  space, not the model parameter space!
- ϕ::TupleVector: `ndraws` from th variational estimated draws from the target distribution
- component_inds::Vector{Int}: Indices `k` of components in `q` from which each element in
                               ϕ is drawn.

# Notes
For a list of potential keyword arguments please the the [Pathfinder.jl](https://sethaxen.github.io/Pathfinder.jl/dev/)
documentation.
"""
function Pathfinder.multipathfinder(tpost::TransformedPosterior, ndraws::Int; init_params=50, kwargs...)
    ℓ(x) = logdensity(tpost, x)

    # if init_params is an integer then we draw init_params from the prior
    if init_params isa Integer
        θ0 = inverse.(Ref(tpost.transform), rand(post.prior, init_params))
    else
        θ0 = init_params
    end

    #println(θ0)
    q, ϕ, inds = multipathfinder(ℓ, θ0, ndraws; kwargs...)
    out = collect(eachcol(ϕ))
    return q, out, inds
end
