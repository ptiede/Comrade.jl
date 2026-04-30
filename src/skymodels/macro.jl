export @sky, @skymodel


_kwname(x) = x isa Symbol ? x : x.args[1]

function _scan_used_names(stmts, candidates::Vector{Symbol})
    used = Set{Symbol}()
    candset = Set(candidates)
    walk(::Any) = nothing
    function walk(node::Symbol)
        if node in candset
            push!(used, node)
        end
        return nothing
    end
    function walk(node::Expr)
        for a in node.args
            walk(a)
        end
        return nothing
    end
    for s in stmts
        walk(s)
    end
    return Symbol[n for n in candidates if n in used]
end

function _contains_tilde(node)
    if Meta.isexpr(node, :call) && length(node.args) ≥ 1 && node.args[1] === :~
        return true
    end
    if node isa Expr
        return any(_contains_tilde, node.args)
    end
    return false
end

function _sky_impl(fexpr)
    Meta.isexpr(fexpr, :function) ||
        error("@sky expects a `function name(grid; kwargs...) ... end` definition")
    sig, body = fexpr.args
    Meta.isexpr(sig, :call) || error("@sky: malformed function signature")

    name = sig.args[1]
    name isa Symbol || error("@sky: expected a function name, got $(name)")

    if length(sig.args) ≥ 2 && Meta.isexpr(sig.args[2], :parameters)
        kw_args = sig.args[2].args
        pos_args = sig.args[3:end]
    else
        kw_args = Any[]
        pos_args = sig.args[2:end]
    end
    length(pos_args) == 1 ||
        error("@sky: expected exactly one positional argument (the grid), got $(length(pos_args))")
    grid_arg = pos_args[1]
    grid_arg isa Symbol ||
        error("@sky: positional argument must be a plain symbol, got $(grid_arg)")

    kw_names = Symbol[_kwname(k) for k in kw_args]
    for k in kw_args
        kn = _kwname(k)
        kn isa Symbol ||
            error("@sky: keyword arguments must be plain names (no type annotations); got $(k)")
    end

    tilde_names = Symbol[]
    tilde_exprs = Any[]
    body_stmts = Any[]
    for stmt in body.args
        if stmt isa LineNumberNode
            push!(body_stmts, stmt)
        elseif Meta.isexpr(stmt, :call) && length(stmt.args) ≥ 3 && stmt.args[1] === :~
            lhs, rhs = stmt.args[2], stmt.args[3]
            lhs isa Symbol || error("@sky: LHS of `~` must be a symbol, got $(lhs)")
            push!(tilde_names, lhs)
            push!(tilde_exprs, rhs)
        else
            _contains_tilde(stmt) && error(
                "@sky: `~` must appear at the top level of the function body, " *
                    "not nested inside another expression"
            )
            push!(body_stmts, stmt)
        end
    end

    isempty(tilde_names) &&
        error("@sky: at least one `name ~ distribution` line is required")
    allunique(tilde_names) ||
        error("@sky: duplicate sampled-parameter names in tilde lines: $(tilde_names)")

    candidate_meta_names = Symbol[kw_names; grid_arg]
    clashes = intersect(Set(tilde_names), Set(candidate_meta_names))
    isempty(clashes) ||
        error("@sky: name clash between sampled parameters and metadata: $(collect(clashes))")
    allunique(kw_names) ||
        error("@sky: duplicate keyword-argument names: $(kw_names)")
    grid_arg in kw_names &&
        error("@sky: positional grid name `$(grid_arg)` collides with a keyword argument")

    used_meta = _scan_used_names(body_stmts, candidate_meta_names)

    sky_fn_name = Symbol("_", name, "_sky")
    prior_fn_name = Symbol("_", name, "_prior")

    θ_destr = Expr(:(=), Expr(:tuple, Expr(:parameters, tilde_names...)), :θ)
    sky_body = Expr(:block, θ_destr)
    if !isempty(used_meta)
        push!(
            sky_body.args,
            Expr(:(=), Expr(:tuple, Expr(:parameters, used_meta...)), :metadata)
        )
    end
    append!(sky_body.args, body_stmts)
    sky_fn = Expr(
        :function,
        Expr(:call, sky_fn_name, :θ, :metadata),
        sky_body
    )

    prior_required_kwargs = Symbol[kw_names; grid_arg]
    prior_sig = Expr(
        :call, prior_fn_name,
        Expr(:parameters, prior_required_kwargs...)
    )
    prior_nt = Expr(
        :tuple,
        Expr(
            :parameters,
            [Expr(:kw, n, e) for (n, e) in zip(tilde_names, tilde_exprs)]...
        )
    )
    prior_fn = Expr(:function, prior_sig, Expr(:block, Expr(:return, prior_nt)))

    ctor_sig = Expr(:call, name, Expr(:parameters, kw_args...), grid_arg)
    metadata_nt = Expr(:tuple, Expr(:parameters, kw_names..., grid_arg))
    prior_call = Expr(
        :call, prior_fn_name,
        Expr(:parameters, [Expr(:kw, n, n) for n in prior_required_kwargs]...)
    )
    skymodel_call = Expr(
        :call, :SkyModel,
        Expr(:parameters, Expr(:kw, :metadata, :metadata)),
        sky_fn_name, prior_call, grid_arg
    )
    ctor_body = Expr(
        :block,
        Expr(:(=), :metadata, metadata_nt),
        Expr(:return, skymodel_call)
    )
    ctor_fn = Expr(:function, ctor_sig, ctor_body)

    return Expr(:block, sky_fn, prior_fn, ctor_fn)
end

"""
    @sky function name(grid; kwargs...)
        param₁ ~ dist₁
        param₂ ~ dist₂
        # ... model body using params and kwargs ...
        return AbstractModel
    end

Define a [`SkyModel`](@ref) constructor in a single block. Each `name ~ dist` line
contributes an entry to the prior `NamedTuple`; everything else becomes the model
body. The macro emits three definitions in the calling module:

  - `_<name>_sky(θ, metadata)` — the model function passed to `SkyModel`.
  - `_<name>_prior(; kwargs..., grid)` — builds the prior `NamedTuple`.
  - `<name>(grid; kwargs...)` — the user-facing constructor, returning a `SkyModel`.

The first positional argument is the image grid. Keyword arguments are stored in
the model `metadata` (along with the grid) and are in scope for both the prior
RHS expressions and the model body. The metadata destructure inside the model
body only unpacks names that the body actually references.

# Hierarchical priors

The RHS of `~` is any Julia expression evaluated in the prior-builder's scope. Use
this to plug in pre-built compound priors:

```julia
@sky function imager(grid; ftot, mimg, cprior)
    c    ~ cprior                                       # e.g. HierarchicalPrior or corr_image_prior
    σimg ~ truncated(Normal(0.0, 0.5); lower = 0.0)
    fg   ~ Uniform(0.0, 1.0)
    ρs   ~ ntuple(Returns(Uniform(0.01, 10.0)), 3)      # tuple-of-IID priors
    rast = apply_fluctuations(CenteredLR(), mimg, σimg .* c.params)
    pimg = parent(rast)
    @. pimg = (ftot * (1 - fg)) * pimg
    return ContinuousImage(rast, BSplinePulse{3}())
end
```

Conditional priors where one tilde RHS references another tilde's drawn value are
not supported directly. Express such cases by wrapping the chain in a
`HierarchicalPrior` on a single tilde RHS — `HierarchicalPrior` nests, so
`HierarchicalPrior(base, HierarchicalPrior(...))` covers arbitrary-depth
hyperprior chains.
"""
macro sky(fexpr)
    return esc(_sky_impl(fexpr))
end

"""
    @skymodel

Alias for [`@sky`](@ref).
"""
macro skymodel(fexpr)
    return esc(_sky_impl(fexpr))
end
