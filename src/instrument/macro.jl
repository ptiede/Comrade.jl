export @instrument, @jones

# `_kwname`, `_contains_tilde`, `_scan_used_names`, and `_split_tildes` are shared with the
# `@sky` macro; see src/macro_utils.jl.

# Name of a macro reference: `Symbol("@jones")` for a bare `@jones`, or the tail symbol for a
# qualified `Comrade.@jones`.
_macroname(s::Symbol) = s
_macroname(e::Expr) = (e.head === :. && length(e.args) == 2 && e.args[2] isa QuoteNode) ?
    e.args[2].value : nothing
_macroname(::Any) = nothing

_is_jones_call(node) =
    Meta.isexpr(node, :macrocall) && _macroname(node.args[1]) === Symbol("@jones")

# Argument names of an anonymous-function signature (`x`, `(x, y)`, `()`, or a do-block header).
_argsig(lhs::Symbol) = (lhs,)
_argsig(lhs::Expr) = lhs.head === :tuple ? Tuple(lhs.args) : (lhs,)
_argsig(::Any) = ()

_asblock(body) = (body isa Expr && body.head === :block) ? body : Expr(:block, body)

# Emit a fresh named function `_<base><suffix><k>(argsig...) = body` into `defs` and return its
# name. The generated functions are spliced *inside* `_<name>_jones`, so they may capture the
# instrument's construction-time keyword arguments while still being named (a named inner
# function has a clean `nameof` and serializes — unlike an anonymous closure).
function _emit_fn!(defs, base, suffix, counter, argsig, body)
    counter[] += 1
    fname = Symbol("_", base, suffix, counter[])
    push!(defs, Expr(:function, Expr(:call, fname, argsig...), _asblock(body)))
    return fname
end

# Expand one `@jones begin ... end` block. Returns a NamedTuple:
#   `jones`  — the expression that builds the `AbstractJonesMatrix` (referencing the lifted
#              param_map by name),
#   `fndefs` — the generated param_map function definition(s),
#   `priors` — the `name => prior_rhs` pairs contributed by this block's `~` lines.
#
# The block's `~` lines are its priors; the rest is its param_map. A parameterized block must
# end in a constructor call `Ctor(posarg; kws...)` with exactly one positional argument: the
# preceding statements plus that argument become `Ctor`'s param_map, a named function of the
# per-site parameter NamedTuple `x` with a leading `(; used...) = x` destructure. No constructor
# name is consulted, so custom `AbstractJonesMatrix` types work automatically. A param-free
# block (no `~` lines) is emitted verbatim — its value is the Jones matrix.
function _jones_impl(args, base::Symbol, counter::Base.RefValue{Int})
    length(args) == 1 ||
        error("@jones: expected `@jones begin ... end`, got $(length(args)) arguments")
    block = args[1]
    Meta.isexpr(block, :block) ||
        error("@jones: expected a `begin ... end` block body, got `$(block)`")

    names, exprs, stmts = _split_tildes(block.args, "@jones")
    allunique(names) ||
        error("@jones: duplicate sampled-parameter names: $(names)")

    realidx = findlast(s -> !(s isa LineNumberNode), stmts)
    realidx === nothing && error("@jones: empty block body")
    finalstmt = stmts[realidx]

    # Param-free block: run the body verbatim (stripping a trailing `return`, which would
    # otherwise return from the enclosing `_<name>_jones`). Its value is the Jones matrix.
    if isempty(names)
        value = Meta.isexpr(finalstmt, :return) ? finalstmt.args[1] : finalstmt
        blk = Expr(:block, stmts[1:(realidx - 1)]..., value)
        return (jones = blk, fndefs = Any[], priors = Pair{Symbol, Any}[])
    end

    # Parameterized block: the final statement builds the Jones from a single positional arg.
    call = Meta.isexpr(finalstmt, :return) ? finalstmt.args[1] : finalstmt
    Meta.isexpr(call, :call) || error(
        "@jones: the final statement of a parameterized block (one with `~` priors) must be a " *
            "Jones constructor call with a single positional argument, e.g. `return JonesG((g1, g2))`; " *
            "got `$(finalstmt)`. Control-flow finals (`if`/`let`/`begin`) are not decomposable into a " *
            "param_map — compute the entries into a single value first, then `return Ctor(entries)`."
    )
    ctor = call.args[1]
    params_expr = nothing
    positionals = Any[]
    for a in call.args[2:end]
        if a isa Expr && a.head === :parameters
            params_expr = a
        else
            push!(positionals, a)
        end
    end
    length(positionals) == 1 || error(
        "@jones: the final Jones constructor must take exactly one positional argument " *
            "(its param_map), got $(length(positionals)) in `$(call)`"
    )
    posarg = positionals[1]

    preceding = stmts[1:(realidx - 1)]
    used = _scan_used_names(Any[preceding..., posarg], names)
    pmbody = Any[]
    isempty(used) ||
        push!(pmbody, Expr(:(=), Expr(:tuple, Expr(:parameters, used...)), :x))
    append!(pmbody, preceding)
    push!(pmbody, posarg)

    fndefs = Any[]
    pmname = _emit_fn!(fndefs, base, "_pmap_", counter, (:x,), Expr(:block, pmbody...))

    newargs = Any[ctor]
    params_expr === nothing || push!(newargs, params_expr)
    push!(newargs, pmname)
    return (
        jones = Expr(:call, newargs...),
        fndefs = fndefs,
        priors = Pair{Symbol, Any}[n => e for (n, e) in zip(names, exprs)],
    )
end

# Replace every `@jones` macrocall in `node` with the Jones-building expression it expands to,
# accumulating the generated param_map definitions and prior pairs.
function _subst_jones(node, base, counter, fndefs, priors)
    if _is_jones_call(node)
        uargs = Any[a for a in node.args[2:end] if !(a isa LineNumberNode)]
        res = _jones_impl(uargs, base, counter)
        append!(fndefs, res.fndefs)
        append!(priors, res.priors)
        return res.jones
    end
    node isa Expr || return node
    return Expr(node.head, Any[_subst_jones(a, base, counter, fndefs, priors) for a in node.args]...)
end

# Lift anonymous functions (the `JonesSandwich` combination function, written as a `do`-block
# or `->`) in the composition to named inner functions, so the model's `jones_map` is a named
# (serializable) function rather than a gensym closure. These functions take the assembled
# Jones matrices as arguments and never reference sampled parameters, so no destructure or
# capture handling is needed.
#
# A lifted function's body is kept *verbatim* — we deliberately do not recurse into it. A nested
# anonymous function inside the body may close over the outer lambda's arguments (e.g.
# `do g, d; map(x -> g * x, d); end`); hoisting that inner lambda to a sibling named function
# would put `g` out of scope and crash at evaluation time. Leaving nested closures inline keeps
# their lexical scope intact, and only the outermost (now named) function matters for
# serialization. Sibling anonymous functions elsewhere in the expression are still lifted by the
# general recursion below.
function _lift_anons(node, base, counter, fndefs)
    node isa Expr || return node
    if node.head === :->
        return _emit_fn!(fndefs, base, "_fn_", counter, _argsig(node.args[1]), node.args[2])
    elseif node.head === :do
        call = node.args[1]
        func = node.args[2]
        fname = _emit_fn!(fndefs, base, "_fn_", counter, _argsig(func.args[1]), func.args[2])
        rest = Any[_lift_anons(a, base, counter, fndefs) for a in call.args[2:end]]
        return Expr(:call, call.args[1], fname, rest...)
    end
    return Expr(node.head, Any[_lift_anons(a, base, counter, fndefs) for a in node.args]...)
end

function _instrument_impl(fexpr)
    Meta.isexpr(fexpr, :function) ||
        error("@instrument expects a `function name(; kwargs...) ... end` definition")
    sig, body = fexpr.args
    Meta.isexpr(sig, :call) || error("@instrument: malformed function signature")

    name = sig.args[1]
    name isa Symbol || error("@instrument: expected a function name, got $(name)")

    if length(sig.args) ≥ 2 && Meta.isexpr(sig.args[2], :parameters)
        kw_args = sig.args[2].args
        pos_args = sig.args[3:end]
    else
        kw_args = Any[]
        pos_args = sig.args[2:end]
    end
    isempty(pos_args) ||
        error("@instrument: expected no positional arguments (the instrument model takes none), got $(length(pos_args))")

    for k in kw_args
        _kwname(k) isa Symbol ||
            error("@instrument: keyword arguments must be plain names (no type annotations); got $(k)")
    end

    # Separate the special `refbasis` kwarg from the user kwargs. The user may supply their own
    # default for `refbasis`; otherwise inject `CirBasis()`.
    refbasis_kw = nothing
    user_kw_args = Any[]
    for k in kw_args
        if _kwname(k) === :refbasis
            refbasis_kw = k
        else
            push!(user_kw_args, k)
        end
    end
    refbasis_kw === nothing && (refbasis_kw = Expr(:kw, :refbasis, :(CirBasis())))
    user_kw_names = Symbol[_kwname(k) for k in user_kw_args]

    # Walk the body: expand `@jones` blocks (collecting their param_maps + priors) and lift the
    # composition's combination function(s) to named inner functions. `~` is only allowed
    # inside a `@jones` block.
    counter = Ref(0)
    fndefs = Any[]
    prior_pairs = Pair{Symbol, Any}[]
    body_stmts = Any[]
    for stmt in body.args
        if stmt isa LineNumberNode
            push!(body_stmts, stmt)
        elseif Meta.isexpr(stmt, :call) && length(stmt.args) ≥ 3 && stmt.args[1] === :~
            error(
                "@instrument: `~` priors must be declared inside a `@jones` block, " *
                    "not at the top level of the instrument body"
            )
        else
            s = _subst_jones(stmt, name, counter, fndefs, prior_pairs)
            s = _lift_anons(s, name, counter, fndefs)
            push!(body_stmts, s)
        end
    end

    tilde_names = Symbol[p.first for p in prior_pairs]
    tilde_exprs = Any[p.second for p in prior_pairs]
    allunique(tilde_names) ||
        error("@instrument: duplicate sampled-parameter name across @jones blocks: $(tilde_names)")
    allunique(user_kw_names) ||
        error("@instrument: duplicate keyword-argument names: $(user_kw_names)")
    clashes = intersect(Set(tilde_names), Set(user_kw_names))
    isempty(clashes) ||
        error("@instrument: name clash between sampled parameters and keyword arguments: $(collect(clashes))")
    :refbasis in tilde_names &&
        error("@instrument: `refbasis` is reserved and cannot be used as a sampled-parameter name")

    jones_fn_name = Symbol("_", name, "_jones")
    prior_fn_name = Symbol("_", name, "_prior")

    # `_<name>_jones(; user_kwargs...)` — defines the lifted param_map / combination functions,
    # then runs the (jones-substituted) body, returning the composed `AbstractJonesMatrix`.
    jones_fn = Expr(
        :function,
        Expr(:call, jones_fn_name, Expr(:parameters, user_kw_args...)),
        Expr(:block, fndefs..., body_stmts...)
    )

    # `_<name>_prior(; user_kwargs...)` — returns the (flat) prior NamedTuple. Empty when there
    # are no `~` lines (a parameter-free instrument).
    prior_nt = Expr(
        :tuple,
        Expr(:parameters, [Expr(:kw, n, e) for (n, e) in zip(tilde_names, tilde_exprs)]...)
    )
    prior_fn = Expr(
        :function,
        Expr(:call, prior_fn_name, Expr(:parameters, user_kw_args...)),
        Expr(:block, Expr(:return, prior_nt))
    )

    # `<name>(; refbasis = CirBasis(), user_kwargs...)` — the user-facing constructor.
    ctor_sig = Expr(:call, name, Expr(:parameters, refbasis_kw, user_kw_args...))
    fwd_kwargs = Expr(:parameters, [Expr(:kw, n, n) for n in user_kw_names]...)
    prior_call = Expr(:call, prior_fn_name, fwd_kwargs)
    jones_call = Expr(:call, jones_fn_name, fwd_kwargs)
    instrument_call = Expr(
        :call, :InstrumentModel,
        Expr(:parameters, Expr(:kw, :refbasis, :refbasis)),
        :jones, :prior
    )
    ctor_body = Expr(
        :block,
        Expr(:(=), :prior, prior_call),
        Expr(:(=), :jones, jones_call),
        Expr(:return, instrument_call)
    )
    ctor_fn = Expr(:function, ctor_sig, ctor_body)

    return Expr(:block, jones_fn, prior_fn, ctor_fn)
end

"""
    @jones begin
        param₁ ~ ArrayPrior(...)
        # ... statements building the matrix entries from the params ...
        return JonesG(entries)
    end

Build one Jones term inside an [`@instrument`](@ref) block. Each `name ~ ArrayPrior(...)` line
is a prior for this term; the rest of the block is how the we paraemeterized the Jones matrix.

A block with no `~` lines is parameter-free and is emitted verbatim (e.g. `JonesR`):

```julia
R = @jones begin
    return JonesR(; add_fr = true)
end
```

`@jones` is only meaningful inside [`@instrument`](@ref) (its `~` priors have nowhere else to
go); using it standalone with `~` lines is an error.
"""
macro jones(args...)
    uargs = Any[a for a in args if !(a isa LineNumberNode)]
    res = _jones_impl(uargs, :jones, Ref(0))
    isempty(res.priors) || error(
        "@jones: `~` priors are only valid inside an `@instrument` block. A standalone " *
            "`@jones` may only build a parameter-free Jones term."
    )
    return esc(Expr(:block, res.fndefs..., res.jones))
end

"""
    @instrument function name(; refbasis = CirBasis(), kwargs...)
        G = @jones begin ... end
        D = @jones begin ... end
        # ... build and compose the Jones terms ...
        return AbstractJonesMatrix
    end

Define an [`InstrumentModel`](@ref) constructor in a single block, mirroring [`@sky`](@ref).
Each Jones term is built with a [`@jones`](@ref) block, which carries its own `~` priors and
`param_map`; the priors of all blocks are collected into the instrument prior. The body then
composes the built Jones terms (e.g. with [`JonesSandwich`](@ref)) and returns the result.

The instrument model takes no positional arguments. Keyword arguments are tunable constants in
scope for the `~` prior expressions and the `@jones` bodies. `refbasis` is treated specially:
it is forwarded only to `InstrumentModel` and defaults to `CirBasis()` if not listed. Priors
are declared *only* inside `@jones` blocks — a top-level `~` is an error. `~` lines are
optional overall, so a parameter-free response model is allowed.

The `JonesSandwich` combination function is lifted to a named function automatically (named
functions, unlike closures, serialize reliably).

Example:
```julia
@instrument function fullpol(; refbasis = CirBasis(), gain_std = 0.1, dterm_std = 0.2)
    G = @jones begin
        lgR   ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, gain_std)))
        gpR   ~ ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2)));
                           refant = SEFDReference(0.0), phase = true)
        lgrat ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, gain_std)))
        gprat ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
        gR = exp(complex(lgR, gpR))
        gL = gR * exp(complex(lgrat, gprat))
        return JonesG((gR, gL))
    end
    D = @jones begin
        dRx ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, dterm_std)))
        dRy ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, dterm_std)))
        dLx ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, dterm_std)))
        dLy ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, dterm_std)))
        return JonesD((complex(dRx, dRy), complex(dLx, dLy)))
    end
    R = @jones begin
        return JonesR(; add_fr = true)
    end
    return JonesSandwich(G, D, R) do g, d, r
        adjoint(r) * g * d * r
    end
end

intm = fullpol(; gain_std = 0.05)
```

A Stokes-I example is a single block:

```julia
@instrument function gainsonly(; refbasis = CirBasis())
    return @jones begin
        lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
        gp ~ ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))))
        return SingleStokesGain(exp(complex(lg, gp)))
    end
end
```
"""
macro instrument(fexpr)
    return esc(_instrument_impl(fexpr))
end
