export @instrument

# `_kwname` and `_contains_tilde` are shared with the `@sky` macro; see src/macro_utils.jl.

# --- Symbol scan ----------------------------------------------------------------------
# Sampled-parameter names are *reserved* within the instrument body: they may be referenced
# only as bare value names inside a *param_map* — a zero-argument function (`JonesG() do ...
# end` or `JonesG(() -> ...)`) passed to a Jones constructor. The macro turns such a function
# into `_<name>_pmap_<k>(x)` by injecting the per-site parameter NamedTuple `x` and
# destructuring the parameters it uses (`(; p...) = x`). The function's *arity* — not its
# constructor's name — is what marks it as a param_map, so any `AbstractJonesMatrix` taking a
# `param_map` (built-in or user-defined) works the same way. Everywhere else a bare parameter
# name is a *leak* and is rejected at macro-expansion time.
#
# `_inst_refs!` collects the bare-symbol references to `names` in `node`. The only positions
# skipped are ones that are not value references: the field of a `.` access (`x.lg`) and the
# key of a keyword argument (`f(lg = 1)`). A broadcast `f.(args)` is NOT a field access (its
# second arg is a tuple of values), so its arguments are descended into.
function _inst_refs!(found, node, names)
    if node isa Symbol
        node in names && push!(found, node)
        return found
    end
    node isa Expr || return found
    h = node.head
    if h === :.
        _inst_refs!(found, node.args[1], names)
        if length(node.args) ≥ 2 && !(node.args[2] isa QuoteNode)
            _inst_refs!(found, node.args[2], names)   # broadcast call args
        end
    elseif h === :kw && length(node.args) == 2
        _inst_refs!(found, node.args[2], names)
    else
        for a in node.args
            _inst_refs!(found, a, names)
        end
    end
    return found
end

# Does `node` reference any of `names`?
_inst_references(node, names) = !isempty(_inst_refs!(Set{Symbol}(), node, names))

# Number of arguments of an anonymous-function signature (`x`, `(x, y)`, `()`, or a do-block
# header). An empty tuple (`() -> ...` / `do; ...; end`) has zero arguments.
function _inst_func_argcount(lhs)
    lhs isa Symbol && return 1
    (lhs isa Expr && lhs.head === :tuple) && return length(lhs.args)
    return 1
end

_inst_argsig(lhs::Symbol) = (lhs,)
_inst_argsig(lhs::Expr) = lhs.head === :tuple ? Tuple(lhs.args) : (lhs,)
_inst_argsig(::Any) = ()

_inst_asblock(body) = (body isa Expr && body.head === :block) ? body : Expr(:block, body)

# Collect the names bound by an assignment LHS: `a`, destructuring `a, b` and `(; a, b)`
# (possibly nested), `a::T` (the variable, not the type), and short-form function defs
# `f(args) = ...` (the function name).
_inst_collect_lhs!(set, lhs::Symbol) = (push!(set, lhs); set)
function _inst_collect_lhs!(set, lhs::Expr)
    if lhs.head === :tuple || lhs.head === :parameters
        for a in lhs.args
            _inst_collect_lhs!(set, a)
        end
    elseif lhs.head === :(::)
        _inst_collect_lhs!(set, lhs.args[1])
    elseif lhs.head === :call && !isempty(lhs.args)
        lhs.args[1] isa Symbol && push!(set, lhs.args[1])
    end
    return set
end
_inst_collect_lhs!(set, ::Any) = set

# Collect every name a top-level body statement binds in the enclosing scope: assignments
# (`=`), `local`/`global` declarations, and long-form function definitions. These names —
# the construction-time locals — together with the user kwargs make up the `guard`: a
# param_map that references one of them must stay a closure so the capture works. Over-
# collecting here is safe: it only keeps more param_maps as closures instead of lifting them.
function _inst_collect_bound!(set, s)
    s isa Expr || return set
    if s.head === :(=)
        _inst_collect_lhs!(set, s.args[1])
    elseif s.head === :local || s.head === :global
        for a in s.args
            (a isa Expr && a.head === :(=)) ? _inst_collect_lhs!(set, a.args[1]) :
                _inst_collect_lhs!(set, a)
        end
    elseif s.head === :function && !isempty(s.args)
        _inst_collect_lhs!(set, s.args[1])
    end
    return set
end

# Names bound by a `let` / `for` / generator binder section (`a = ...`, destructuring, or a
# bare `a`), descending through binding blocks. Used to widen the `guard` as the body walk
# descends into a `let` / `for` / comprehension, so a param_map nested inside one treats the
# binder as a captured local (and is kept a closure).
function _inst_loop_binders!(set, spec)
    if spec isa Symbol
        push!(set, spec)
    elseif spec isa Expr
        if spec.head === :block || spec.head === :filter
            for a in spec.args
                _inst_loop_binders!(set, a)
            end
        elseif spec.head === :(=)
            _inst_collect_lhs!(set, spec.args[1])
        end
    end
    return set
end

# Context threaded through the body walk.
mutable struct _InstPmapCtx
    names::Set{Symbol}    # sampled-parameter names (reserved; destructured from `x`)
    base::Symbol          # instrument name, used to name generated functions
    genfns::Vector{Any}   # generated top-level function definitions
    counter::Int
    captured::Set{Symbol} # guard names that forced a closure (for the serialization warning)
    leaked::Set{Symbol}   # param names reached outside a param_map (the leak check)
end

# Emit a fresh named top-level function `_<base><suffix><k>(argsig...) = body` and return its
# name. Shared by the param_map and combination-function lifters.
function _inst_emit_fn!(ctx::_InstPmapCtx, suffix, argsig, body)
    ctx.counter += 1
    fname = Symbol("_", ctx.base, suffix, ctx.counter)
    push!(ctx.genfns, Expr(:function, Expr(:call, fname, argsig...), _inst_asblock(body)))
    return fname
end

# Build a param_map body that binds the sampled parameters it uses as locals via a leading
# `(; p₁, p₂) = x` destructure, then runs the user body verbatim. The body is ordinary Julia
# in a scope where the parameters are plain locals, so Julia itself resolves every
# shadowing/binding form (nested closures, `let`, `for`, comprehensions, broadcast, `where`,
# ...) with no AST surgery on our part.
function _inst_destructure_body(bodyexpr, names)
    used = sort!(collect(_inst_refs!(Set{Symbol}(), bodyexpr, names)))
    blk = _inst_asblock(bodyexpr)
    isempty(used) && return blk
    destr = Expr(:(=), Expr(:tuple, Expr(:parameters, used...)), :x)
    return Expr(:block, destr, blk.args...)
end

# Process an anonymous function (`->` or the function of a `do`-block) appearing in the
# instrument body, returning either a Symbol (the name of a lifted top-level function) or an
# `Expr(:->, ...)` (a closure kept inline because it captures a construction-time local).
#
# The function's arity decides its role — no constructor name is consulted:
#   * zero arguments   -> a *param_map*. Inject the per-site parameter NamedTuple `x`,
#     destructure the sampled parameters it uses, and lift it to `_<name>_pmap_<k>(x)`.
#   * one or more args -> a *plain function* (a `JonesSandwich` combination function, or an
#     explicit `x -> f(x.lg)` closure for a Jones type that bare syntax cannot express). Its
#     arguments are supplied by the caller, so it is lifted verbatim to `_<name>_fn_<k>`.
#     Sampled parameters are *not* in scope here, so a bare param reference is a leak.
# A function that captures a construction-time kwarg or body-local cannot become a top-level
# function (it could not see the captured value), so it is kept as an inline closure and
# recorded for the serialization warning.
function _inst_process_func(sig, body, ctx::_InstPmapCtx, guard)
    captures = _inst_references(body, guard)
    if _inst_func_argcount(sig) == 0
        destructured = _inst_destructure_body(body, ctx.names)
        if captures
            _inst_refs!(ctx.captured, body, guard)
            return Expr(:->, Expr(:tuple, :x), destructured)
        end
        return _inst_emit_fn!(ctx, "_pmap_", (:x,), destructured)
    end
    _inst_refs!(ctx.leaked, body, ctx.names)   # bare params are reserved; not in scope here
    if captures
        _inst_refs!(ctx.captured, body, guard)
        return Expr(:->, sig, body)
    end
    return _inst_emit_fn!(ctx, "_fn_", _inst_argsig(sig), body)
end

# Walk the Jones-building body, lifting every anonymous function (param_maps and combination
# functions alike) to a named top-level function. A bare sampled-parameter name reached by
# this walk — i.e. anywhere outside a zero-argument param_map — is recorded as a leak. `guard`
# is the set of construction-time local names in effect at `node` (kwargs + body-locals,
# widened by any `let` / `for` / generator binder we descend through).
function _inst_xform(node, ctx::_InstPmapCtx, guard)
    if node isa Symbol
        node in ctx.names && push!(ctx.leaked, node)
        return node
    end
    node isa Expr || return node
    if node.head === :->
        return _inst_process_func(node.args[1], node.args[2], ctx, guard)
    elseif node.head === :do
        # A `do`-block passes its function as the *first* positional argument, so
        # `f(args...) do ... end` is `f(do_fn, args...)`. Lift the function and rebuild as a
        # plain call (which also handles the case where it is kept as an inline closure).
        call = node.args[1]
        func = node.args[2]
        pf = _inst_process_func(func.args[1], func.args[2], ctx, guard)
        rest = Any[_inst_xform(a, ctx, guard) for a in call.args[2:end]]
        return Expr(:call, call.args[1], pf, rest...)
    end
    # Descend, widening `guard` with any binders introduced by a nested `let` / `for` /
    # generator so a param_map inside such a block sees them as captured (and is kept a
    # closure).
    child = guard
    if node.head === :let || node.head === :for
        child = _inst_loop_binders!(copy(guard), node.args[1])
    elseif node.head === :generator
        child = copy(guard)
        for a in node.args[2:end]
            _inst_loop_binders!(child, a)
        end
    end
    return Expr(node.head, Any[_inst_xform(a, ctx, child) for a in node.args]...)
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

    # Separate the special `refbasis` kwarg from the user kwargs. The user may
    # supply their own default for `refbasis`; otherwise inject `CirBasis()`.
    refbasis_kw = nothing
    user_kw_args = Any[]
    for k in kw_args
        if _kwname(k) === :refbasis
            refbasis_kw = k
        else
            push!(user_kw_args, k)
        end
    end
    if refbasis_kw === nothing
        refbasis_kw = Expr(:kw, :refbasis, :(CirBasis()))
    end
    user_kw_names = Symbol[_kwname(k) for k in user_kw_args]

    tilde_names = Symbol[]
    tilde_exprs = Any[]
    body_stmts = Any[]
    for stmt in body.args
        if stmt isa LineNumberNode
            push!(body_stmts, stmt)
        elseif Meta.isexpr(stmt, :call) && length(stmt.args) ≥ 3 && stmt.args[1] === :~
            lhs, rhs = stmt.args[2], stmt.args[3]
            lhs isa Symbol || error("@instrument: LHS of `~` must be a symbol, got $(lhs)")
            push!(tilde_names, lhs)
            push!(tilde_exprs, rhs)
        else
            _contains_tilde(stmt) && error(
                "@instrument: `~` must appear at the top level of the function body, " *
                    "not nested inside another expression"
            )
            push!(body_stmts, stmt)
        end
    end

    allunique(tilde_names) ||
        error("@instrument: duplicate sampled-parameter names in tilde lines: $(tilde_names)")
    allunique(user_kw_names) ||
        error("@instrument: duplicate keyword-argument names: $(user_kw_names)")
    clashes = intersect(Set(tilde_names), Set(user_kw_names))
    isempty(clashes) ||
        error("@instrument: name clash between sampled parameters and keyword arguments: $(collect(clashes))")
    :refbasis in tilde_names &&
        error("@instrument: `refbasis` is reserved and cannot be used as a sampled-parameter name")

    # Lift every anonymous function (param_maps and combination functions) to named top-level
    # functions so the body references the sampled-parameter names directly (no dummy `x`) and
    # avoids closures. The `guard` is the set of construction-time locals (kwargs + body-
    # locals): a param_map that references one of them is kept as a closure so the capture
    # works.
    locals = Set{Symbol}()
    for s in body_stmts
        _inst_collect_bound!(locals, s)
    end
    guard = union(Set(user_kw_names), locals)
    ctx = _InstPmapCtx(Set(tilde_names), name, Any[], 0, Set{Symbol}(), Set{Symbol}())
    body_stmts = Any[_inst_xform(s, ctx, guard) for s in body_stmts]

    # Sampled parameters are reserved within the instrument body: they may be referenced only
    # inside a zero-argument param_map. `_inst_xform` records every param name it reaches
    # outside such a param_map (the Jones composition, a `JonesSandwich` combination function,
    # or a ≥1-arg explicit closure); we also reject params used in a `~` prior RHS, where they
    # would leak — unbound — into the generated prior function. A leaked name would fail with
    # a confusing `UndefVarError` (or silently read a same-named global), so we reject it here,
    # at macro-expansion time.
    for e in tilde_exprs
        _inst_refs!(ctx.leaked, e, ctx.names)
    end
    if !isempty(ctx.leaked)
        leaked = join(sort!(collect(ctx.leaked)), ", ")
        error(
            "@instrument: sampled parameter(s) $(leaked) are used outside a param_map. " *
                "Parameter names are reserved in the instrument body: reference them only " *
                "inside a zero-argument param_map, e.g. `JonesG() do; ...; end`. Do not use " *
                "them in the Jones composition, a `JonesSandwich` combination function, a " *
                "`~` prior RHS, or as a local / loop / argument name."
        )
    end

    # A param_map that captures a construction-time kwarg or a body-local is kept as a
    # closure (so the capture works), but anonymous closures have gensym'd types that do not
    # serialize reliably — the very thing this macro tries to avoid. Warn so the user can
    # restructure if they need to serialize the model.
    if !isempty(ctx.captured)
        @warn "@instrument: in `$(name)`, a param_map captures $(join(sort!(collect(ctx.captured)), ", ")) " *
            "and was kept as a closure instead of a named function; such closures may not " *
            "serialize reliably. Inline the captured value(s) as literals in the param_map " *
            "if you need to serialize the instrument model."
    end

    jones_fn_name = Symbol("_", name, "_jones")
    prior_fn_name = Symbol("_", name, "_prior")

    # `_<name>_jones(; user_kwargs...)` — runs the (tilde-stripped) body verbatim,
    # which builds and returns the composed `AbstractJonesMatrix`.
    jones_fn = Expr(
        :function,
        Expr(:call, jones_fn_name, Expr(:parameters, user_kw_args...)),
        Expr(:block, body_stmts...)
    )

    # `_<name>_prior(; user_kwargs...)` — returns the prior NamedTuple. Each RHS is
    # evaluated with the user kwargs in scope. Empty when there are no `~` lines.
    prior_nt = Expr(
        :tuple,
        Expr(
            :parameters,
            [Expr(:kw, n, e) for (n, e) in zip(tilde_names, tilde_exprs)]...
        )
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

    return Expr(:block, ctx.genfns..., jones_fn, prior_fn, ctor_fn)
end

"""
    @instrument function name(; refbasis = CirBasis(), kwargs...)
        param₁ ~ ArrayPrior(...)
        param₂ ~ ArrayPrior(...)
        # ... body building Jones terms and composing them ...
        return AbstractJonesMatrix
    end

Define an [`InstrumentModel`](@ref) constructor in a single block, mirroring [`@sky`](@ref).
Each `name ~ ArrayPrior(...)` line contributes an entry to the instrument prior.

The instrument model takes no positional arguments. Keyword arguments are tunable constants
in scope for both the prior RHS expressions and the Jones body. `refbasis` is treated
specially: it is forwarded only to `InstrumentModel` and defaults to `CirBasis()` if not
listed. Any number of Jones terms is supported: build and compose them (e.g. with
[`JonesSandwich`](@ref)) and `return` the result. Unlike `@sky`, `~` lines are optional, so a
parameter-free response model is allowed.

# param_maps

A parameterized Jones term (`SingleStokesGain`, `JonesG`, `JonesD`, `GenericJones`, or any
`AbstractJonesMatrix` that takes a `param_map`) is given its parameter map as a
**zero-argument function** (a `do`-block) that references the sampled-parameter names
directly:

```julia
G = JonesG() do
    gR = exp(complex(lgR, gpR))
    gL = gR * exp(complex(lgrat, gprat))
    return gR, gL
end
D = JonesD() do
    return (complex(dRx, dRy), complex(dLx, dLy))
end
```

The macro lifts each zero-argument function into a named top-level function of the per-site
parameter NamedTuple, destructuring the parameters it uses. Named functions, unlike closures,
serialize reliably. The arity is what marks a function as a param_map, so the macro never
inspects constructor names — a user-defined `AbstractJonesMatrix` works exactly like the
built-ins. The [`JonesSandwich`](@ref) combination function takes the assembled matrices, so
it is written *with* arguments and is lifted verbatim:

```julia
return JonesSandwich(G, D, R) do g, d, r
    return adjoint(r) * g * d * r
end
```

!!! note
    Sampled-parameter names are **reserved** within the instrument body: they may be used
    only as bare value references inside a zero-argument param_map. Using a parameter name
    anywhere else — in the Jones composition, a `JonesSandwich` combination function, a `~`
    prior RHS, or as a local / loop / argument name — is a macro-expansion error.

Example:
```julia
@instrument function fullpol(; refbasis = CirBasis(), gain_std = 0.1, dterm_std = 0.2)
    lgR   ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, gain_std)))
    gpR   ~ ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2)));
                       refant = SEFDReference(0.0), phase = true)
    lgrat ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, gain_std)))
    gprat ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
    dRx   ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, dterm_std)))
    dRy   ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, dterm_std)))
    dLx   ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, dterm_std)))
    dLy   ~ ArrayPrior(IIDSitePrior(TrackSeg(), VLBIGaussian(0.0, dterm_std)))

    G = JonesG() do                       # zero-arg do-block, params used directly
        gR = exp(complex(lgR, gpR))
        gL = gR * exp(complex(lgrat, gprat))
        return gR, gL
    end
    D = JonesD() do                       # zero-arg do-block
        return complex(dRx, dRy), complex(dLx, dLy)
    end
    R = JonesR(; add_fr = true)
    return JonesSandwich(G, D, R) do g, d, r
        adjoint(r) * g * d * r
    end
end

intm = fullpol(; gain_std = 0.05)
```

A Stokes-I example:

```julia
@instrument function gainsonly(; refbasis = CirBasis())
    lg ~ ArrayPrior(IIDSitePrior(ScanSeg(), VLBIGaussian(0.0, 0.1)))
    gp ~ ArrayPrior(IIDSitePrior(ScanSeg(), DiagonalVonMises(0.0, inv(π^2))))
    return SingleStokesGain() do
        exp(complex(lg, gp))
    end
end
```
"""
macro instrument(fexpr)
    return esc(_instrument_impl(fexpr))
end
