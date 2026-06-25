export @instrument


_inst_kwname(x) = x isa Symbol ? x : x.args[1]

function _inst_contains_tilde(node)
    if Meta.isexpr(node, :call) && length(node.args) ≥ 1 && node.args[1] === :~
        return true
    end
    if node isa Expr
        return any(_inst_contains_tilde, node.args)
    end
    return false
end

# The param-bearing Jones matrices: their constructor takes a single `param_map`
# function `x -> ...` mapping the per-site NamedTuple `x` to the matrix entries. The
# macro rewrites these so the user can reference the sampled-parameter names directly
# (without the dummy `x`). The parameter-free matrices (`JonesR`, `JonesF`, `JonesConst`)
# are not touched.
const _INST_PARAM_JONES = (:SingleStokesGain, :JonesG, :JonesD, :GenericJones)

_inst_jonesname(f::Symbol) = f
function _inst_jonesname(f::Expr)
    return (f.head === :. && length(f.args) == 2 && f.args[2] isa QuoteNode) ?
        f.args[2].value : nothing
end
_inst_jonesname(::Any) = nothing

# Does `node` reference any of `names` as a bare symbol? Field accesses (`x.lg`) are
# skipped so explicit `x -> ...` closures are recognised as not using bare param names.
function _inst_references(node, names)
    node isa Symbol && return node in names
    if node isa Expr
        node.head === :. && return _inst_references(node.args[1], names)
        return any(a -> _inst_references(a, names), node.args)
    end
    return false
end

# Replace every bare reference to a sampled-parameter name with `x.<name>`.
_inst_rewrite_params(node, names) = node
_inst_rewrite_params(node::Symbol, names) = node in names ? Expr(:., :x, QuoteNode(node)) : node
function _inst_rewrite_params(node::Expr, names)
    node.head === :. && return Expr(:., _inst_rewrite_params(node.args[1], names), node.args[2])
    return Expr(node.head, Any[_inst_rewrite_params(a, names) for a in node.args]...)
end

# Number of arguments of an anonymous-function signature (`x`, `(x, y)`, or `()`).
function _inst_func_argcount(lhs)
    lhs isa Symbol && return 1
    (lhs isa Expr && lhs.head === :tuple) && return length(lhs.args)
    return 1
end

_inst_argsig(lhs::Symbol) = (lhs,)
_inst_argsig(lhs::Expr) = lhs.head === :tuple ? Tuple(lhs.args) : (lhs,)
_inst_argsig(::Any) = ()

_inst_asblock(body) = (body isa Expr && body.head === :block) ? body : Expr(:block, body)

# Collect the names bound by a top-level assignment LHS (`a`, `a, b`, `a::T`).
_inst_collect_lhs!(set, lhs::Symbol) = (push!(set, lhs); set)
function _inst_collect_lhs!(set, lhs::Expr)
    if lhs.head === :tuple || lhs.head === :(::)
        for a in lhs.args
            a isa Symbol && push!(set, a)
        end
    end
    return set
end
_inst_collect_lhs!(set, ::Any) = set

# Context threaded through the body walk.
mutable struct _InstPmapCtx
    names::Set{Symbol}   # sampled-parameter names (rewritten to `x.<name>`)
    guard::Set{Symbol}   # kwargs + body-locals; referencing one keeps a closure
    base::Symbol         # instrument name, used to name generated functions
    genfns::Vector{Any}  # generated top-level function definitions
    counter::Int
end

# Closures get gensym'd types that do not serialize reliably, so we lift anonymous
# functions to *named* top-level functions (mirroring how `@sky` emits a named
# `_<name>_sky`). A function that captures a construction-time kwarg or a body-local is
# left as a closure so the capture still works.

# `param_map` body written without the dummy `x`: rewrite params to `x.<name>`, then lift
# to a named `_<name>_pmap_<k>(x)` (or keep a closure when it captures guarded names).
function _inst_make_pmap(bodyexpr, ctx::_InstPmapCtx)
    rewritten = _inst_rewrite_params(bodyexpr, ctx.names)
    _inst_references(bodyexpr, ctx.guard) && return Expr(:->, :x, rewritten)
    ctx.counter += 1
    fname = Symbol("_", ctx.base, "_pmap_", ctx.counter)
    push!(ctx.genfns, Expr(:function, Expr(:call, fname, :x), _inst_asblock(rewritten)))
    return fname
end

# A non-`param_map` anonymous function (e.g. the `JonesSandwich` combination function):
# lift verbatim, keeping its own argument list, to a named `_<name>_fn_<k>`.
function _inst_lift_plain(func, ctx::_InstPmapCtx)
    (func isa Expr && func.head === :->) || return func
    _inst_references(func.args[2], ctx.guard) && return func
    ctx.counter += 1
    fname = Symbol("_", ctx.base, "_fn_", ctx.counter)
    push!(
        ctx.genfns,
        Expr(:function, Expr(:call, fname, _inst_argsig(func.args[1])...), _inst_asblock(func.args[2]))
    )
    return fname
end

# `param_map` written as a call argument (bare expression/tuple, or anonymous function).
# Explicit ≥1-arg closures and param-free arguments (e.g. `JonesG(fgain)`) are untouched.
function _inst_wrap_pmap(arg, ctx::_InstPmapCtx)
    if arg isa Expr && arg.head === :->
        if _inst_func_argcount(arg.args[1]) == 0
            body = arg.args[2]
            return _inst_references(body, ctx.names) ? _inst_make_pmap(body, ctx) : arg
        end
        return arg
    end
    return _inst_references(arg, ctx.names) ? _inst_make_pmap(arg, ctx) : arg
end

# `param_map` written as a `do`-block (`JonesG() do ... end`).
function _inst_wrap_func(func, ctx::_InstPmapCtx)
    (func isa Expr && func.head === :->) || return func
    if _inst_func_argcount(func.args[1]) == 0
        body = func.args[2]
        return _inst_references(body, ctx.names) ? _inst_make_pmap(body, ctx) : func
    end
    return func
end

# Walk the Jones-building body and lift `param_map`s and combination functions.
function _inst_transform_pmaps(node, ctx::_InstPmapCtx)
    node isa Expr || return node
    if node.head === :do
        call = node.args[1]
        if call isa Expr && call.head === :call
            nm = _inst_jonesname(call.args[1])
            if nm in _INST_PARAM_JONES
                newcall = Expr(:call, Any[_inst_transform_pmaps(a, ctx) for a in call.args]...)
                wrapped = _inst_wrap_func(node.args[2], ctx)
                return wrapped isa Symbol ? Expr(:call, newcall.args..., wrapped) :
                    Expr(:do, newcall, wrapped)
            elseif nm === :JonesSandwich
                # The `do`-function is the combination function (passed first). Lift it
                # and turn `JonesSandwich(mats...) do ... end` into `JonesSandwich(fn, mats...)`.
                newcall = Expr(:call, Any[_inst_transform_pmaps(a, ctx) for a in call.args]...)
                lifted = _inst_lift_plain(node.args[2], ctx)
                return lifted isa Symbol ?
                    Expr(:call, newcall.args[1], lifted, newcall.args[2:end]...) :
                    Expr(:do, newcall, lifted)
            end
        end
    elseif node.head === :call
        nm = _inst_jonesname(node.args[1])
        if nm in _INST_PARAM_JONES
            newargs = Any[node.args[1]]
            for a in node.args[2:end]
                push!(newargs, (a isa Expr && a.head === :parameters) ? a : _inst_wrap_pmap(a, ctx))
            end
            return Expr(:call, newargs...)
        elseif nm === :JonesSandwich
            newargs = Any[node.args[1]]
            for (i, a) in enumerate(node.args[2:end])
                if i == 1 && a isa Expr && a.head === :->
                    push!(newargs, _inst_lift_plain(a, ctx))   # combination function
                else
                    push!(newargs, _inst_transform_pmaps(a, ctx))
                end
            end
            return Expr(:call, newargs...)
        end
    end
    return Expr(node.head, Any[_inst_transform_pmaps(a, ctx) for a in node.args]...)
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
        _inst_kwname(k) isa Symbol ||
            error("@instrument: keyword arguments must be plain names (no type annotations); got $(k)")
    end

    # Separate the special `refbasis` kwarg from the user kwargs. The user may
    # supply their own default for `refbasis`; otherwise inject `CirBasis()`.
    refbasis_kw = nothing
    user_kw_args = Any[]
    for k in kw_args
        if _inst_kwname(k) === :refbasis
            refbasis_kw = k
        else
            push!(user_kw_args, k)
        end
    end
    if refbasis_kw === nothing
        refbasis_kw = Expr(:kw, :refbasis, :(CirBasis()))
    end
    user_kw_names = Symbol[_inst_kwname(k) for k in user_kw_args]

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
            _inst_contains_tilde(stmt) && error(
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

    # Lift `param_map`s and combination functions to named top-level functions so the body
    # references the sampled-parameter names directly (no dummy `x`) and avoids closures.
    locals = Set{Symbol}()
    for s in body_stmts
        s isa Expr && s.head === :(=) && _inst_collect_lhs!(locals, s.args[1])
    end
    guard = union(Set(user_kw_names), locals)
    ctx = _InstPmapCtx(Set(tilde_names), guard, name, Any[], 0)
    body_stmts = Any[_inst_transform_pmaps(s, ctx) for s in body_stmts]

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
Each `name ~ ArrayPrior(...)` line contributes an entry to the instrument prior
`NamedTuple`; everything else becomes the body that builds the Jones matrix. The macro
emits the following definitions in the calling module:

  - `_<name>_pmap_<k>(x)` / `_<name>_fn_<k>(...)` — named functions for the `param_map`s
    and combination functions found in the body (see below).
  - `_<name>_jones(; kwargs...)` — builds and returns the composed `AbstractJonesMatrix`.
  - `_<name>_prior(; kwargs...)` — builds the prior `NamedTuple` of [`ArrayPrior`](@ref)s.
  - `<name>(; refbasis = CirBasis(), kwargs...)` — the user-facing constructor, returning
    an `InstrumentModel`.

The instrument model takes no positional arguments. Keyword arguments are tunable
constants in scope for both the prior RHS expressions and the Jones body. `refbasis` is
treated specially: it is forwarded only to `InstrumentModel` and defaults to `CirBasis()`
if not listed. Any number of Jones terms is supported: just build and compose them (e.g.
with [`JonesSandwich`](@ref)) and `return` the result. Unlike `@sky`, `~` lines are
optional, so a parameter-free response model is allowed.

# Writing `param_map`s

For the param-bearing Jones matrices ([`SingleStokesGain`](@ref), [`JonesG`](@ref),
[`JonesD`](@ref), [`GenericJones`](@ref)) you can reference the sampled-parameter names
directly. Both a bare expression and a zero-argument `do` block work. The macro lifts
each `param_map` (and the `JonesSandwich` combination function) to a *named* top-level
function — anonymous closures have gensym'd types that do not serialize reliably:

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
    D = JonesD((complex(dRx, dRy), complex(dLx, dLy)))   # bare expression
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
    return SingleStokesGain(exp(complex(lg, gp)))
end
```

A function that captures a construction-time `kwarg` or a body-local is kept as a closure
(so the capture still works) rather than lifted. The explicit `x -> ...` form is also
accepted and is required for custom `AbstractJonesMatrix` types the macro does not know
about, e.g. `SingleStokesGain(x -> exp(complex(x.lg, x.gp)))`.
"""
macro instrument(fexpr)
    return esc(_instrument_impl(fexpr))
end
