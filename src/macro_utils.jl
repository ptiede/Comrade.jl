# Shared AST helpers for the `@sky` (src/skymodels/macro.jl), `@instrument`, and `@jones`
# (src/instrument/macro.jl) macros. Included before all of them so they can share a single
# definition instead of keeping private copies.

# Name of a keyword-argument signature entry: `x` for `x` and `x = default` for `x = ...`.
_kwname(x) = x isa Symbol ? x : x.args[1]

# Does `node` contain a nested *prior* `~` declaration (`lhs ~ rhs`)? Used to reject tildes
# nested inside other expressions in the macro bodies. Only the binary form (`Expr(:call, :~,
# lhs, rhs)`, i.e. ≥ 3 args) is a prior — this matches the detection in `_split_tildes`. The
# unary form `~x` (bitwise NOT, 2 args) is an ordinary operator and must not be flagged, so a
# body statement like `mask = ~flags` is left alone.
function _contains_tilde(node)
    if Meta.isexpr(node, :call) && length(node.args) ≥ 3 && node.args[1] === :~
        return true
    end
    if node isa Expr
        return any(_contains_tilde, node.args)
    end
    return false
end

# Of `candidates`, which names appear as a bare symbol anywhere in `stmts`? Returns them in
# `candidates` order. Used to build the minimal destructure of a parameter / metadata
# NamedTuple inside a generated function body.
function _scan_used_names(stmts, candidates::Vector{Symbol})
    used = Set{Symbol}()
    candset = Set(candidates)
    walk(::Any) = nothing
    function walk(node::Symbol)
        node in candset && push!(used, node)
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

# Split the statements of a macro body into `~` (tilde) declarations and the rest. Returns
# `(tilde_names, tilde_exprs, body_stmts)`: a `name ~ rhs` line contributes `name` and `rhs`;
# `LineNumberNode`s and every other statement are kept (in order) in `body_stmts`. A `~`
# nested inside another expression is an error. `macroname` is used in error messages.
function _split_tildes(stmts, macroname::AbstractString)
    tilde_names = Symbol[]
    tilde_exprs = Any[]
    body_stmts = Any[]
    for stmt in stmts
        if stmt isa LineNumberNode
            push!(body_stmts, stmt)
        elseif Meta.isexpr(stmt, :call) && length(stmt.args) ≥ 3 && stmt.args[1] === :~
            lhs, rhs = stmt.args[2], stmt.args[3]
            lhs isa Symbol || error("$macroname: LHS of `~` must be a symbol, got $(lhs)")
            push!(tilde_names, lhs)
            push!(tilde_exprs, rhs)
        else
            _contains_tilde(stmt) && error(
                "$macroname: `~` must appear at the top level of the body, " *
                    "not nested inside another expression"
            )
            push!(body_stmts, stmt)
        end
    end
    return tilde_names, tilde_exprs, body_stmts
end
