# Shared AST helpers for the `@sky` (src/skymodels/macro.jl) and `@instrument`
# (src/instrument/macro.jl) macros. Included before both so the two macros can share a
# single definition instead of keeping private copies.

# Name of a keyword-argument signature entry: `x` for `x` and `x = default` for `x = ...`.
_kwname(x) = x isa Symbol ? x : x.args[1]

# Does `node` contain a top-level-or-nested `~` (tilde) call? Used to reject tildes nested
# inside other expressions in the macro bodies.
function _contains_tilde(node)
    if Meta.isexpr(node, :call) && length(node.args) ≥ 1 && node.args[1] === :~
        return true
    end
    if node isa Expr
        return any(_contains_tilde, node.args)
    end
    return false
end
