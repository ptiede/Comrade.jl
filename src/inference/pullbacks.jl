"""
    make_pullback(ℓ, autodiff::AD.AbstractBackend)

Create the pullback function using the autodiff backend `autodiff`.

# Note

This is an internal function and not part of the public API.
"""
function make_pullback(ℓ, autodiff::AD.AbstractBackend)
    function ∇ℓ(x)
        res = AD.value_and_gradient(autodiff, ℓ, x)
        return (first(res), first(last(res)))
    end
end

# This is custom since Zygote and AD have some weird performanc regression currently
function make_pullback(ℓ, ::AD.ReverseRuleConfigBackend)
    function ∇ℓ(x)
        f, b = AD.Zygote.pullback(ℓ, x)
        return (f, first(b(1.0)))
    end
end

"""
    make_pullback(ℓ, grad::Function)

Create the pullback function using the function `grad` which should return the gradient
of ℓ.

# Note

This is an internal function and shouldn't be typically used by an end-user.
"""
function make_pullback(ℓ, grad::Function)
    function ∇ℓ(x)
        return (ℓ(x), grad(x))
    end
end
