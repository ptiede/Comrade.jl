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

function make_pullback(ℓ, autodiff::Function)
    function ∇ℓ(x)
        return (ℓ(x), autodiff(x))
    end
end
