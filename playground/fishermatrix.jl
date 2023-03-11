export fishermatrix


"""
    $(SIGNATURES)

Compute the Fisher information matrix of a `model` transformed by `t` using the array configuration
`ac`. You can also optionally pass the `ad_type` as the last argument which is which
backend to `AbstractDifferentiation` you want to use to compute derivatives.

`model` is assumed to be a function that takes the `NamedTuple`, `θ` and returns
a `<:Comrade.AbstractModel` that coincides with the visibility model you want to consider.

This returns a `Tuple` where the first entry is Fisher information metric and the second
is distribution assuming that the mean is the transformed `θ` and the covariance matrix
is the inverse of the Fisher information.
"""
function fishermatrix(model, t, θ::NamedTuple, ac::ArrayConfiguration)
    v = Base.Fix2(visibilities, ac)
    tr = Base.Fix1(transform, t)

    # split into real and imaginary since ForwardDiff struggles with complex functions
    fr = real∘v∘model∘tr
    fi = imag∘v∘model∘tr
    x0 = inverse(t, θ)
    vr = first(AD.jacobian(ad_type, fr, x0))
    vi = first(AD.jacobian(ad_type, fi, x0))
    v1 = complex.(vr, vi)
    M = Symmetric(real.(v1'*v1) .+ eps())
    h = M*x0
    return M, Dists.MvNormalCanon(h, M)
end
