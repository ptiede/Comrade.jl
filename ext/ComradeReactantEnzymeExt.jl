module ComradeReactantEnzymeExt

using Comrade
using ComradeBase: ReactantEx
using Enzyme
using Reactant

# Gradient of the flat log-density on the device. `Enzyme.gradient` returns one
# derivative per differentiated argument and `tpost` is `Const`, so the derivative for
# `x` is the last entry. `set_strong_zero` sends 0*Inf and 0*NaN in the reverse pass to
# 0, without which a stiff image model has no finite gradient at all.
function _device_flat_grad(tpost, x)
    derivs, val = Enzyme.gradient(
        Enzyme.set_strong_zero(Enzyme.ReverseWithPrimal),
        Comrade.logdensityof, Enzyme.Const(tpost), x
    )
    return last(derivs), val
end

# Gradient wall times by the device array the cost depends on (the preconditioner's `V`, or
# `tpost` itself without one), held weakly. An in-place refit keeps the shapes, so the time
# stays valid.
const _GRADIENT_SECONDS = Dict{UInt, Tuple{WeakRef, Float64}}()

function Comrade._gradient_seconds(tpost, x::Reactant.AbstractConcreteArray)
    pre = Comrade._transport_pre(tpost)
    obj = isnothing(pre) ? tpost : pre.V
    e = get(_GRADIENT_SECONDS, objectid(obj), nothing)
    (!isnothing(e) && e[1].value === obj) && return e[2]
    g = Reactant.@compile sync = true _device_flat_grad(tpost, x)
    g(tpost, x)
    t = sort([@elapsed(g(tpost, x)) for _ in 1:5])[3]
    filter!(kv -> !isnothing(kv[2][1].value), _GRADIENT_SECONDS)
    _GRADIENT_SECONDS[objectid(obj)] = (WeakRef(obj), t)
    return t
end

function Comrade._compiled_score(post::Comrade.VLBIPosterior, x0::AbstractVector, space)
    dpost = Comrade.prepare_device(post, ReactantEx())
    tflat = Comrade.maybe_transport(dpost, space)
    vg = Reactant.@compile sync = true _device_flat_grad(
        tflat, Reactant.to_rarray(collect(x0))
    )
    return function (x)
        g, _ = vg(tflat, Reactant.to_rarray(collect(x)))
        return Array(g)
    end
end

end
