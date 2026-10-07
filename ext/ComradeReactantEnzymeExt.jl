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

# The median wall time of `f()` once the device runs at its working speed. A device left idle
# (for example during a host-side refit) runs its first calls far slower than it does while
# sampling, so the calls before timing span `warmup` seconds rather than a fixed count.
function _steady_seconds(f; warmup = 0.5, mincalls = 10, ncalls = 21)
    t0 = time()
    n = 0
    while n < mincalls || time() - t0 < warmup
        f()
        n += 1
    end
    ts = [@elapsed(f()) for _ in 1:ncalls]
    return sort!(ts)[cld(ncalls, 2)]
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
    t = _steady_seconds(() -> g(tpost, x))
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
