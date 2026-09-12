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

function Comrade._compiled_flat_score(post::Comrade.VLBIPosterior, x0::AbstractVector)
    dpost = Comrade.prepare_device(post, ReactantEx())
    tflat = Comrade.asflat(dpost)
    vg = Reactant.@compile sync = true _device_flat_grad(
        tflat, Reactant.to_rarray(collect(x0))
    )
    return function (x)
        g, _ = vg(tflat, Reactant.to_rarray(collect(x)))
        return Array(g)
    end
end

end
