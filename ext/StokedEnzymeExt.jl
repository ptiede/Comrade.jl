module StokedEnzymeExt

using Enzyme
using Stoked
using LogDensityProblems

LogDensityProblems.dimension(d::Stoked.TransformedVLBIPosterior) = dimension(d)
LogDensityProblems.capabilities(::Type{<:Stoked.TransformedVLBIPosterior}) = LogDensityProblems.LogDensityOrder{1}()


function LogDensityProblems.logdensity_and_gradient(d::Stoked.TransformedVLBIPosterior, x::AbstractArray)
    mode = Enzyme.EnzymeCore.WithPrimal(Stoked.admode(d))
    dx = zero(x)
    (_, y) = autodiff(mode, Stoked.logdensityof, Active, Const(d), Duplicated(x, dx))
    return y, dx
end



end