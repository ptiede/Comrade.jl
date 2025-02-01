module ComradeEnzymeExt

using Enzyme
using Comrade
using LogDensityProblems

LogDensityProblems.dimension(d::Comrade.TransformedVLBIPosterior) = dimension(d)
LogDensityProblems.capabilities(::Type{<:Comrade.TransformedVLBIPosterior}) = LogDensityProblems.LogDensityOrder{1}()


function LogDensityProblems.logdensity_and_gradient(d::Comrade.TransformedVLBIPosterior, x::AbstractArray)
    mode = Enzyme.EnzymeCore.WithPrimal(Comrade.admode(d))
    dx = zero(x)
    (_, y) = autodiff(mode, Comrade.logdensityof, Active, Const(d), Duplicated(x, dx))
    return y, dx
end



end