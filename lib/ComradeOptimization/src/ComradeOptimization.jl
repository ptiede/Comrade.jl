module ComradeOptimization

using Comrade: VLBIPosterior, asflat, ascube, transform
using Reexport
using Distributions
using ForwardDiff
using LinearAlgebra
import SciMLBase

function __init__()
    @warn "ComradeOptimization is deprecated. Optimization.jl is now an extension."
end

end
