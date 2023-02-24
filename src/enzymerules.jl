using Enzyme
using Enzyme: EnzymeRules
using ArgCheck

function EnzymeRules.inactive(::typeof(ArgCheck.check), args...)
    return nothing
end
