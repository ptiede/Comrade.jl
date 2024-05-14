module ComradeNested

using Comrade
using AbstractMCMC
using Reexport
using Random

export resample_equal

function __init__()
    @warn "ComradeNested is deprecated. Dynesty.jl is now an extension."
end

end
