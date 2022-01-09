if lowercase(get(ENV, "CI", "false")) == "true"
    include("install_pycall.jl")
    using PyCall
    PyCall.Conda.pip_interop(true)
    PyCall.Conda.pip("install", "ehtim")
end

using ROSEx
using Statistics
using Test

include("test_util.jl")

@testset "ROSEx.jl" begin

    include("observation.jl")

end
