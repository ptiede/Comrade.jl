using Pkg; Pkg.activate(@__DIR__)
# Really only need this the first time
Pkg.instantiate()
using NFFT
using LinearAlgebra
import Comrade
using Zygote


"""
    ObservedNUFT

Container type for a non-uniform Fourier transform (NUFT).
This stores the uv-positions that the model will be sampled at in the Fourier domain,
allowing certain transformtion factors (e.g., NUFT matrix) to be cached.

This is an internal type, an end user should instead create this using [`NFFTAlg`](@ref NFFTAlg)
or [`DFTAlg`](@ref DFTAlg).
"""
struct ObservedNUFT{A<:Comrade.NUFT, T} <: Comrade.NUFT
    """
    Which NUFT algorithm to use (e.g. NFFTAlg or DFTAlg)
    """
    alg::A
    """
    uv positions of the NUFT transform. This is used for precomputation.
    """
    uv::T
end


struct DFTAlg <: Comrade.NUFT end

"""
    DFTAlg(u::AbstractArray, v::AbstractArray)

Create an algorithm object using the direct Fourier transform object using the uv positions
`u`, `v` allowing for a more efficient transform.
"""
function DFTAlg(u::AbstractArray, v::AbstractArray)
    uv = Matrix{eltype(u)}(undef, 2, length(u))
    uv[1,:] .= u
    uv[2,:] .= v
    return ObservedNUFT(DFTAlg(), uv)
end



# Creates the DFT matrix, this needs to be moved to the gpu (after construction probably)
function plan_nuft(alg::ObservedNUFT{<:DFTAlg}, img)
    uv = alg.uv
    xitr, yitr = Comrade.imagepixels(img)
    dft = similar(img, Complex{eltype(img)}, (size(uv, 2), size(img)...))
    @fastmath for i in axes(img,2), j in axes(img,1), k in axes(uv,2)
        u = uv[1,k]
        v = uv[2,k]
        dft[k, j, i] = cispi(-2(u*xitr[i] + v*yitr[j]))
    end
    # reshape to a matrix so we can take advantage of an easy BLAS call
    return reshape(dft, size(uv,2), :)
end


"""
    create_cache(alg::ObservedNUFT, plan , phases, img)

Create a cache for the DFT algorithm with precomputed `plan`, `phases` and `img`.
This is an internal version.
"""
function create_cache(alg::ObservedNUFT{<:DFTAlg}, plan, phases, img)
    return NUFTCache(alg, plan, phases, img.pulse, reshape(img.img, :))
end


function dft(dft_matrix::AbstractMatrix, img::AbstractMatrix)
    dft_matrix*reshape(img, :)
end

# For gpu yoy first want to make all arrays a gpu array.
nvis = 500
u = randn(nvis)
v = randn(nvis)

# Make some synthetic data
# Also should move to gpu array
sigma = fill(0.01, nvis)
vis = complex.(exp.(-2π^2. *(u.^2 .+ v.^2))) .+ sigma.*randn(nvis)

struct TestPost{A, V, F, S}
    dft_mat::A
    vis::V
    fov::F
    sigma::S
end

function TestPost(u, v, vis, sigma, fov, npix)
    alg = DFTAlg(u, v)
    img = Comrade.IntensityMap(similar(u, (npix, npix)), fov, fov)
    dft_mat = plan_nuft(alg, img)
    return TestPost(dft_mat, vis, fov, sigma)
end


function (m::TestPost)(I)
    vmodel = dft(m.dft_mat, I)
    # Now let's do a likelihood
    return -sum(abs2, (vmodel .- m.vis)./m.sigma)/2
end


post = TestPost(u, v, vis, sigma, 10.0, 32)

post(randn(32, 32))

f(x) = post(x)

using Zygote

using BenchmarkTools

@benchmark $(post)($rand(32,32))
@benchmark Zygote.gradient($post, $(rand(32,32)))
