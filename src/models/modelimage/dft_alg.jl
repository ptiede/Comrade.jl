padfac(::DFTAlg) = 1

# We don't pad a DFT since it is already correct
padimage(::DFTAlg, img) = img

function DFTAlg(u::AbstractArray, v::AbstractArray)
    @argcheck length(u) == length(v)
    uv = Matrix{eltype(u)}(undef, 2, length(u))
    uv[1,:] .= u
    uv[2,:] .= v
    return ObservedNUFT(DFTAlg(), uv)
end

function plan_nuft(alg::ObservedNUFT{<:DFTAlg}, img, dx, dy)
    uv = alg.uv
    xitr, yitr = imagepixels(img)
    dft = Array{Complex{eltype(uv)},3}(undef, size(uv,2), size(img)...)
    @fastmath for i in axes(img,2), j in axes(img,1), k in axes(uv,2)
        u = uv[1,k]
        v = uv[2,k]
        dft[k, j, i] = cispi(2(u*xitr[i] + v*yitr[j]))
    end
    # reshape to a matrix so we can take advantage of an easy BLAS call
    return reshape(dft, size(uv,2), :)
end

function make_phases(alg::ObservedNUFT{<:DFTAlg}, img)
    u = @view alg.uv[1,:]
    v = @view alg.uv[2,:]
    dx, dy = pixelsizes(img)
    visibility.(Ref(img.pulse), u*dx, v*dy).*dx.*dy
end

@inline function create_cache(alg::ObservedNUFT{<:DFTAlg}, plan, phases, img)
    return NUFTCache(alg, plan, phases, reshape(img.im, :))
end
