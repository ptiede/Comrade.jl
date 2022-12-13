padfac(::DFTAlg) = 1

# We don't pad a DFT since it is already correct
padimage(::DFTAlg, img::IntensityMap) = img

"""
    DFTAlg(obs::EHTObservation)

Create an algorithm object using the direct Fourier transform object from the observation
`obs`. This will extract the uv positions from the observation to allow for a more efficient
FT cache.
"""
function DFTAlg(obs::EHTObservation)
    (;U,V) = getuv(obs.config)
    return DFTAlg(U, V)
end

"""
    DFTAlg(ac::ArrayConfiguration)

Create an algorithm object using the direct Fourier transform object from the array configuration
`ac`. This will extract the uv positions from the observation to allow for a more efficient
FT cache.
"""
function DFTAlg(ac::ArrayConfiguration)
    (;U,V) = getuv(ac)
    return DFTAlg(U, V)
end

"""
    DFTAlg(u::AbstractArray, v::AbstractArray)

Create an algorithm object using the direct Fourier transform object using the uv positions
`u`, `v` allowing for a more efficient transform.
"""
function DFTAlg(u::AbstractArray, v::AbstractArray)
    @argcheck length(u) == length(v)
    uv = Matrix{eltype(u)}(undef, 2, length(u))
    uv[1,:] .= u
    uv[2,:] .= v
    return ObservedNUFT(DFTAlg(), uv)
end


# internal function that creates an DFT matrix/plan to use used for the img.
function plan_nuft(alg::ObservedNUFT{<:DFTAlg}, img::Union{IntensityMap{T,2}, StokesIntensityMap{T,2}}) where {T}
    uv = alg.uv
    xitr, yitr = imagepixels(img)
    dft = similar(parent(img), Complex{eltype(uv)}, size(uv,2), size(img)...)
    @fastmath for i in axes(img,2), j in axes(img,1), k in axes(uv,2)
        u = uv[1,k]
        v = uv[2,k]
        # - sign is taken care of in _visibilities
        dft[k, j, i] = cispi(-2(u*xitr[j] + v*yitr[i]))
    end
    # reshape to a matrix so we can take advantage of an easy BLAS call
    return reshape(dft, size(uv,2), :)
end

# internal function to make the phases to phase center the image.
function make_phases(alg::ObservedNUFT{<:DFTAlg}, img::IntensityMapTypes, pulse::Pulse)
    u = @view alg.uv[1,:]
    v = @view alg.uv[2,:]
    # We don't need to correct for the phase offset here since that
    # is taken care of in plan_nuft for DFTAlg
    dx, dy = pixelsizes(img)
    return _visibilities(stretched(pulse, dx, dy), u, v, 0.0, 0.0)
end

"""
    create_cache(alg::ObservedNUFT, plan , phases, img)

Create a cache for the DFT algorithm with precomputed `plan`, `phases` and `img`.
This is an internal version.
"""
function create_cache(alg::ObservedNUFT{<:DFTAlg}, plan, phases, img, pulse::Pulse)
    return NUFTCache(alg, plan, phases, pulse, reshape(img, :))
end

function _visibilities(m::ModelImage{M,I,<:NUFTCache{A}},
    u, v, time, freq) where {M,I<:StokesIntensityMap,A<:ObservedNUFT{<:DFTAlg}}
checkuv(m.cache.alg.uv, u, v)
    visI =  nuft(m.cache.plan, complex.(reshape((ComradeBase.AxisKeys.keyless_unname(m.cache.img).I), :)))
    visQ =  nuft(m.cache.plan, complex.(reshape((ComradeBase.AxisKeys.keyless_unname(m.cache.img).Q), :)))
    visU =  nuft(m.cache.plan, complex.(reshape((ComradeBase.AxisKeys.keyless_unname(m.cache.img).U), :)))
    visV =  nuft(m.cache.plan, complex.(reshape((ComradeBase.AxisKeys.keyless_unname(m.cache.img).V), :)))
    r = StructArray{StokesParams{eltype(visI)}}((I=visI, Q=visQ, U=visU, V=visV)).*m.cache.phases
    return r
end
