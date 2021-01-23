abstract type AbstractStokesMatrix{T} <: AbstractMatrix{T} end


mutable struct StokesMatrix{T,S<:AbstractMatrix{T}} <: AbstractStokesMatrix{T}
    im::S
    fovx::T
    fovy::T
end


Base.size(im::StokesMatrix{T,S}) where {T,S} = size(im.im)
Base.IndexStyle(::StokesMatrix{T,S}) where {T,S} = IndexLinear()
Base.getindex(im::StokesMatrix, i::Int) = getindex(im.im, i)
Base.setindex!(im::StokesMatrix, x, i::Int) = setindex!(im.im, x, i)


function flux(im::StokesMatrix{T,S}) where {T,S}
    sum = zero(T)
    @avx for i in eachindex(im)
        sum += im[i]
    end
    ny,nx = size(im)
    psizex,psizey = im.fovx/max(nx-1,1), im.fovy/max(ny-1,1)
    sum *= psizex*psizey

    return sum
end

@inline function pixelsizes(im::StokesMatrix)
    ny,nx = size(im)
    return im.fovx/max(nx-1,1), im.fovy/max(ny-1,1)
end

@inline function pixel_iterator(im::StokesMatrix)
    px,py = pixelsizes(im)
    return range(-im.fovx/2, im.fovx/2, step=px),
           range(-im.fovy/2, im.fovy/2, step=py)
end


function stokesmatrix!(im::StokesMatrix{T,S},
                      m::AbstractModel
                     ) where {T,S}
    ny,nx = size(im)
    psizex = im.fovx/max(nx-1,1)
    psizey = im.fovy/max(ny-1,1)
    flux = zero(T)
    @inbounds @simd for I in CartesianIndices(im)
        iy,ix = Tuple(I)
        x = -im.fovx/2 + psizex*(ix-1)
        y = -im.fovy/2 + psizey*(iy-1)
        tmp = intensity(m, x, y)
        im[I] = tmp
    end
    return nothing
end

function stokesmatrix(m::AbstractModel{T}, nx, ny, fovx::T, fovy::T) where {T}
    psizex = fovx/max(nx-1,1)
    psizey = fovy/max(ny-1,1)
    im = zeros(T, ny, nx)
    @inbounds @simd for I in CartesianIndices(im)
        iy,ix = Tuple(I)
        x = -fovx/2 + psizex*(ix-1)
        y = -fovy/2 + psizey*(iy-1)
        tmp = intensity(m, x, y)
        im[I] = tmp
    end
    return StokesMatrix(im, fovx, fovy)
end
