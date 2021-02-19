abstract type AbstractStokesImage{T,S} end


mutable struct StokesImage{T,S<:AbstractMatrix} <: AbstractStokesImage{T,S}
    im::S
    fovx::T
    fovy::T
    psizex::T
    psizey::T
end

Base.size(im::StokesImage) = size(im.im)
Base.getindex(im::StokesImage, i::Int) = getindex(im.im, i)
Base.setindex!(im::StokesImage, x, i::Int) = setindex!(im.im, x, i)


function flux(im::StokesImage{T,S}) where {T,S}
    sum = zero(T)
    @avx for i in eachindex(im)
        sum += im[i]
    end
    ny,nx = size(im)
    psizex,psizey = im.fovx/max(nx-1,1), im.fovy/max(ny-1,1)
    sum *= psizex*psizey
    return sum
end

@inline function pixelsizes(im::StokesImage)
    ny,nx = size(im)
    return im.fovx/max(nx-1,1), im.fovy/max(ny-1,1)
end

@inline function pixel_iterator(im::StokesImage)
    px,py = pixelsizes(im)
    return range(-im.fovx/2, im.fovx/2, step=px),
           range(-im.fovy/2, im.fovy/2, step=py)
end


function stokesimage!(im::StokesImage{T,S},
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

function stokesimage(m::AbstractModel{T}, nx, ny, fovx, fovy) where {T}
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
    return StokesImage(im, fovx, fovy, psizex, psizey)
end
