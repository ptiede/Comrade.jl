abstract type AbstractStokesImage{T,S} <: AbstractMatrix{T} end

function fov(im::AbstractStokesImage)
    return im.fovx, im.fovy
end

mutable struct StokesImage{T,S<:AbstractMatrix} <: AbstractStokesImage{T,S}
    im::S
    fovx::T
    fovy::T
    psizex::T
    psizey::T
end
function StokesImage(im, fovx, fovy)
    ny,nx = size(im)
    psizex=fovx/(nx-1)
    psizey=fovy/(nx-1)
    return StokesImage(im,
                       convert(typeof(psizex),fovx),
                       convert(typeof(psizey),fovy),
                       psizex,
                       psizey)
end


# Define the array interface
Base.IndexStyle(::Type{StokesImage{T,S}}) where {T,S} = Base.IndexStyle(S)
Base.size(im::AbstractStokesImage) = size(im.im)
Base.getindex(im::AbstractStokesImage, i::Int) = getindex(im.im, i)
Base.setindex!(im::AbstractStokesImage, x, i::Int) = setindex!(im.im, x, i)

function Base.similar(im::AbstractStokesImage, ::Type{T}) where{T}
    sim = similar(im.im, T)
    return StokesImage(sim, im.fovx, im.fovy, im.psizex, im.psizey)
end

#function Base.similar(im::AbstractStokesImage, ::Type{T}, dims::Dims) where {T}
#    fovx = im.psizex*last(dims)
#    fovy = im.psizey*first(dims)
#    sim = similar(im.im, T, dims)
#    return StokesImage(sim, fovx, fovy, im.psizex, im.psizey)
#end

#Define the broadcast interface
struct StokesImageStyle <: Broadcast.AbstractArrayStyle{2} end
StokesImageStyle(::Val{2}) = StokesImageStyle()

Base.BroadcastStyle(::Type{<:AbstractStokesImage}) = StokesImageStyle()
function Base.similar(bc::Broadcast.Broadcasted{StokesImageStyle}, ::Type{ElType}) where ElType
    #Scan inputs for StokesImage
    #print(bc.args)
    Im = _find_sim(bc)
    #print(Im)
    #fovxs = getproperty.(Ims, Ref(:fovx))
    #fovys = getproperty.(Ims, Ref(:fovy))
    #@assert all(i->i==first(fovxs), fovxs) "StokesImage fov must be equal to add"
    #@assert all(i->i==first(fovys), fovys) "StokesImage fov must be equal to add"
    return StokesImage(similar(Array{ElType}, axes(bc)), Im.fovx, Im.fovy, Im.psizex, Im.psizey)
end

#Finds the first StokesImage and uses that as the base
#TODO: If multiple StokesImages maybe I should make sure they are consistent?
_find_sim(bc::Base.Broadcast.Broadcasted) = _find_sim(bc.args)
_find_sim(args::Tuple) = _find_sim(_find_sim(args[1]), Base.tail(args))
_find_sim(x) = x
_find_sim(::Tuple{}) = nothing
_find_sim(a::StokesImage, rest) = a
_find_sim(::Any, rest) = _find_sim(rest)

#Guards to prevent someone from adding two Images with different FOV's
function Base.:+(x::AbstractStokesImage, y::AbstractStokesImage)
    @assert fov(x) == fov(y) "StokesImages must share same field of view"
    return x .+ y
end


function flux(im::AbstractStokesImage{T,S}) where {T,S}
    sum = zero(T)
    @avx for i in eachindex(im.im)
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

@inline function imagepixels(im::StokesImage)
    px,py = pixelsizes(im)
    return range(-im.fovx/2, im.fovx/2, step=px),
           range(-im.fovy/2, im.fovy/2, step=py)
end
