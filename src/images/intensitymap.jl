export IntensityMap, fov, imagepixels, pixelsizes


mutable struct IntensityMap{T,S<:AbstractMatrix, K<:Pulse} <: AbstractIntensityMap{T,S}
    im::S
    fovx::T
    fovy::T
    psizex::T
    psizey::T
    pulse::K
end

fov(m::AbstractIntensityMap) = (m.fovx, m.fovy)


function IntensityMap(im, fovx, fovy, pulse=DeltaPulse())
    ny,nx = size(im)
    psizex=fovx/(nx-1)
    psizey=fovy/(ny-1)
    return IntensityMap(im,
                       convert(typeof(psizex),fovx),
                       convert(typeof(psizey),fovy),
                       psizex,
                       psizey,
                       pulse)
end


"""
    $(SIGNATURES)
Computes the flux of a intensity map
"""
function flux(im::AbstractIntensityMap{T,S}) where {T,S}
    sum = zero(T)
    x,y = imagepixels(im)
    @inbounds for i in axes(im,1), j in axes(im, 2)
        xx = x[i]
        yy = y[j]
        sum += im[j,i]*intensity_point(im.pulse, xx, yy)
    end
    return sum*prod(pixelsizes(im))
end



# Define the array interface
Base.IndexStyle(::Type{<: IntensityMap{T,S,K}}) where {T,S,K} = Base.IndexStyle(S)
Base.size(im::AbstractIntensityMap) = size(im.im)
Base.getindex(im::AbstractIntensityMap, i::Int) = getindex(im.im, i)
Base.getindex(im::AbstractIntensityMap, I...) = getindex(im.im, I...)
Base.setindex!(im::AbstractIntensityMap, x, i::Int) = setindex!(im.im, x, i)

function Base.similar(im::IntensityMap, ::Type{T}) where{T}
    sim = similar(im.im, T)
    return IntensityMap(sim, im.fovx, im.fovy, im.psizex, im.psizey, im.pulse)
end

#function Base.similar(im::AbstractIntensityMap, ::Type{T}, dims::Dims) where {T}
#    fovx = im.psizex*last(dims)
#    fovy = im.psizey*first(dims)
#    sim = similar(im.im, T, dims)
#    return IntensityMap(sim, fovx, fovy, im.psizex, im.psizey)
#end

#Define the broadcast interface
struct IntensityMapStyle <: Broadcast.AbstractArrayStyle{2} end
IntensityMapStyle(::Val{2}) = IntensityMapStyle()

Base.BroadcastStyle(::Type{<:AbstractIntensityMap}) = IntensityMapStyle()
function Base.similar(bc::Broadcast.Broadcasted{IntensityMapStyle}, ::Type{ElType}) where ElType
    #Scan inputs for IntensityMap
    #print(bc.args)
    Im = _find_sim(bc)
    #fovxs = getproperty.(Im, Ref(:fovx))
    #fovys = getproperty.(Im, Ref(:fovy))
    #@assert all(i->i==first(fovxs), fovxs) "IntensityMap fov must be equal to add"
    #@assert all(i->i==first(fovys), fovys) "IntensityMap fov must be equal to add"
    return IntensityMap(similar(Array{ElType}, axes(bc)), Im.fovx, Im.fovy, Im.psizex, Im.psizey, Im.pulse)
end

#Finds the first IntensityMap and uses that as the base
#TODO: If multiply IntensityMaps maybe I should make sure they are consistent?
_find_sim(bc::Base.Broadcast.Broadcasted) = _find_sim(bc.args)
_find_sim(args::Tuple) = _find_sim(_find_sim(args[1]), Base.tail(args))
_find_sim(x) = x
_find_sim(::Tuple{}) = nothing
_find_sim(a::AbstractIntensityMap, rest) = a
_find_sim(::Any, rest) = _find_sim(rest)

#Guards to prevent someone from adding two Images with different FOV's
function Base.:+(x::AbstractIntensityMap, y::AbstractIntensityMap)
    @assert fov(x) == fov(y) "IntensityMaps must share same field of view"
    return x .+ y
end



@inline function pixelsizes(im::AbstractIntensityMap)
    ny,nx = size(im)
    return im.fovx/max(nx-1,1), im.fovy/max(ny-1,1)
end

@inline function imagepixels(im::AbstractIntensityMap)
    px,py = pixelsizes(im)
    return range(-im.fovx/2, im.fovx/2, step=px),
           range(-im.fovy/2, im.fovy/2, step=py)
end
