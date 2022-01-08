export RImage





@doc raw"""
    $(TYPEDEF)
An image model given by a set of coefficients and a kernel response or basis function.
This corresponds to a continous image defined by a finite set of points. The defined
intensity is given by
```math
    I(x,y) = \sum_{ij} c_{ij}κ(x-x_i)κ(y-y_i).
```
An important thing to note is that the ``c_{ij}`` do not represent pixel intensities, i.e.
the κ doesn't have to be an interpolating kernel.

## Example
```julia
samples = rand(10,10)
model = RImage(samples, BSplineKernel{3})
```

## Notes

This is defined in terms of pixel response, so the image size is 1μas. To resize the image
use the scale function like with other models.

## Fields

$(FIELDS)

"""
struct RImage{S,B<:Pulse,M<:AbstractMatrix{S}} <: AbstractModel
    """ Image coefficients cᵢⱼ in expansion """
    coeff::M
    """ Image kernel/basis κ that defined the delta image response """
    kernel::B
    # pixel size in 1/pixels
    psizex::S
    # pixel size in 1/pixels
    psizey::S
    function RImage(coeff::M, basis::B) where {S,M<:AbstractMatrix{S},B}
        ny, nx = size(coeff)
        psizex = one(S)/max(nx-1, 1)
        psizey = one(S)/max(ny-1, 1)
        new{S,B,M}(coeff, basis, psizex, psizey)
    end
end
@inline visanalytic(::Type{<:RImage}) = IsAnalytic()
@inline isprimitive(::Type{<:RImage}) = IsPrimitive()

#=
struct FourierCache{C} <: ObservationCache
    cache::C
end

function FourierCache(rimage::I, obs::Observation) where {S,I<:AbstractRadioImage{S}}
    cache = zeros(Complex{S}, size(rimage)..., nsamples(obs))
    ny,nx = size(rimage)
    u = getdata(obs, :u)
    v = getdata(obs, :v)
    dx = 1/max(nx-1,1)
    dy = 1/max(ny-1,1)
    startx = -0.5
    starty = -0.5
    x = range(startx, length=nx, step=dx)
    y = range(starty, length=ny, step=dy)
    for i in eachindex(u,v)
        cache[:,:,i] .= exp.(2im*π*(u[i].*x' .+ v[i].*y))
    end
    FourierCache(cache)
end
=#


@inline function flux(model::RImage{S,B,M}) where {S,B,M}
    sum = zero(S)
    @turbo for i in eachindex(model.coeff)
        sum += model.coeff[i]
    end
    # Divide by pixel number to convert properly
    ny,nx = size(model.coeff)
    dx = 1/max(nx-1,1)
    dy = 1/max(ny-1,1)
    return sum*κflux(model.kernel)*dx*dy
end




"""
    $(SIGNATURES)
return the size of the coefficient matrix for `model`.
"""
@inline Base.size(model::RImage) = size(model.coeff)

@inline function intensity_point(model::RImage{S,M,B}, x, y, args...) where {S,M,B}
    sum = zero(S)
    ny,nx = size(model)
    dx = 1/(max(nx-1,1))
    dy = 1/(max(ny-1,1))
    #The kernel is written in terms of pixel number so we convert x to it
    @inbounds @simd for I in CartesianIndices(model.coeff)
        iy,ix = Tuple(I)
        xx = x - (-0.5 + dx*(ix-1))
        yy = y - (-0.5 + dy*(iy-1))
        sum += model.coeff[I]* κ(model.kernel, xx/dx)*κ(model.kernel, yy/dy)
    end
    # Note this will be intensity per uas
    return sum
end


@inline function visibility_point(model::RImage{S,M,B}, u, v, args...) where {S,M,B}
    sum = zero(Complex{S})
    ny,nx = size(model)
    dx = 1/max(nx-1,1)
    dy = 1/max(ny-1,1)
    startx = -0.5
    starty = -0.5
    upx = u*dx
    vpx = v*dy
    phasecenter = exp(2im*π*(u*startx + v*starty))
    @inbounds for i in axes(model.coeff,2), j in axes(model.coeff,1)
        sum += model.coeff[j,i]*exp(2im*π*(upx*(i-1) + vpx*(j-1)))
    end
    return sum*dx*dy*ω(model.kernel, u*dx)*ω(model.kernel, v*dy)*phasecenter
end
