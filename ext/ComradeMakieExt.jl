module ComradeMakieExt

using Comrade
if isdefined(Base, :get_extension)
    using Makie
else
    using ..Makie
end

function Makie.convert_arguments(::SurfaceLike, img::IntensityMap{T, 2}) where {T}
    (;X, Y) = img
    return rad2μas(X), rad2μas(Y), Comrade.baseimage(img)./rad2μas(step(X))./rad2μas(step(Y))
end

function Makie.convert_arguments(::SurfaceLike, x::AbstractVector, y::AbstractVector, m::Comrade.AbstractModel)
    img = intensitymap(m, GriddedKeys((X=x, Y=y)))
    return x, y, Comrade.baseimage(img)./rad2μas(step(X))./rad2μas(step(Y))
end


end
