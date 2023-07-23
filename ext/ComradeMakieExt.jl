module ComradeMakieExt

using Comrade
if isdefined(Base, :get_extension)
    using Makie
else
    using ..Makie
end

function Makie.convert_arguments(::SurfaceLike, img::IntensityMap{T, 2}) where {T}
    (;X, Y) = img
    return rad2μas(X), rad2μas(Y), Comrade.baseimage(img)
end

function Makie.convert_arguments(::SurfaceLike, x::AbstractVector, y::AbstractVector, m::Comrade.AbstractModel)
    img = intensitymap(m, GriddedKeys((X=x, Y=y)))
    return x, y, Comrade.baseimage(img)
end

@Makie.recipe(ImDisplay, img) do scene
    Makie.Attributes(;
        # Unit conversion
        unit_conv   = Comrade.rad2μas,

        # Plotting stuff
        color         = :black,
        colormap      = :magma,
        dynamic_range = 100,
        alpha         = 1.0,
        pointsize     = 12,
        segmentsize   = 1.5
    )
end

function Makie.plot!(plot::ImDisplay)
    img = plot[:img][]
    bimg = Comrade.baseimage(img)
    (;X, Y) = img
    uc = plot[:unit_conv][]
    cmap = plot[:colormap]
    alpha= plot[:alpha]
    dr = plot[:dynamic_range]
    cmax = maximum(bimg)
    colorrange = (Observable(cmax/dr[]), Observable(cmax))
    Makie.image!(plot, uc.(X), uc.(Y), bimg; colorrange, colormap=cmap, alpha=alpha)
    plot
end

const PolImDisplay = ImDisplay{Tuple{<:Union{<:IntensityMap{<:StokesParams}}}}
function Makie.plot!(plot::PolImDisplay)
    img = plot[:img][]
    bimg = Comrade.baseimage(img)
    (;X, Y) = stokes(img, :I)
    uc = plot[:unit_conv][]
    cmap = plot[:colormap]
    alpha= plot[:alpha]
    dr = plot[:dynamic_range]
    cmax = maximum(bimg)
    colorrange = (Observable(cmax/dr[]), Observable(cmax))
    Makie.image!(plot, uc.(X), uc.(Y), bimg; colorrange, colormap=cmap, alpha=alpha)

    # Now we plot the ticks

    plot
end


end
