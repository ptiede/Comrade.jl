using RecipesBase
using Printf


"""
    plot(image::IntensityMap)

where `image` is templated off of EHTImage struct.

# Details
This was created to be close to the ehtim display object. It takes an
EHTImage object and plots it according to EHT conventions.

Note that is does not save the figure.
"""
@recipe function f(image::IntensityMap; uvscale=rad2μas)

    #Define some constants
    #Construct the image grid in μas
    fovx, fovy = uvscale.(fov(image))
    xitr, yitr = imagepixels(image)
    x0, x1 = uvscale.(extrema(xitr))
    y0, y1 = uvscale.(extrema(yitr))

    tickfontsize --> 11
    guidefontsize --> 14
    size --> (500,400)
    xaxis --> "ΔRA  (μas)"
    yaxis --> "ΔDEC (μas)"
    seriescolor --> :afmhot
    aspect_ratio --> 1
    bar_width --> 0
    xlims --> (x0, x1)
    ylims --> (y0, y1)
    #left_margin --> -2mm
    #right_margin --> 5
    z = image
    seriestype := :heatmap
    #fontfamily --> "sans serif"
    colorbar_title --> "Jy/px"
    xflip --> true
    widen := false
    framestyle --> :box
    title --> @sprintf("flux = %.2f Jy", flux(image))
    linecolor-->:black
    tick_direction --> :out
    uvscale.(collect(xitr)),uvscale.(collect(yitr)),z
end


"""
    plot(image::AbstractModel)

where `image` is templated off of EHTImage struct.

# Details
This was created to be close to the ehtim display object. It takes an
EHTImage object and plots it according to EHT conventions.

Note that is does not save the figure.
"""
@recipe function f(m::AbstractModel; uvscale=rad2μas,
                   fovx = 2*radialextent(m), fovy=2*radialextent(m), dims=(512, 512),
                   phasecenter = (0.0, 0.0), pulse=DeltaPulse())

    ny, nx = dims
    image = intensitymap(m, (fovx, fovy), dims; phasecenter, pulse)
    xitr, yitr = imagepixels(image)
    psizex,psizey = pixelsizes(image)
    x0, x1 = uvscale.(extrema(xitr))
    y0, y1 = uvscale.(extrema(yitr))

    #Define some constants
    #Construct the image grid in μas
    fovx, fovy = uvscale.(fov(image))
    xitr, yitr = uvscale.(imagepixels(image))
    tickfontsize --> 11
    guidefontsize --> 14
    size --> (500,400)
    xaxis --> "ΔRA  (μas)"
    yaxis --> "ΔDEC (μas)"
    seriescolor --> :afmhot
    aspect_ratio --> 1
    bar_width --> 0
    xlims --> (x0, x1)
    ylims --> (y0, y1)
    #left_margin --> -2mm
    #right_margin --> 5mm
    z = image
    seriestype := :heatmap
    #fontfamily --> "sans serif"
    colorbar_title --> "Jy/px²"
    xflip --> true
    widen := false
    #framestyle --> :box
    title --> @sprintf("flux = %.2f Jy", flux(image))
    linecolor-->:black
    tick_direction --> :out
    collect(xitr),collect(yitr),z
end
