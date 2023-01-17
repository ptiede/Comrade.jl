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
@recipe function f(image::Union{StokesIntensityMap, IntensityMap}; uvscale=rad2μas)

    #Define some constants
    #Construct the image grid in μas
    xitr, yitr = imagepixels(image)
    x0, x1 = uvscale.(extrema(xitr))
    y0, y1 = uvscale.(extrema(yitr))

    xitr, yitr = uvscale.(values(imagepixels(image)))
    tickfontsize --> 11
    guidefontsize --> 14

    tickfontsize --> 11
    guidefontsize --> 14
    if image isa StokesIntensityMap

        # get the mean linear pol
        maxI = maximum(stokes(image, :I))

        layout --> (2, 2)
        size --> (500*2,400*2)
        @series begin
            subplot := 1
            seriestype := :heatmap
            seriescolor --> :afmhot
            aspect_ratio --> 1
            bar_width --> 0
            xlims --> (x0, x1)
            ylims --> (y0, y1)
            z = Comrade.baseimage(stokes(image, :I))'
            title --> "Stokes I"
            seriestype := :heatmap
            #fontfamily --> "sans serif"
            xflip --> true
            widen := false
            linecolor-->:black
            tick_direction --> :out
            #colorrange-->(0.0, maxI)

            collect(xitr),collect(yitr),z
        end
        @series begin
            subplot := 2
            seriestype := :heatmap
            seriescolor --> :seismic
            aspect_ratio --> 1
            bar_width --> 0
            xlims --> (x0, x1)
            ylims --> (y0, y1)
            z = Comrade.baseimage(stokes(image, :Q))'
            title --> "Stokes Q"
            #fontfamily --> "sans serif"
            xflip --> true
            widen := false
            linecolor-->:black
            tick_direction --> :out
            clims-->(-maxI/2, maxI/2)

            collect(xitr),collect(yitr),z
        end
        @series begin
            subplot := 3
            seriestype := :heatmap
            xaxis --> "ΔRA  (μas)"
            yaxis --> "ΔDEC (μas)"
            seriescolor --> :seismic
            aspect_ratio --> 1
            bar_width --> 0
            xlims --> (x0, x1)
            ylims --> (y0, y1)
            z = Comrade.baseimage(stokes(image, :U))'
            title --> "Stokes U"
            seriestype := :heatmap
            #fontfamily --> "sans serif"
            xflip --> true
            widen := false
            linecolor-->:black
            tick_direction --> :out
            clims-->(-maxI/2, maxI/2)

            collect(xitr),collect(yitr),z
        end
        @series begin
            subplot := 4
            seriestype := :heatmap
            xaxis --> "ΔRA  (μas)"
            yaxis --> "ΔDEC (μas)"
            seriescolor --> :seismic
            aspect_ratio --> 1
            bar_width --> 0
            xlims --> (x0, x1)
            ylims --> (y0, y1)
            z = Comrade.baseimage(stokes(image, :V))'
            title --> "Stokes V"
            seriestype := :heatmap
            #fontfamily --> "sans serif"
            colorbar_title --> "Jy/px²"
            xflip --> true
            widen := false
            linecolor-->:black
            tick_direction --> :out
            clims-->(-maxI/2, maxI/2)

            collect(xitr),collect(yitr),z
        end
    else
        seriestype := :heatmap
        xaxis --> "ΔRA  (μas)"
        yaxis --> "ΔDEC (μas)"
        seriescolor --> :afmhot
        aspect_ratio --> 1
        bar_width --> 0
        xlims --> (x0, x1)
        ylims --> (y0, y1)
        z = Comrade.baseimage(image)'
        title --> "Stokes I"
        seriestype := :heatmap
        #fontfamily --> "sans serif"
        colorbar_title --> "Jy/px²"
        xflip --> true
        widen := false
        linecolor-->:black
        tick_direction --> :out

        collect(xitr),collect(yitr),z
    end
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

    image = intensitymap(m, fovx, fovy, dims[1], dims[2], phasecenter[1], phasecenter[2])
    xitr, yitr = values(imagepixels(image))
    x0, x1 = uvscale.(extrema(xitr))
    y0, y1 = uvscale.(extrema(yitr))

    #Define some constants
    #Construct the image grid in μas
    fovx, fovy = uvscale.(values(fieldofview(image)))
    xitr, yitr = uvscale.(values(imagepixels(image)))
    tickfontsize --> 11
    guidefontsize --> 14
    if ispolarized(typeof(m)) === IsPolarized()

        # get the mean linear pol
        maxI = maximum(stokes(image, :I))

        layout --> (2, 2)
        size --> (500*2,400*2)
        @series begin
            subplot := 1
            seriestype := :heatmap
            seriescolor --> :afmhot
            aspect_ratio --> 1
            bar_width --> 0
            xlims --> (x0, x1)
            ylims --> (y0, y1)
            z = ComradeBase.baseimage(stokes(image, :I))'
            title --> "Stokes I"
            seriestype := :heatmap
            #fontfamily --> "sans serif"
            xflip --> true
            widen := false
            linecolor-->:black
            tick_direction --> :out
            colorrange-->(0.0, maxI)

            collect(xitr),collect(yitr),z
        end
        @series begin
            subplot := 2
            seriestype := :heatmap
            seriescolor --> :afmhot
            aspect_ratio --> 1
            bar_width --> 0
            xlims --> (x0, x1)
            ylims --> (y0, y1)
            z = ComradeBase.baseimage(stokes(image, :Q))'
            title --> "Stokes Q"
            seriestype := :heatmap
            #fontfamily --> "sans serif"
            xflip --> true
            widen := false
            linecolor-->:black
            tick_direction --> :out
            clims-->(-maxI/2, maxI/2)

            collect(xitr),collect(yitr),z
        end
        @series begin
            subplot := 3
            seriestype := :heatmap
            xaxis --> "ΔRA  (μas)"
            yaxis --> "ΔDEC (μas)"
            seriescolor --> :afmhot
            aspect_ratio --> 1
            bar_width --> 0
            xlims --> (x0, x1)
            ylims --> (y0, y1)
            z = ComradeBase.baseimage(stokes(image, :U))'
            title --> "Stokes U"
            seriestype := :heatmap
            #fontfamily --> "sans serif"
            xflip --> true
            widen := false
            linecolor-->:black
            tick_direction --> :out
            clims-->(-maxI/2, maxI/2)

            collect(xitr),collect(yitr),z
        end
        @series begin
            subplot := 4
            seriestype := :heatmap
            xaxis --> "ΔRA  (μas)"
            yaxis --> "ΔDEC (μas)"
            seriescolor --> :afmhot
            aspect_ratio --> 1
            bar_width --> 0
            xlims --> (x0, x1)
            ylims --> (y0, y1)
            z = ComradeBase.baseimage(stokes(image, :V))'
            title --> "Stokes V"
            seriestype := :heatmap
            #fontfamily --> "sans serif"
            colorbar_title --> "Jy/px²"
            xflip --> true
            widen := false
            linecolor-->:black
            tick_direction --> :out
            clims-->(-maxI/2, maxI/2)

            collect(xitr),collect(yitr),z
        end
    else
        seriestype := :heatmap
        xaxis --> "ΔRA  (μas)"
        yaxis --> "ΔDEC (μas)"
        seriescolor --> :afmhot
        aspect_ratio --> 1
        bar_width --> 0
        xlims --> (x0, x1)
        ylims --> (y0, y1)
        z = ComradeBase.AxisKeys.keyless_unname(image)'
        title --> "Stokes I"
        seriestype := :heatmap
        #fontfamily --> "sans serif"
        colorbar_title --> "Jy/px²"
        xflip --> true
        widen := false
        linecolor-->:black
        tick_direction --> :out

        collect(xitr),collect(yitr),z
    end
end
