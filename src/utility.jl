export centroid_mean, center_image, convolve!, convolve

function centroid_mean(imgs::AbstractVector{<:IntensityMap})
    return mapreduce(+, imgs) do img
        center_image(img)
    end
    return mimg./length(imgs)
end

"""
    center_image(img::IntensityMap)

centers the `img` such that the centroid of the image is approximately at the origin.
"""
function center_image(img::IntensityMap)
    x, y = centroid(img)
    return modify(img, Shift(-x, -y))
end

# This is an internal struct that is use to modify IntensityMaps so that we can hook into
# Comrade image modifier interface.
struct InterpolatedImage{I,P} <: AbstractModel
    img::I
    itp::P
    function InterpolatedImage(img::IntensityMap)
        itp = BicubicInterpolator(img.X, img.Y, img, StrictBoundaries())
        return new{typeof(img), typeof(itp)}(img, itp)
    end
end


imanalytic(::Type{<:InterpolatedImage})  = IsAnalytic()
visanalytic(::Type{<:InterpolatedImage}) = NotAnalytic()

function intensity_point(m::InterpolatedImage, p)
    (m.img.X[begin] > p.X || p.X > m.img.X[end]) && return zero(p.X)
    (m.img.Y[begin] > p.Y || p.Y > m.img.Y[end]) && return zero(p.X)
    return m.itp(p.X, p.Y)/(step(m.img.X)*step(m.img.Y))
end

"""
    modify(img::IntensityMap, transforms...)

This modifies the `img` by applying the `transforms...` returning a transformed `IntensityMap`

!!! note
Unlike when `modify` is applied to a `<:AbstractModel` this returns an already modified image.
"""
function ModifiedModel(img::IntensityMap, transforms)
    ms = ModifiedModel(InterpolatedImage(img), transforms)
    return intensitymap(ms, axiskeys(img))
end
modify(img::IntensityMap, transforms...) = ModifiedModel(img, transforms)


"""
    convolve!(img::IntensityMap, m::AbstractModel)

Convolves an `img` with a given analytic model `m`. This is useful for blurring the
image with some model. For instance to convolve a image with a Gaussian you would do
```julia
convolve!(img, Gaussian())
```

# Notes
This method does not automatically pad your image. If there is substantial flux at the boundaries
you will start to see artifacts.
"""
function convolve!(img::IntensityMap{<:Real}, m::AbstractModel)
    @assert visanalytic(typeof(m)) isa IsAnalytic "Convolving model must have an analytic Fourier transform currently"
    p = plan_rfft(baseimage(img))

    (;X, Y) = imagepixels(img)
    # plan_rfft uses just the positive first axis to respect real conjugate symmetry
    U = rfftfreq(size(img, 1), inv(step(X)))
    V = fftfreq(size(img, 2), inv(step(Y)))

    # TODO maybe ask a user to pass a vis buffer as well?
    vis = p*baseimage(img)

    # Conjugate because Comrade uses +2πi exponent
    vis .*= conj(visibility_point.(Ref(m), U, V', 0, 0))
    pinv = plan_irfft(vis, size(img, 1))
    mul!(baseimage(img), pinv, vis)
    return img
end

"""
    convolve(img::IntensityMap, m::AbstractModel)

Convolves an `img` with a given analytic model `m`. This is useful for blurring the
image with some model. For instance to convolve a image with a Gaussian you would do
```julia
convolve(img, Gaussian())
```

For the inplace version of the function see [`convolve!`](@ref)

# Notes
This method does not automatically pad your image. If there is substantial flux at the boundaries
you will start to see artifacts.
"""
function convolve(img::IntensityMap, m::AbstractModel)
    cimg = copy(img)
    return convolve!(cimg, m)
end

function convolve!(img::IntensityMap{<:StokesParams}, m)
    convolve!(stokes(img, :I), m)
    convolve!(stokes(img, :Q), m)
    convolve!(stokes(img, :U), m)
    convolve!(stokes(img, :V), m)
    return img
end




export dirty_image, dirty_beam
function reflect_vis(obs)
    meas = obs.data
    con  = copy(obs.config.data)
    conn = copy(con)
    conn.U .= -con.U
    conn.V .= -con.V
    con2 = vcat(con, conn)

    meas  = copy(obs.data)
    measn = copy(meas)
    measn.U .= -meas.U
    measn.V .= -meas.V
    measn.measurement .= conj.(meas.measurement)

    meas2 = vcat(meas, measn)

    conf2 = EHTArrayConfiguration(obs.config.bandwidth, obs.config.tarr, obs.config.scans, con2)
    vis2  = EHTObservation(; data = meas2, config=conf2,
                             mjd=obs.mjd, ra=obs.ra, dec=obs.dec,
                             bandwidth=obs.bandwidth, source=obs.source,
                             timetype=obs.timetype)
    return vis2

end


"""
    dirty_image(fov::Real, npix::Int, obs::EHTObservation{T,<:EHTVisibilityDatum}) where T

Computes the dirty image of the complex visibilities assuming a field of view of `fov`
and number of pixels `npix` using the complex visibilities found in the observation `obs`.

The `dirty image` is the inverse Fourier transform of the measured visibilties assuming every
other visibility is zero.

"""
function dirty_image(fov::Real, npix::Int, obs::EHTObservation{T,D}) where {T, D<:EHTVisibilityDatum}
    # First we double the baselines, i.e. we reflect them and conjugate the measurements
    # This ensures a real NFFT
    vis2 = reflect_vis(obs)

    # Get the number of pixels
    alg = NFFTAlg(vis2)
    img = IntensityMap(zeros(npix, npix), imagepixels(fov, fov, npix, npix))
    cache = create_cache(alg, img)
    m = cache.plan'*conj.(vis2[:measurement])
    return IntensityMap(real.(m)./npix^2, imagepixels(fov, fov, npix, npix))
end


"""
    dirty_beam(fov::Real, npix::Int, obs::EHTObservation{T,<:EHTVisibilityDatum}) where T

Computes the dirty beam of the complex visibilities assuming a field of view of `fov`
and number of pixels `npix` using baseline coverage found in `obs`.

The `dirty beam` is the inverse Fourier transform of the (u,v) coverage assuming every
visibility is unity and everywhere else is zero.

"""
function dirty_beam(fov, npix, obs::EHTObservation{T,D}) where {T, D<:EHTVisibilityDatum}
    vis2 = reflect_vis(obs)
    vis2.data.measurement .= complex(1.0, 0.0)
    return dirty_image(fov, npix, vis2)
end
