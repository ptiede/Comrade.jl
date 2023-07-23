export dirty_image, dirty_beam, load_clean_components, MultiComponentModel

using DelimitedFiles

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
    img = IntensityMap(zeros(npix, npix), imagepixels(fov, fov, npix, npix))
    dx, dy = pixelsizes(img)
    vis2 = reflect_vis(obs)

    # Get the number of pixels
    alg = NFFTAlg(vis2)
    cache = create_cache(alg, img)
    u = vis2[:U]
    v = vis2[:V]
    m = cache.plan'*(conj.(vis2[:measurement].*cispi.(-(u.*(dx) .+ v.*(dy)))))
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
    vis2.data.measurement .= complex(one(T), zero(T))
    return dirty_image(fov, npix, vis2)
end
