export dirty_image, dirty_beam

function reflect_vis(obs::EHTObservationTable{<:EHTVisibilityDatum})
    ac = arrayconfig(obs)
    con  = datatable(arrayconfig(obs))
    conn = copy(con)
    conn.U .= -con.U
    conn.V .= -con.V
    con2 = vcat(con, conn)

    meas  = measurement(obs)
    noise = Comrade.noise(obs)
    measn = conj.(copy(meas))

    meas2 = vcat(meas, measn)
    noise2= vcat(noise, noise)

    conf2 = EHTArrayConfiguration(
                ac.bandwidth, ac.tarr, ac.scans,
                ac.mjd, ac.ra, ac.dec, ac.source,
                ac.timetype, con2)

    return EHTObservationTable{datumtype(obs)}(meas2, noise2, conf2)
end


"""
    dirty_image(fov::Real, npix::Int, obs::EHTObservation{<:EHTVisibilityDatum}) where T

Computes the dirty image of the complex visibilities assuming a field of view of `fov`
and number of pixels `npix` using the complex visibilities found in the observation `obs`.

The `dirty image` is the inverse Fourier transform of the measured visibilties assuming every
other visibility is zero.

"""
function dirty_image(fov::Real, npix::Int, obs::EHTObservationTable{D}) where {D<:EHTVisibilityDatum}
    # First we double the baselines, i.e. we reflect them and conjugate the measurements
    # This ensures a real NFFT
    img = IntensityMap(zeros(npix, npix), imagepixels(fov, fov, npix, npix))
    vis2 = reflect_vis(obs)

    # Get the number of pixels
    gfour = FourierDualDomain(axisdims(img), arrayconfig(vis2), NFFTAlg())
    plr = VLBISkyModels.reverse_plan(gfour)
    m = plr.plan*(conj.(vis2[:measurement]).*plr.phases)
    return IntensityMap(real.(m)./npix^2, axisdims(img))
end


"""
    dirty_beam(fov::Real, npix::Int, obs::EHTObservation{<:EHTVisibilityDatum})

Computes the dirty beam of the complex visibilities assuming a field of view of `fov`
and number of pixels `npix` using baseline coverage found in `obs`.

The `dirty beam` is the inverse Fourier transform of the (u,v) coverage assuming every
visibility is unity and everywhere else is zero.

"""
function dirty_beam(fov, npix, obs::EHTObservationTable{D}) where {D<:EHTVisibilityDatum}
    vis2 = reflect_vis(obs)
    vis2.measurement .= complex(oneunit(fov), zero(fov))
    return dirty_image(fov, npix, vis2)
end
