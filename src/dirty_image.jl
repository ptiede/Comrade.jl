export dirty_image, dirty_beam

"""
    dirty_image(fov::Real, npix::Int, t::VisTable)

Computes the dirty image of the complex visibilities in the table `t` assuming a field of
view of `fov` and number of pixels `npix`.

The `dirty image` is the inverse Fourier transform of the measured visibilties assuming every
other visibility is zero.
"""
function dirty_image(fov::Real, npix::Int, t::VisTable)
    keep = _dirtyvalues(t)
    return _dirty(fov, npix, t, keep, baseimage(measurement(t))[keep])
end

"""
    dirty_beam(fov::Real, npix::Int, t::VisTable)

Computes the dirty beam of the baseline coverage of the table `t` assuming a field of view of
`fov` and number of pixels `npix`.

The `dirty beam` is the inverse Fourier transform of the (u,v) coverage assuming every
visibility is unity and everywhere else is zero.
"""
function dirty_beam(fov::Real, npix::Int, t::VisTable)
    keep = _dirtyvalues(t)
    return _dirty(fov, npix, t, keep, fill(complex(oneunit(fov), zero(fov)), length(keep)))
end

function _dirtyvalues(t::VisTable)
    eltype(measurement(t)) <: Complex ||
        throw(ArgumentError("a dirty image is formed from a table of complex visibilities, got one of $(eltype(measurement(t)))"))
    return _measuredvalues(t)
end

# Each baseline and its reflection, with the conjugate visibility, so the transform is real.
function _dirty(fov, npix, t, keep, vis)
    d = domain(t)
    U, V = d.U[keep], d.V[keep]
    (; Ti, Fr) = _measuredpoints(t)
    points = UnstructuredDomain((U = vcat(U, -U), V = vcat(V, -V), Ti = vcat(Ti, Ti), Fr = vcat(Fr, Fr)))
    g = spatialgrid(fov, fov, npix, npix)
    plr = VLBISkyModels.reverse_plan(FourierDualDomain(g, points, NFFTAlg()))
    m = plr.plan * (conj.(vcat(vis, conj.(vis))) .* plr.phases)./npix^2
    return IntensityMap(real.(m), g)
end
