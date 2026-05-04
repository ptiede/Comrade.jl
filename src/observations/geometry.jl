export ecef_to_geodetic, elevation_parallactic

# WGS84 ellipsoid constants
const _WGS84_A = 6378137.0                # semi-major axis (m)
const _WGS84_F = 1 / 298.257223563        # flattening
const _WGS84_B = _WGS84_A * (1 - _WGS84_F) # semi-minor axis
const _WGS84_E2  = _WGS84_F * (2 - _WGS84_F)
const _WGS84_EP2 = (_WGS84_A^2 - _WGS84_B^2) / _WGS84_B^2


"""
    ecef_to_geodetic(x, y, z) -> (lat, lon, h)

Convert WGS84 Earth-centered Earth-fixed coordinates `(x, y, z)` in metres to
geodetic latitude (rad), longitude (rad), and ellipsoidal height (m).

Uses the closed-form Bowring formula. Result is accurate to better than a
millimetre for any point on or near the Earth's surface.
"""
function ecef_to_geodetic(x, y, z)
    a = _WGS84_A
    b = _WGS84_B
    e2 = _WGS84_E2
    ep2 = _WGS84_EP2

    p = hypot(x, y)
    θ = atan(z * a, p * b)
    lon = atan(y, x)
    lat = atan(z + ep2 * b * sin(θ)^3, p - e2 * a * cos(θ)^3)
    N = a / sqrt(1 - e2 * sin(lat)^2)
    h = p / cos(lat) - N
    return (lat, lon, h)
end


"""
    elevation_parallactic(antenna_xyz, ra_deg, dec_deg, jd)
        -> (elevation_rad, parallactic_rad)

Compute the source elevation and parallactic angle at a single antenna and time.

Arguments:
- `antenna_xyz`: ECEF coordinates of the antenna in metres (3-tuple, vector, or
  `SVector{3}`).
- `ra_deg`, `dec_deg`: J2000 right ascension and declination of the source in
  degrees.
- `jd`: Julian date (UT) of the observation.

Returns the elevation angle (radians, above horizon) and parallactic angle
(radians, atan2 form) at that antenna.

Uses [`AstroLib.ct2lst`](https://github.com/JuliaAstro/AstroLib.jl) for local
mean sidereal time. The parallactic angle uses the standard atan2 form

    p = atan2(sin(HA) cos(lat), sin(lat) cos(δ) - cos(lat) sin(δ) cos(HA))

which matches the convention used by both eht-imaging and Astropy.
"""
function elevation_parallactic(antenna_xyz, ra_deg, dec_deg, jd;
                                precession::Bool = false,
                                nutate::Bool = false,
                                aberration::Bool = false,
                                refract::Bool = false)
    x, y, z = antenna_xyz[1], antenna_xyz[2], antenna_xyz[3]
    lat, lon, alt_m = ecef_to_geodetic(x, y, z)

    # Defer to AstroLib's eq2hor for the equatorial -> horizon conversion. By
    # default we disable precession/nutation/aberration/refraction so the
    # result matches eht-imaging's `utc_to_gmst`-based pipeline; flip the
    # flags to get full Astropy-style apparent-place corrections.
    alt_deg, _, ha_deg = AstroLib.eq2hor(
        ra_deg, dec_deg, jd,
        rad2deg(lat), rad2deg(lon), alt_m;
        precession, nutate, aberration, refract
    )
    el = deg2rad(alt_deg)
    ha = deg2rad(ha_deg)
    dec = deg2rad(dec_deg)
    par = atan(
        sin(ha) * cos(lat),
        sin(lat) * cos(dec) - cos(lat) * sin(dec) * cos(ha)
    )
    return (el, par)
end
