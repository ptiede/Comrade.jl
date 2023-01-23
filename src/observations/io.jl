"""
    Comrade.load(fitsfile::String, IntensityMap)

This loads in a fits file that is more robust to the various imaging algorithms
in the EHT, i.e. is works with clean, smili, eht-imaging.
The function returns an tuple with an intensitymap and a second named tuple with ancillary
information about the image, like the source name, location, mjd, and radio frequency.
"""
function load(file, T::Type{<:IntensityMapTypes})
    if !endswith(file, ".fits")
        @warn "File does not end with FITS trying to load anyways"
    end
    return _load_fits(file, T)
end

function _load_fits(fname, ::Type{IntensityMap})
    img = FITS(fname, "r") do f
        if length(f) > 1
            @warn "Currently only loading stokes I. To load polarized quantities\n"*
                  "please call `Comrade.load(filename, StokesIntensityMap)`"
        end
        # assume that the first element is stokes I
        return _extract_fits_image(f[1])
    end
    return img
end

function try_loading(f, stokes, imgI)
    try
        return _extract_fits_image(f[stokes])
    catch
        @warn "No stokes $(stokes) found creating a zero array"
        imgQ = zeros(imgI)
        return imgQ

    end
end


function _load_fits(fname, ::Type{StokesIntensityMap})
    img = FITS(fname, "r") do f
        # assume that the first element is stokes I
        imgI = _extract_fits_image(f[1])
        imgQ = try_loading(f, "Q", imgI)
        imgU = try_loading(f, "U", imgI)
        imgV = try_loading(f, "V", imgI)

        return StokesIntensityMap(imgI, imgQ, imgU, imgV), head
    end
    return img
end



function _extract_fits_image(f::FITSIO.ImageHDU{T,2}) where {T}
    image = read(f)[end:-1:begin,:]
    header = read_header(f)
    nx = Int(header["NAXIS1"])
    ny = Int(header["NAXIS2"])

    psizex = abs(float(header["CDELT1"]))*π/180
    psizey = abs(float(header["CDELT2"]))*π/180

    ra = float(header["OBSRA"])
    dec = float(header["OBSDEC"])

    #Get frequency
    freq = 0.0
    if haskey(header, "FREQ")
        freq = parse(Float64, string(header["FREQ"]))
    elseif "CRVAL3" in keys(header)
        freq = float(header["CRVAL3"])
    end
    mjd = 0.0
    if haskey(header, "MJD")
        mjd = parse(Float64, string(header["MJD"]))
    end
    source = "NA"
    if haskey(header,"OBJECT")
        source = string(header["OBJECT"])
    end
    stokes = "NA"
    if haskey(header, "STOKES")
        stokes = Symbol(header["STOKES"])
    end
    bmaj = 1.0 #Nominal values
    bmin = 1.0
    if haskey(header, "BUNIT")
        if header["BUNIT"] == "JY/BEAM"
            @info "Converting Jy/Beam => Jy/pixel"
            bmaj = header["BMAJ"]*π/180
            bmin = header["BMIN"]*π/180
            beamarea = (2.0*π*bmaj*bmin)/(8*log(2))
            image .= image.*(psizex*psizey/beamarea)
        end
    end
    info = (source=source, RA=ra, DEC=dec, mjd=mjd, ν=freq, stokes=stokes)
    imap = IntensityMap(image, psizex*nx, psizey*ny; header=info)
    return imap
end

"""
    Comrade.save(file::String, img::IntensityMap, obs)

Saves an image to a fits file. You can optionally pass an EHTObservation so that ancillary information
will be added.
"""
function save(fname::String, img::IntensityMapTypes, obs = nothing)
    head = make_header(obs)
    _save_fits(fname, img, head)
end

function make_header(obs)
    if isnothing(obs)
        return (source="NA", RA=0.0, DEC=0.0, mjd=0, freq=0.0)
    else
        return (source=String(obs.source), RA=obs.ra, DEC=obs.dec, mjd=obs.mjd, freq=first(obs[:F]))
    end
end

function _prepare_header(image, stokes="I")
    head = header(image)
    headerkeys = ["SIMPLE",
                  "BITPIX",
                  "NAXIS",
                  "NAXIS1",
                  "NAXIS2",
                  "EXTEND",
                  "OBJECT",
                  "CTYPE1",
                  "CTYPE2",
                  "CDELT1",
                  "CDELT2",
                  "OBSRA",
                  "OBSDEC",
                  "FREQ",
                  "CRPIX1",
                  "CRPIX2",
                  "MJD",
                  "TELESCOP",
                  "BUNIT",
                  "STOKES"]

    psizex, psizey = pixelsizes(image)
    values = [true,
              -64,
              2,
              size(image, 1),
              size(image, 2),
              true,
              head.source,
              "RA---SIN",
              "DEC---SIN",
              rad2deg(psizex),
              rad2deg(psizey),
              head.RA,
              head.DEC,
              head.freq,
              size(image,1)/2+0.5,
              size(image,2)/2+0.5,
              head.mjd,
              "VLBI",
              "JY/PIXEL",
              stokes]
    comments = ["conforms to FITS standard",
                "array data type",
                "number of array dimensions",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                ""]

    return headerkeys, values, comments
end

function _save_fits(fname::String, image::IntensityMap{T}, head) where {T<:Number}
    FITS(fname, "w") do hdu
        write_stokes(hdu, image)
    end
end

function write_stokes(f, image, stokes="I", innername="")
    headerkeys, values, comments = _prepare_header(image, stokes)
    hdeheader = FITSHeader(headerkeys, values, comments)
    img = ComradeBase.AxisKeys.keyless_unname(image[end:-1:1, :])
    FITSIO.write(f, img; header=hdeheader, name=innername)
end

function _save_fits(fname::String, image::Union{StokesIntensityMap, IntensityMap{T}}) where {T<:StokesParams}
    FITS(fname, "w") do fits
        write_stokes(fits, ComradeBase.stokes(image, :I), "I")
        write_stokes(fits, ComradeBase.stokes(image, :Q), "Q", "Q")
        write_stokes(fits, ComradeBase.stokes(image, :U), "U", "U")
        write_stokes(fits, ComradeBase.stokes(image, :V), "V", "V")
    end
end


"""
    $(SIGNATURES)

Load a ThemisPy style ascii EHT observation file.
"""
function load_tpy(file)
    data = readdlm(file, skipstart=1)
    bs1 = Symbol.(getindex.(data[:,5], Ref(1:2)))
    bs2 = Symbol.(getindex.(data[:,5], Ref(3:4)))
    baselines = tuple.(bs1, bs2)
    edata = StructArray{EHTVisibilityDatum{Float64}}(
                visr=float.(data[:,8]),
                visi=float.(data[:,10]),
                error=float.(data[:,9]),
                u=float.(data[:,6])*1e6,
                v=float.(data[:,7])*1e6,
                time=float.(data[:,4]),
                frequency=fill(227e9, size(data,1)),
                bandwidth=fill(4e6, size(data,1)),
                baseline=baselines
            )
    mjd =  Int(modified_julian(UTCEpoch(Int(data[1,2]), 1,1,0)).Δt+float(data[1,3]))
    return EHTObservation(data=edata, mjd=mjd, ra=180.0, dec=0.0, source=Symbol(data[1,1]))
end
