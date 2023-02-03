# # Self-calibrating data with Comrade

# Here we will self-cal the data with Comrade. This means we take the results of a
# Comrade run, and the use eht-imaging to self-calibrate the fitted data and save the
# new observation object.

# This script assumes that you have the eht-imaging Obsdata object, the
# Comrade MCMC chain, and the model you fit to the data.

using Pkg;Pkg.activate(@__DIR__)
using Comrade
using Tables
using DataFrames
using CSV
using PyCall
using Printf
using Serialization

load_ehtim()

@pyimport numpy as np

"""
    selfcal(obs, model::GainModel)

Takes in a `eht-imaging` observation object that you fit the data to, and
a `Comrade` gain model and returns the `self-calibrated` data and the calibration table
object from eht-imaging.
"""
function selfcal(obs::PyCall.PyObject, model::GainModel)
    obscpy = obs.copy()
    ctable = caltable(model)

    # Now make the caltable for eht-imaging
    sites = stations(ctable)
    time = ctable.time
    ctabcol = map(sites) do s
        mask = Base.:!.(ismissing.(getproperty(ctable, s)))
        tmask = time[mask]
        gmask = ctable[:,s][mask]
        # make the ndarray that eht-imaging expects
        arr = np.array(PyObject(tmask), dtype = ehtim.DTCAL)
        set!(arr, "time", tmask)
        set!(arr, "rscale", inv.(gmask))
        set!(arr, "lscale", inv.(gmask))
        return Pair(s, arr)
    end
    datatable = PyDict(Dict(ctabcol))
    ctab_eh = ehtim.caltable.Caltable(
            obs.ra, obs.dec,
            obs.rf, obs.bw,
            datatable,
            obs.tarr,
            source=obs.source, mjd=obs.mjd)
    o2 = ctab_eh.applycal(obscpy)
    return o2, ctab_eh
end


using Serialization

struct GModel{C,G}
    cache::C
    gcache::G
    fovx::Float64
    fovy::Float64
    npixx::Int
    npixy::Int
    function GModel(obs::Comrade.EHTObservation, fovx, fovy, npixx, npixy)
        buffer = IntensityMap(zeros(npixx, npixy), fovx, fovy, (0.0,0.0), BSplinePulse{3}())
        cache = create_cache(DFTAlg(obs), buffer)
        gcache = GainCache(scantable(obs))
        return new{typeof(cache), typeof(gcache)}(cache, gcache, fovx, fovy, npixx, npixy)
    end
end

function (model::GModel)(θ)
    (;c, f, lgamp) = θ
    # Construct the image model
    img = IntensityMap(f*c, model.fovx, model.fovy, (0.0, 0.0), BSplinePulse{3}())
    m = modelimage(img, model.cache)
    # Now corrupt the model with Gains
    g = exp.(lgamp)
    Comrade.GainModel(model.gcache, g, m)
end

"""
    selfcal_submission(results, obsfile, outdir, nsamples, nburn)

Creates a directory with the selfcalibrated data sets and their corresponding images.

# Arguments
    - `results::String`: The .jls serialized file from a run
    - `obsfile::String`: The data you fit
    - `outdir::String`: The directory you wish to save the images and uvfits files to
    - `nsamples::Int`: The number of samples from the posterior you want to create self-cal data sets for
    - `nburn::Int`: The number of adaptation steps used in the sample. This is required because you need to remove these

# Warning
We have assumed a number of things about how the data was processed, i.e. we cut zero baselines. If this
is not the case then this function will need to be modified.
"""
function selfcal_submission(results::String, obsfile::String, outdir::String, nsamples::Int, nburn::Int)
    res = deserialize(results)
    obs = ehtim.obsdata.load_uvfits(obsfile)
	obs.add_scans()
	# make scan-average data #.flag_uvdist(uv_min=0.1e9)
	obs = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true).add_fractional_noise(0.01)

	# extract amplitudes and closure phases
	damp = extract_amp(obs)

    mms = GModel(damp, res[:fovx], res[:fovy], res[:npixx], res[:npixy])

    mkpath(outdir)

    fov = max(res[:fovx], res[:fovy])*1.1

    chain = sample(res[:chain][nburn:end], nsamples)
    for i in 1:nsamples
        model = mms(chain[i])
        img = intensitymap(model, fov, fov, 256, 256)
        obscal, _ = selfcal(obs, model)

        outim = @sprintf "image_%04d.fits" i
        Comrade.save(joinpath(outdir, outim), img, damp)

        outcal = @sprintf "selfcal_data_comrade_%04d.uvfits" i
        obscal.save_uvfits(joinpath(outdir, outcal))
    end
    return nothing
end
