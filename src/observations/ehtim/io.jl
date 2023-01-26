"""
    load_ehtim_uvfits(uvfile, arrayfile=nothing; kwargs...)

Load a uvfits file with eht-imaging and returns a eht-imaging `Obsdata`
object. You can optionally pass an array file as well that will load
additional information such at the telescopes field rotation information
with the arrayfile. This is expected to be an eht-imaging produced array
or antenna file.
"""
function load_ehtim_uvfits(uvfile, arrayfile=nothing; kwargs...)
    ehtim isa PyNULL && load_ehtim()
    obs = ehtim.obsdata.load_uvfits(uvfile; kwargs...)
    if arrayfile !== nothing
        tarr = ehtim.io.load.load_array_txt(arrayfile).tarr
        obs.tarr = tarr
    end
    return obs
end
