"""
    DFTAlg(obs::EHTObservation)

Create an algorithm object using the direct Fourier transform object from the observation
`obs`. This will extract the uv positions from the observation to allow for a more efficient
FT cache.
"""
function VLBISkyModels.DFTAlg(obs::EHTObservationTable)
    (;U,V) = getuv(obs.config)
    return DFTAlg(U, V)
end

"""
    DFTAlg(ac::ArrayConfiguration)

Create an algorithm object using the direct Fourier transform object from the array configuration
`ac`. This will extract the uv positions from the observation to allow for a more efficient
FT cache.
"""
function VLBISkyModels.DFTAlg(ac::EHTArrayConfiguration)
    (;U,V) = getuv(ac)
    return DFTAlg(U, V)
end


"""
    NFFTAlg(ac::ArrayConfiguration; kwargs...)

Create an algorithm object using the non-unform Fourier transform object from the array
configuration `ac`. This will extract the uv positions from the observation to allow
for a more efficient FT cache.

The optional arguments are: `padfac` specifies how much to pad the image by, and `m`
is an internal variable for `NFFT.jl`.
"""
function VLBISkyModels.NFFTAlg(ac::EHTArrayConfiguration; kwargs...)
    (;U, V) = getuv(ac)
    return NFFTAlg(U, V; kwargs...)
end


"""
    NFFTAlg(obs::EHTObservation; kwargs...)

Create an algorithm object using the non-unform Fourier transform object from the observation
`obs`. This will extract the uv positions from the observation to allow for a more efficient
FT cache.

The possible optional arguments are given in the [`NFFTAlg`](@ref) struct.
"""
function VLBISkyModels.NFFTAlg(obs::EHTObservationTable; kwargs...)
    (;U, V) = getuv(arrayconfig(obs))
    return NFFTAlg(U, V; kwargs...)
end
