"""
    $(TYPEDEF)

Stores all the non-visibility data products for an EHT array. This is useful when evaluating
model visibilities.

# Fields
$(FIELDS)
"""
Base.@kwdef struct EHTArrayConfiguration{F,T,S,D<:AbstractArray} <: ArrayConfiguration
    """
    Observing bandwith (Hz)
    """
    bandwidth::F
    """
    Telescope array file
    """
    tarr::T
    """
    scan times of the observation
    """
    scans::S
    """
    Time zone used.
    """
    timetype::Symbol = :UTC
    """
    A struct array of `ArrayBaselineDatum` holding time, freq, u, v, baselines.
    """
    data::D
end


Tables.istable(::Type{<:EHTArrayConfiguration}) = true
Tables.columnaccess(::Type{<:EHTArrayConfiguration}) = true
Tables.columns(t::EHTArrayConfiguration) = getfield(t, :data)

Tables.getcolumn(t::EHTArrayConfiguration, ::Type{T}, col::Int, nm::Symbol) where {T} = data(t, nm)
Tables.getcolumn(t::EHTArrayConfiguration, nm::Symbol) = data(t, nm)
Tables.getcolumn(t::EHTArrayConfiguration, i::Int) = Tables.getcolumn(t, Tables.columnames(t)[i])
Tables.columnnames(t::EHTArrayConfiguration) = propertynames(getfield(t, :data))

Base.getindex(data::EHTArrayConfiguration, s::Symbol) = Tables.getcolumn(data, s)
Base.getindex(data::EHTArrayConfiguration, i::Int) = data.data[i]
Base.getindex(data::EHTArrayConfiguration, I...) = getindex(data.data, I...)
Base.length(data::EHTArrayConfiguration) = length(data.data)


"""
    $(SIGNATURES)

Get the u, v, time, freq of the array as a tuple.
"""
function getuvtimefreq(ac::EHTArrayConfiguration)
    u,v = getuv(ac)
    t = ac[:T]
    ν = ac[:F]
    return (U=u, V=v, T=t, F=ν)
end


function ehtarrayconfig(obs::EHTIMObs, pol=nothing)
    u = get(obs.obsdata.data, "u")
    v = get(obs.obsdata.data, "v")
    t = get(obs.obsdata.data, "time")
    f = Fill(obs.rf, length(t))
    t1 = get(obs.data, Vector{Symbol}, "t1")
    t2 = get(obs.data, Vector{Symbol}, "t2")
    baseline = tuple.(t1, t2)
    if pol = nothing
        polbasis = Fill((CirBasis(), CirBasis(), length(t)))
    else
        polbasis = Fill(pol, length(t))
    end

    angles = get_fr_angles(obs)
    tarr = make_array_table(obsc)
    obs.add_scans()
    scans = Table((start=obs.scans[:,1], stop=obs.scans[:,2]))


    uvsamples = StructArray{EHTBaselineDatum}(
                                        T=t,
                                        U=u,
                                        V=v,
                                        F = f,
                                        baseline=baseline,
                                        polbasis=polbasis,
                                        elevation = StructArray(angles[1]),
                                        parallactic  = StructArray(angles[2])
                                    )
    return EHTArrayConfiguration(;bandwidth, tarr, scans,
                                  data = uvsamples)
end
