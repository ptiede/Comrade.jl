"""
    $(TYPEDEF)

This defined the abstract type for an array configuration. Namely, baseline
times, SEFD's, bandwidth, observation frequencies, etc.
"""
abstract type ArrayConfiguration end

scans(d::ArrayConfiguration) = d.scans
telescope_array(d::ArrayConfiguration) = d.tarr
bandwidth(d::ArrayConfiguration) = d.bandwidth

"""
    getuv

Get the u, v positions of the array.
"""
function getuv(ac::ArrayConfiguration)
    return (U=ac.data.U, V=ac.data.V)
end

data(obs::ArrayConfiguration) = obs.data
data(obs::ArrayConfiguration, v::Symbol) = getproperty(obs.data, v)


"""
    $(TYPEDEF)
Array config file for closure quantities. This stores the design matrix `designmat`
that transforms from visibilties to closure products.

# Fields
$(FIELDS)
"""
struct ClosureConfig{A,D} <: ArrayConfiguration
    """Array configuration for visibilities"""
    ac::A
    """Closure design matrix"""
    designmat::D
    function ClosureConfig(ac, dmat)
        A = typeof(ac)
        sdmat = blockdiag(sparse.(dmat)...)
        D = typeof(sdmat)
        return new{A,D}(ac, sdmat)
    end
end

function getuvtimefreq(ac::ClosureConfig)
    return getuvtimefreq(ac.ac)
end
