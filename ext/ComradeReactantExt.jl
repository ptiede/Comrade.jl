module ComradeReactantExt
using Comrade
using Reactant

@inline function Comrade._apply_instrument!(
    vout::Reactant.AnyTracedRArray,
    vis,
    J::Comrade.ObservedInstrumentModel, 
    xint
)
    Reactant.@allowscalar @trace track_numbers=false for i in eachindex(vis)
        vout[i] = Comrade.apply_jones(vis[i], i, J, xint)
    end
end

@inline function Comrade.site_sum(y::Reactant.AnyTracedRArray, site_map::Comrade.SiteLookup)
    yout = similar(y)
    vals = values(lookup(site_map))
    #unroll this?
    for site in vals
        ys = y[site]
        yout[site] = cumsum(ys)
    end
    return yout
end

function Comrade.branchcut(x::Reactant.TracedRNumber)
    xmod = mod(x, oftype(x, 2π))
    return ifelse(xmod > π, xmod - oftype(x, 2π), xmod)
end




end