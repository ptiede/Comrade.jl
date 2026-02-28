module ComradeReactantExt
using Comrade
using Reactant
using StructArrays
using StaticArraysCore

@inline function Comrade._apply_instrument!(
        vout::Union{Reactant.AnyTracedRArray, AbstractArray{<:SMatrix{2, 2, <:Reactant.TracedRNumber}}},
        vis,
        J::Comrade.ObservedInstrumentModel,
        xint
    )
    return Reactant.@allowscalar @trace track_numbers = false for i in eachindex(vis)
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

# @inline function Comrade.site_sum(y::Reactant.AnyTracedRArray, site_map::Comrade.SiteLookup)
#     yout = similar(y)
#     vals = values(lookup(site_map))
#     #unroll this?
#     @trace track_numbers=false for i in eachindex(vals)
#         site = vals[i]
#         ys = y[site]
#         acc = zero(eltype(y))
#         @inbounds ys = @view y[site]
#         for i in eachindex(ys)
#             acc += ys[i]
#             ys[i] = acc
#         end
#         yout[site] = ys
#     end
#     return yout
# end


function Comrade.branchcut(x::Reactant.TracedRNumber)
    xmod = mod(x, oftype(x, 2π))
    return ifelse(xmod > π, xmod - oftype(x, 2π), xmod)
end


end
