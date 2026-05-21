module ComradeReactantExt
using Comrade
using ComradeBase: ReactantEx
using Reactant
using StructArrays
using StaticArraysCore
using Accessors
using PolarizedTypes

@inline function Comrade._apply_instrument!(
        vout::Union{Reactant.AnyTracedRArray, AbstractArray{<:SMatrix{2, 2, <:Reactant.TracedRNumber}}},
        vis,
        J::Comrade.ObservedInstrumentModel,
        xint
    )
    Reactant.@allowscalar @trace track_numbers = false for i in eachindex(vis)
        vout[i] = Comrade.apply_jones(vis[i], i, J, xint)
    end
    return vout
end

function StructArrays.createinstance(::Type{<:StokesParams}, args...)
    return StokesParams(args...)
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

function Comrade.prepare_device(post::VLBIPosterior, ex::ReactantEx)
    ## TODO add Sharding and other options here
    for p in propertynames(post)
        p == :prior && continue
        post = Accessors.set(post, PropertyLens(p), Comrade.prepare_device(getproperty(post, p), ex))
    end
    return post
end

function Comrade.prepare_device(m, ex::ReactantEx)
    return Reactant.to_rarray(m)
end

function Comrade.prepare_device(m::Comrade.ObservedSkyModel, ex::ReactantEx)
    grid = Comrade.prepare_device(m.grid, ex)
    return Comrade.ObservedSkyModel(m.f, grid, Reactant.to_rarray(m.metadata))
end

function Comrade.prepare_device(m::Comrade.ObservedInstrumentModel, ex::ReactantEx)
    return Comrade.ObservedInstrumentModel(Reactant.to_rarray(m.instrument), m.refbasis, m.bsitelookup)
end

function Comrade.prepare_device(grid::VLBISkyModels.FourierDualDomain, ex::ReactantEx)
    gimr = @jit identity(grid.imgdomain)
    guvr = Reactant.to_rarray(grid.visdomain)
    if grid.algorithm isa VLBISkyModels.ReactantNUFFTAlg
        algr = grid.algorithm
    else
        algr = VLBISkyModels.ReactantNUFFTAlg(eltype(grid.imgdomain); eps=1e-16)
    end
    return VLBISkyModels.FourierDualDomain(gimr, guvr, algr)
end
end
