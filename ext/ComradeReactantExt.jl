module ComradeReactantExt
using Comrade
using Reactant

@inline function Comrade._apply_instrument(
    vis::Reactant.AnyTracedRArray, 
    J::Comrade.ObservedInstrumentModel, 
    xint
)
    Reactant.@allowscalar @trace track_numbers=false for i in eachindex(vis)
        vis[i] = apply_jones(vis[i], i, J, xint)
    end
    return vis
end


end