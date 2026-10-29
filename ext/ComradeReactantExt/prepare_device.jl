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
    return Comrade.ObservedSkyModel(Comrade.prepare_device(m.f, ex), grid, Reactant.to_rarray(m.metadata))
end

# The cell layout stays on the host, so indexing with it traces to constant gathers.
function Comrade.prepare_device(m::Comrade.ObservedInstrumentModel, ex::ReactantEx)
    return Comrade.ObservedInstrumentModel(_ondevice(m.instrument), m.refbasis, m.layout, m.pol, Comrade.corrected(m))
end

_ondevice(J) = Reactant.to_rarray(J)
_ondevice(J::Comrade._StationGroups) = Comrade._StationGroups(Reactant.to_rarray(J.terms), J.cells)

function Comrade.prepare_device(m::Comrade.ObservedSelfCalInstrumentModel, ex::ReactantEx)
    stacks = map(st -> merge(st, (; eye = Reactant.to_rarray(st.eye))), m.layout.stacks)
    layout = merge(m.layout, (; stacks))
    return Comrade.ObservedSelfCalInstrumentModel(
        Comrade.prepare_device(m.model, ex), m.phase,
        map(Reactant.to_rarray, m.wdata), map(Reactant.to_rarray, m.weight), layout, m.squarings,
        m.keep, m.entry1, m.entry2
    )
end

function Comrade.prepare_device(grid::VLBISkyModels.FourierDualDomain, ex::ReactantEx)
    gimr = @jit identity(grid.imgdomain)
    guvr = Reactant.to_rarray(grid.visdomain)
    algr = if grid.algorithm isa VLBISkyModels.ReactantNUFFTAlg
        grid.algorithm
    else
        VLBISkyModels.ReactantNUFFTAlg(eltype(grid.imgdomain))
    end
    return VLBISkyModels.FourierDualDomain(gimr, guvr, algr)
end
