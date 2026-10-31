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

function Comrade.prepare_device(m::Comrade.ObservedInstrumentModel, ex::ReactantEx)
    return Comrade.ObservedInstrumentModel(Reactant.to_rarray(m.instrument), m.refbasis, m.bsitelookup)
end

# The index maps stay on the host so that indexing with them traces to constant gathers.
function Comrade.prepare_device(m::Comrade.ObservedSelfCalInstrumentModel, ex::ReactantEx)
    layout = merge(m.layout, (; eye = Reactant.to_rarray(m.layout.eye)))
    return Comrade.ObservedSelfCalInstrumentModel(
        Comrade.prepare_device(m.model, ex), m.phase,
        map(Reactant.to_rarray, m.wdata), map(Reactant.to_rarray, m.weight), layout, m.iterations
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
