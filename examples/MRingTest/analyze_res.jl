using CSV, DataFrames, StatsPlots, Plots


function load_chains()
    files = filter(endswith(".csv"), readdir("examples/", join=true))

    dfs = CSV.File.(files) .|> DataFrame

    rings = getindex.(dfs, !, :ring)

    return dfs, process_rings.(rings), files

end

function process_rings(ring)
    tmp = strip.(ring, Ref(['(', ')']))
    tmp2 = split.(tmp, ",")
    names = Symbol.(strip.(first.(split.(tmp2[1], "=")), Ref([' ', '\n'])))
    vals = zeros(length(tmp2), length(names))
    for (i,a) in pairs(tmp2)
        vals[i, :] = parse.(Float64, last.(split.(a, "=")))
    end
    df = DataFrame(vals, names)
    df.davg = df.R .* (2 .- df.ψ)
    return df
end

function plot_ravg(dfrings, labels)
    p = plot()
    xlabel!("Average Radius (μas)")
    for i in eachindex(labels)
        @df dfrings[i] density!(p, :davg, label=labels[i])
    end
    vline!([53.5], label="Truth", linestyle=:dash, color=:black)
    return p
end

function plot_ring(dfrings, quant, labels)
    p = plot()
    xlabel!(String(quant))
    for i in eachindex(labels)
        density!(p, dfrings[i][!,quant], label=labels[i])
    end
    return p
end
