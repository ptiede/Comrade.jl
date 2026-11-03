module ComradeMakieExt

using Makie
using Comrade

function Comrade.residualplot(post::Comrade.AbstractVLBIPosterior, p; figure = (;), kwargs...)
    res = Comrade.residuals(post, p)
    fig = Figure(; size = (600, 350 * length(res)), figure...)
    for (i, r) in enumerate(res)
        fieldplot(fig[i, 1], r, :uvdist, :res; kwargs...)
    end
    return fig
end

end
