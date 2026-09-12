function Comrade._device_pre(
        pre::Comrade.LowRankPreconditioner; rank_cap::Int = max(size(pre.V, 2), 1)
    )
    n, m = size(pre.V)
    m <= rank_cap || throw(ArgumentError("rank_cap $rank_cap below fitted rank $m"))
    V = zeros(n, rank_cap)
    V[:, 1:m] = pre.V
    s = ones(rank_cap)
    s[1:m] = pre.s
    return Comrade.LowRankPreconditioner(
        Reactant.to_rarray(copy(pre.b)), Reactant.to_rarray(copy(pre.d)),
        Reactant.to_rarray(V), Reactant.to_rarray(s)
    )
end

function Comrade._update_device_pre!(
        dev::Comrade.LowRankPreconditioner, h::Comrade.LowRankPreconditioner
    )
    n, cap = size(dev.V)
    m = size(h.V, 2)
    m <= cap || throw(ArgumentError("refit rank $m exceeds the device rank cap $cap"))
    V = zeros(n, cap)
    V[:, 1:m] = h.V
    s = ones(cap)
    s[1:m] = h.s
    copyto!(dev.b, h.b)
    copyto!(dev.d, h.d)
    copyto!(dev.V, V)
    copyto!(dev.s, s)
    return dev
end
