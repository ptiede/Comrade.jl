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
    Comrade._drop_host_mirror!(dev.V)
    copyto!(dev.b, h.b)
    copyto!(dev.d, h.d)
    copyto!(dev.V, V)
    copyto!(dev.s, s)
    return dev
end

Comrade._device_space(pre::Comrade.LowRankPreconditioner) =
    Comrade._devicebuffers(pre) ? pre : Comrade._device_pre(pre)
Comrade._device_space(p::Comrade.Preconditioned) =
    Comrade.Preconditioned(p.space, Comrade._device_space(p.pre))

# `tpost` sampling through device buffers of its preconditioner.
function _device_transport(tpost)
    pre = Comrade._transport_pre(tpost)
    (isnothing(pre) || Comrade._devicebuffers(pre)) && return tpost
    space = Comrade._in_space(Comrade._base_space(tpost), Comrade._device_pre(pre))
    return Comrade.maybe_transport(tpost.lpost, space)
end
