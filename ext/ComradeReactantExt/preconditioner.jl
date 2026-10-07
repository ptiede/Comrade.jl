# `V` with zero columns appended up to `cap`.
function _padded(V::AbstractMatrix, cap::Int)
    P = zeros(size(V, 1), cap)
    P[:, 1:size(V, 2)] = V
    return P
end
_padded(V::Comrade.RowSupportedMatrix, cap::Int) = Comrade.RowSupportedMatrix(V.n, V.rows, _padded(V.M, cap))

_on_device(V::AbstractMatrix) = Reactant.to_rarray(V)
_on_device(V::Comrade.RowSupportedMatrix) = Comrade.RowSupportedMatrix(V.n, V.rows, Reactant.to_rarray(V.M))

_copy_columns!(dev::Reactant.AbstractConcreteArray, h::Matrix) = copyto!(dev, h)
function _copy_columns!(dev::Comrade.RowSupportedMatrix, h::Comrade.RowSupportedMatrix)
    dev.rows == h.rows || throw(
        ArgumentError("a refit changed the rows of the preconditioner directions; they must stay fixed")
    )
    copyto!(dev.M, h.M)
    return dev
end
_copy_columns!(dev, h) = throw(
    ArgumentError("a refit changed the storage of the preconditioner directions ($(nameof(typeof(dev))) to $(nameof(typeof(h))))")
)

function Comrade._device_pre(
        pre::Comrade.LowRankPreconditioner; rank_cap::Int = max(size(pre.V, 2), 1)
    )
    m = size(pre.V, 2)
    m <= rank_cap || throw(ArgumentError("rank_cap $rank_cap below fitted rank $m"))
    s = ones(rank_cap)
    s[1:m] = pre.s
    return Comrade.LowRankPreconditioner(
        Reactant.to_rarray(copy(pre.b)), Reactant.to_rarray(copy(pre.d)),
        _on_device(_padded(pre.V, rank_cap)), Reactant.to_rarray(s)
    )
end

function Comrade._update_device_pre!(
        dev::Comrade.LowRankPreconditioner, h::Comrade.LowRankPreconditioner
    )
    cap = size(dev.V, 2)
    m = size(h.V, 2)
    m <= cap || throw(ArgumentError("refit rank $m exceeds the device rank cap $cap"))
    s = ones(cap)
    s[1:m] = h.s
    Comrade._drop_host_mirror!(dev.V)
    copyto!(dev.b, h.b)
    copyto!(dev.d, h.d)
    _copy_columns!(dev.V, _padded(h.V, cap))
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
