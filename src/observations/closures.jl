"""
Helpers for constructing closure design matrices and `ClosureConfig`s
from a visibility `EHTObservationTable` plus a closure-spec table.

These helpers are used by `extract_cphase` / `extract_lcamp` in extensions
that load data from external formats (e.g. ehtim, VLBIFiles).
"""


"""
    closure_designmat(type::Symbol, closures, scanvis::Scan)

Dispatch to `closurephase_designmat` or `closureamp_designmat` based on `type`
(`:cphase` or `:lcamp`).
"""
function closure_designmat(type, closures, scanvis)
    return if type == :cphase
        closurephase_designmat(closures, scanvis)
    elseif type == :lcamp
        closureamp_designmat(closures, scanvis)
    else
        throw(ArgumentError("Not a valid type of closure"))
    end
end


function closurephase_designmat(cphase, scanvis)
    antvis1, antvis2 = baseline(scanvis)
    design_mat = zeros(length(cphase), length(scanvis))
    for i in axes(design_mat, 2), j in axes(design_mat, 1)
        a1, a2, a3 = cphase[j].baseline
        # leg 1
        ((antvis1[i] == a1) & (antvis2[i] == a2)) && (design_mat[j, i] = 1.0)
        ((antvis1[i] == a2) & (antvis2[i] == a1)) && (design_mat[j, i] = -1.0)
        # leg 2
        ((antvis1[i] == a2) & (antvis2[i] == a3)) && (design_mat[j, i] = 1.0)
        ((antvis1[i] == a3) & (antvis2[i] == a2)) && (design_mat[j, i] = -1.0)
        # leg 3
        ((antvis1[i] == a3) & (antvis2[i] == a1)) && (design_mat[j, i] = 1.0)
        ((antvis1[i] == a1) & (antvis2[i] == a3)) && (design_mat[j, i] = -1.0)
    end
    return design_mat
end


function closureamp_designmat(lcamp, scanvis)
    antvis1, antvis2 = baseline(scanvis)
    design_mat = zeros(length(lcamp), length(scanvis))
    for i in axes(design_mat, 2), j in axes(design_mat, 1)
        a1, a2, a3, a4 = lcamp[j].baseline

        av1 = antvis1[i]
        av2 = antvis2[i]
        # leg 1
        (((av1 == a1) & (av2 == a2)) || (av1 == a2) & (av2 == a1)) && (design_mat[j, i] = 1.0)
        # leg 2
        (((av1 == a3) & (av2 == a4)) || (av1 == a4) & (av2 == a3)) && (design_mat[j, i] = 1.0)
        # leg 3
        (((av1 == a1) & (av2 == a4)) || (av1 == a4) & (av2 == a1)) && (design_mat[j, i] = -1.0)
        # leg 4
        (((av1 == a2) & (av2 == a3)) || (av1 == a3) & (av2 == a2)) && (design_mat[j, i] = -1.0)
    end
    return design_mat
end


function build_dmats(type::Symbol, closure, st)
    S = eltype(closure.T)
    dmat = Matrix{S}[]
    for i in 1:length(st)
        scanvis = st[i]
        inds = findall(==(scanvis.time), closure.T)
        if isnothing(inds) || isempty(inds)
            inow = length(dmat)
            # construct block diagonal matrices but be careful with scans that
            # cannot form closures: fold them into the previous scan with zeros
            if i > 1
                dmat[inow] = hcat(dmat[inow], zeros(S, size(dmat[inow], 1), length(scanvis.scan)))
            else
                push!(dmat, zeros(S, 1, length(scanvis.scan)))
            end
            continue
        end
        scancl = closure[inds]
        dmatscan = closure_designmat(type, scancl, scanvis)
        push!(dmat, dmatscan)
    end

    if iszero(dmat[1])
        dmat[2] = hcat(zeros(S, size(dmat[2], 1), size(dmat[1], 2)), dmat[2])
        dmat = dmat[2:end]
    end
    return dmat
end


function _minimal_closure(type, closures, st)
    S = eltype(closures.T)
    dmat = Matrix{S}[]
    for i in 1:length(st)
        scanvis = st[i]
        inds = findall(==(scanvis.time), closures.T)
        if length(inds) == 0
            inow = length(dmat)
            if i > 1
                dmat[inow] = hcat(dmat[inow], zeros(S, size(dmat[inow], 1), length(scanvis.scan)))
            else
                push!(dmat, zeros(S, 1, length(scanvis.scan)))
            end
            continue
        end

        scancl = closures[inds]

        # sort by closure noise so we form a nice minimal set
        snr = inv.(scancl.noise)
        ind_snr = sortperm(snr)
        scancl = scancl[ind_snr]

        design_mat = closure_designmat(type, scancl, scanvis)
        nmin = rank(design_mat)

        dmat_min = minimal_closure_scan(type, scancl, scanvis, nmin)
        push!(dmat, dmat_min)
    end

    if iszero(dmat[1])
        dmat[2] = hcat(zeros(S, size(dmat[2], 1), size(dmat[1], 2)), dmat[2])
        dmat = dmat[2:end]
    end
    return dmat
end


function minimal_closure_scan(type, closures, scanvis, nmin::Int)
    scancl = deepcopy(closures)
    keep = fill(true, length(closures))

    nmin0 = nmin

    count = 1
    keep = fill(true, length(closures))
    good = true
    while good
        closurekeep = closures[keep]
        design_mat = closure_designmat(type, closurekeep, scanvis)
        nmin = rank(design_mat)

        if (sum(keep) == nmin0) && (nmin == nmin0)
            good = false
        else
            if nmin == nmin0
                keep[count] = false
            else
                keep[count - 1] = true
                count -= 1
            end
        end
        count += 1
        count + 1 > length(scancl) && break
    end

    closurekeep = closures[keep]
    dmat = closure_designmat(type, closurekeep, scanvis)
    if length(closurekeep) != nmin
        @error "minimal set not found $(length(closurekeep)) $(nmin)"
        throw("No minimal set found at time $(scanvis.time)")
    end
    return dmat
end


"""
    build_closure_config(dvis::EHTObservationTable, closures::StructArray;
                         type::Symbol, count::String="min")

Build a [`ClosureConfig`](@ref) from a visibility observation table `dvis` and
a closure-spec table `closures` whose rows have at least the fields
`(T, F, noise, baseline)` (with `baseline` a `NTuple{3}` for closure phases or
`NTuple{4}` for log-closure amplitudes).

Arguments
- `type`: `:cphase` or `:lcamp`
- `count`: `"min"` (compute a minimal closure set per scan) or `"max"`
  (use the supplied closures as-is).
"""
function build_closure_config(
        dvis::EHTObservationTable, closures::StructArray;
        type::Symbol, count::String = "min"
    )
    st = timetable(dvis)
    dmat = if count == "min"
        _minimal_closure(type, closures, st)
    elseif count == "max"
        build_dmats(type, closures, st)
    else
        throw(ArgumentError("count=$count not valid; use \"min\" or \"max\""))
    end
    return ClosureConfig(arrayconfig(dvis), dmat, measurement(dvis), noise(dvis))
end
