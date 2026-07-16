abstract type ReferencingScheme end

export NoReference, SingleReference, SEFDReference

struct NoReference <: ReferencingScheme end

"""
    SingleReference(site::Symbol, val)

Selects a single reference site for all scans. The value of the site is set to `val`.
"""
struct SingleReference{T} <: ReferencingScheme
    site::Symbol
    value::T
end


struct SEFDReference{T} <: ReferencingScheme
    value::T
    offset::Int
end

"""
    SEFDReference(val::Number, sefd_index = 1)

Selects the reference site based on the SEFD of each telescope, where the smallest SEFD
is preferentially selected. The reference gain is set to `val` and the user can select to
use the `n` lowest SEFD site by passing `sefd_index = n`.

## Notes
This is done on a per-scan basis so if a site is missing from a scan the next highest SEFD
site will be used.
"""
SEFDReference(val::Number) = SEFDReference(val, 0)

reference_indices(::AbstractArrayConfiguration, st::SiteLookup, ::NoReference) = [], nothing
function reference_indices(::AbstractArrayConfiguration, st::SiteLookup, p::SingleReference)
    inds = findall(==(p.site), sites(st))
    return inds, fill(p.value, length(inds))
end

function reference_indices(array::AbstractArrayConfiguration, st::SiteLookup, r::SEFDReference)
    tarr = array.tarr
    sefd = NamedTuple{Tuple(tarr.sites)}(Tuple(tarr.SEFD1 .+ tarr.SEFD2))
    # A reference group is a (time entry, frequency entry) cell of the lookup's axes.
    # Time membership is by code; frequency membership is by channel *overlap* with the
    # cell's frequency entry, so band-spanning parameters (`freqseg = FullBand()`)
    # compete in every per-channel cell they cover instead of forming their own group.
    # For uniform frequency segmentations overlap reduces to code equality — the classic
    # per-(scan, channel) fixing. A parameter that wins several cells is fixed once.
    best = Dict{Tuple{Int, Int}, Int}()
    order = Tuple{Int, Int}[] # first-seen (storage) order, for deterministic output
    for i in eachindex(st.tcode)
        s = sefd[st.saxis[st.scode[i]]]
        qlo = st.flo[st.fcode[i]]
        qhi = st.fhi[st.fcode[i]]
        for fe in eachindex(st.faxis)
            (st.flo[fe] < qhi && st.fhi[fe] > qlo) || continue
            key = (Int(st.tcode[i]), fe)
            b = get(best, key, 0)
            if b == 0
                best[key] = i
                push!(order, key)
            elseif s < sefd[st.saxis[st.scode[b]]]
                best[key] = i
            end
        end
    end
    fixedinds = unique!([best[k] for k in order])
    return fixedinds, fill(r.value, length(fixedinds))
end
