struct BaselineSiteLookup{V<:AbstractArray{<:Integer}}
    indices_1::V
    indices_2::V
end

function _construct_baselinemap(array::EHTArrayConfiguration, x::SiteArray)
    T = array[:T]
    F = array[:F]
    bl = array[:baseline]

    tcal = times(x)
    scal = sites(x)
    fcal = frequencies(x)
    tsf = StructArray((tcal, scal, fcal))
    ind1 = similar(Int, T)
    ind2 = similar(Int, T)

    for i in eachindex(T, F, bl, ind1, ind2)
        t = T[i]
        f = F[i]
        s1, s2 = bl[i]
        i1 = findfirst(x->(x[1]∈t)&&(x[2]∈f)&&(x[3]==s1), tsf)
        i2 = findfirst(x->(x[1]∈t)&&(x[2]∈f)&&(x[3]==s2), tsf)
        isnothing(i1) && throw(AssertionError("$t, $f, $((s1)) not found in SiteArray"))
        isnothing(i2) && throw(AssertionError("$t, $f, $((s2)) not found in SiteArray"))
        ind1[i] = i1
        ind2[i] = i2
    end
    BaselineSiteLookup(ind1, ind2)
end





struct RIMEModel{J, B<:BaselineSiteLookup, A<:AbstractArrayConfiguration, Ba} <: AbstractRIMEModel
    jones::J
    bsitemaps::B
    array::A
    refbasis::Ba
end

function RIMEModel(jones::JonesModel, prior::InstrumentPrior; refbasis=CirBasis())
    arr = array(prior)
    # First pre_allocate arrays if needed (i.e. feed rotations)
    jones2 = preallocate_jones(jones, arr, refbasis)
    # Construct the baseline site map for each prior
    x = rand(prior)
    bsitemaps = map(x->_construct_baselinemap(arr, x), x)
    # build the jones model
    return RIMEModel(jones2, bsitemaps, array, refbasis)
end



intout(vis::AbstractArray{<:StokesParams{T}}) where {T<:Real} = similar(vis, SMatrix{2,2, Complex{T}, 4})
intout(vis::AbstractArray{T}) where {T<:Real} = similar(vis, Complex{T})
intout(vis::AbstractArray{<:CoherencyMatrix{A,B,T}}) where {A,B,T<:Real} = similar(vis, SMatrix{2,2, Complex{T}, 4})

intout(vis::AbstractArray{<:StokesParams{T}}) where {T<:Complex} = similar(vis, SMatrix{2,2, T, 4})
intout(vis::AbstractArray{T}) where {T<:Complex} = similar(vis, T)
intout(vis::AbstractArray{<:CoherencyMatrix{A,B,T}}) where {A,B,T<:Complex} = similar(vis, SMatrix{2,2, T, 4})


function apply_instrument(vis, J::AbstractRIMEModel, x)
    vout = intout(vis)
    return apply_instrument!(vout, vis, J, x)
end

function apply_instrument!(vout, vis, J::RIMEModel, x)
    vout .= apply_jones.(vis, eachindex(vis), Ref(J), x)
    return nothing
end

get_indices(bsitemaps, index) = map(x->getindex(x, index), bsitemaps)
get_params(x::NamedTuple{N}, indices::NamedTuple{N}) where {N} = NamedTuple{N}(map((xx, ii)->getindex(xx, ii), x, indices))

function apply_jones(v, index::Int, J::RIMEModel, x)
    indices1 = get_indices(J.bsitemaps.indices_1, index)
    indices2 = get_indices(J.bsitemaps.indices_2, index)
    params1 = get_params(x, indices1)
    params2 = get_params(x, indices2)
    j1 = jonesmatrix(J.jones, params1, index, Val(1))
    j2 = jonesmatrix(J.jones, params2, index, Val(2))
    return _apply_jones(v, j1, j2, J.refbasis)
end

_apply_jones(v::Number, j1, j2, ::B) where {B} = j1*v*conj(j2)
_apply_jones(v::CoherencyMatrix, j1, j2, ::B) where {B} = j1*CoherencyMatrix{B,B}(v)*j2'
_apply_jones(v::StokesParams, j1, j2, ::B) where {B} = j1*CoherencyMatrix{B,B}(v)*j2'



struct IdealRIMEModel <: AbstractRIMEModel end
apply_instrument(vis, ::IdealRIMEModel, x) = vis
