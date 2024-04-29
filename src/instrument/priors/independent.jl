export IIDSitePrior

abstract type AbstractSitePrior end

segmentation(d::AbstractSitePrior) = getfield(d, :seg)

struct IIDSitePrior{S<:Segmentation, D} <: AbstractSitePrior
    seg::S
    dist::D
end
