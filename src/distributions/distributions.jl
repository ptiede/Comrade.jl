using Distributions:  ContinuousUnivariateDistribution,
                      @distr_support,
                      params, partype,
                      mean, var,
                      pdf, logpdf, cdf, quantile,
                      sampler

struct Rice{T<:Real} <: ContinuousUnivariateDistribution
    ν::T
    σ::T
    Rice{T}(μ::T, σ::T) where {T} = new{T}(μ,σ)
end

function Rice(ν::T, σ::T; check_args=true) where {T<:Real}
    check_args && @check_args(Rice, ν≥zero(ν), σ≥zero(σ))
    return Rice{T}(μ, σ)
end

Rice(ν::Real, σ::Real) = Rice(promote(ν, σ)...)
Rice(ν::Integer, σ::Integer) = Rice(float(ν), float(σ))

@distr_support Rice 0.0 Inf

params(d::Rice) = (d.ν, d.σ,)
partype(::Rice{T}) where {T<:Real} = T

#### Statistics
@memoize L12(x) = exp(x/2)*( (1-x)*besseli(0,-x/2) - x*besseli(1,-x/2) )
mean(d::Rice) = d.σ*Distributions.sqrt2π/2*L12(-(d.ν/(2d.σ))^2)
var(d::Rice) = 2*d.σ^2 + d.ν^2 - π*d.σ^2/2*L12(-(ν/(2d.σ))^2)^2

pdf(d::Rice{T}, x) where {T} =
        x ≥ zero(T) ? x/d.σ^2*exp(-(x^2+d.ν^2)/(2d.σ^2))*besseli(0, x*d.ν/d.σ^2) : zero(T)
logpdf(d::Rice{T}, x) where {T} =
        x ≥ zero(T) ? log(x/d.σ^2*besseli(0, x*d.ν/d.σ^2)) - (x^2 + ν^2)/(2*d.σ^2) : -T(Inf)
