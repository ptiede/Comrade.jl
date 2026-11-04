# Traits that decide the domain of the data. We include this to prevent additional computation if
# e.g., on visibility data is considered.
abstract type ComradeDataType end
struct NoData <: ComradeDataType end
struct DualData <: ComradeDataType end
struct VisData <: ComradeDataType end
struct ImgData <: ComradeDataType end

datatype(::Type, ::Type{<:Nothing}) = VisData()
datatype(::Type{<:Nothing}, ::Type) = ImgData()
datatype(::Type{<:Nothing}, ::Type{<:Nothing}) = throw(ArgumentError("No data in uv plane or image plane is provided"))
datatype(::Type, ::Type) = DualData()

"""
    datadomain(::Type)

The domain a data product of a [`VLBIPosterior`](@ref) constrains: `VisDomain()` for products
compared with the model visibilities through the instrument model ([`VisTable`](@ref),
[`ClosureTable`](@ref)), `ImageDomain()` for products compared with the sky model's image
([`ImgNormalData`](@ref)).

A new image-domain data type `T` defines `Comrade.datadomain(::Type{<:T}) = ImageDomain()` and
`Comrade.makelikelihood(d::T)`, which returns an object `ℓ` whose `logdensityof(ℓ, img)` is the
log-likelihood of the model image `img`.
"""
datadomain(T::Type) = throw(
    ArgumentError(
        "$T is not a data product of a posterior: define `Comrade.datadomain(::Type{<:$(nameof(T))})` " *
            "and `Comrade.makelikelihood`"
    )
)
datadomain(d) = datadomain(typeof(d))

"""
    VisDomain()

The [`datadomain`](@ref) of data products compared with the model visibilities.
"""
struct VisDomain end

"""
    ImageDomain()

The [`datadomain`](@ref) of data products compared with the model image.
"""
struct ImageDomain end
