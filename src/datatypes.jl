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
