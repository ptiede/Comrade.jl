export modelimage

abstract type AbstractModelImage{M} <: AbstractModel end

struct ModelImage{M,I,C} <: AbstractModelImage{M}
    model::M
    image::I
    cache::C
end
@inline visanalytic(::Type{<:ModelImage{M}}) where {M} = visanalytic(M)
@inline imanalytic(::Type{<:ModelImage{M}}) where {M} = imanalytic(M)
@inline isprimitive(::Type{<:ModelImage{M}}) where {M} = isprimitive(M)

function Base.show(io::IO, mi::ModelImage)
   #io = IOContext(io, :compact=>true)
   #s = summary(mi)
   ci = first(split(summary(mi.cache), "{"))
   println(io, "ModelImage")
   println(io, "\tmodel: ", summary(mi.model))
   println(io, "\timage: ", summary(mi.image))
   print(io, "\tcache: ", ci)
end

intensitymap!(mimg::ModelImage) = intensitymap!(mimg.image, mimg.model)

flux(mimg::ModelImage) = flux(intensitymap!(mimg.image, mimg.model))


@inline visibility_point(m::AbstractModelImage, u, v) = visibility_point(m.model, u, v)

@inline intensity_point(m::ModelImage, x, y) = intensity_point(m.model, x, y)

"""
    $(SIGNATURES)
Construct a `ModelImage` from a `model`, `image` and the optionally
specified visibility algorithm `alg`

# Notes
For analytic models this is a no-op and just return the model.
For non-analytic models this wraps the model in a object with an image
and precomputes the fourier transform using `alg`.
"""
@inline function modelimage(model::M, image; alg=FFT()) where {M}
    return modelimage(visanalytic(M), model, image; alg)
end

@inline function modelimage(::IsAnalytic, model, image; alg=FFT())
    return model
end

@inline function modelimage(::NotAnalytic, model, image; alg=FFT())
    cache = create_cache(alg, model, image)
    return ModelImage(model, image, cache)
end

"""
    modelimage(m)
Construct a `ModelImage` where just the model `m` is specified

# Notes
Currently this is only defined for analytic models. In the future
this will *guess* a reasonable image to use.
"""
function modelimage(m::M) where {M}
    if visanalytic(M) == IsAnalytic()
        return ModelImage(m, Matrix{Float64}(undef, 1, 1), NoCache())
    else
        throw(ArgumentError("$m is not an analytic model a image must be specified"))
    end
end
