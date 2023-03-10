export stretched, shifted, rotated, renormed, modify, Stretch, Renormalize, Shift, Rotate


"""
    $(TYPEDEF)

General type for a model modifier. These transform any model
using simple Fourier transform properties. To modify a model
you can use the [`ModifiedModel`](@ref) constructor or the [`modify`](@ref)
function.

```julia-repl
julia> visanalytic(stretched(Disk(), 2.0, 2.0)) == visanalytic(Disk())
true
```



To implement a model transform you need to specify the following methods:
- [`transform_uv`](@ref)
- [`transform_image`](@ref)
- [`scale_uv`](@ref)
- [`scale_image`](@ref)
- [`radialextent`](@ref)
See thee docstrings of those methods for guidance on implementation details.

Additionally these methods assume the modifiers are of the form

I(x,y) -> fᵢ(x,y)I(gᵢ(x,y))
V(u,v) -> fᵥ(u,v)V(gᵥ(u,v))

where `g` are the transform_image/uv functions and `f` are the scale_image/uv
function.

"""
abstract type ModelModifier end




"""
    scale_image(model::AbstractModifier, x, y)

Returns a number of how to to scale the image intensity at `x` `y` for an modified `model`
"""
function scale_image end


"""
    transform_image(model::AbstractModifier, x, y)

Returns a transformed `x` and `y` according to the `model` modifier
"""
function transform_image end

"""
    scale_image(model::AbstractModifier, u, u)

Returns a number on how to scale the image visibility at `u` `v` for an modified `model`
"""
function scale_uv end

"""
    transform_uv(model::AbstractModifier, u, u)

Returns a transformed `u` and `v` according to the `model` modifier
"""
function transform_uv end

unitscale(T, ::Type{M}) where {M} = unitscale(T, ispolarized(M))
unitscale(T, ::NotPolarized) = one(T)
unitscale(T, ::IsPolarized) = I

struct ModifiedModel{M<:AbstractModel, T<:Tuple}<: AbstractModel
    model::M
    transform::T
end

function Base.show(io::IO, mi::ModifiedModel)
    #io = IOContext(io, :compact=>true)
    #s = summary(mi)
    println(io, "ModifiedModel")
    println(io, "  base model: ", summary(mi.model))
    println(io, "  Modifiers:")
    for i in eachindex(mi.transform)
        println(io, "    $i. ", summary(mi.transform[i]))
    end
 end


"""
    unmodified(model::ModifiedModel)

Returns the un-modified model

### Example
```julia-repl
julia> m = stretched(rotated(Gaussian(), π/4), 2.0, 1.0)
julia> umodified(m) == Gaussian()
true
```
"""
unmodified(model::ModifiedModel) = model.model
unmodified(model::AbstractModel) = model

"""
    basemodel(model::ModifiedModel)

Returns the ModifiedModel with the last transformation stripped.

# Example
```julia-repl
julia> basemodel(stretched(Disk(), 1.0, 2.0)) == Disk()
true
```
"""
basemodel(model::ModifiedModel) = ModifiedModel(model.model, Base.front(model.transform))
basemodel(model::ModifiedModel{M, Tuple{}}) where {M} = model


flux(m::ModifiedModel) = flux(m.model)*mapreduce(flux, Base.:*, m.transform)
flux(t::ModelModifier) = 1


radialextent(m::ModifiedModel) = radialextent_modified(radialextent(m.model), m.transform)

@inline visanalytic(::Type{ModifiedModel{M, T}}) where {M,T} = visanalytic(M)
@inline imanalytic(::Type{ModifiedModel{M, T}}) where {M,T} = imanalytic(M)
@inline ispolarized(::Type{ModifiedModel{M, T}}) where {M,T} = ispolarized(M)



@inline function ModifiedModel(m::AbstractModel, t::ModelModifier)
    return ModifiedModel(m, (t,))
end

@inline function ModifiedModel(m::ModifiedModel, t::ModelModifier)
    model = m.model
    t0 = m.transform
    return ModifiedModel(model, (t0..., t))
end


function modify_uv(model, t::Tuple, u, v, scale)
    ut, vt = transform_uv(model, last(t), u, v)
    scalet = scale_uv(model, last(t), u, v)
    return modify_uv(model, Base.front(t), ut, vt, scalet*scale)
end
modify_uv(model, ::Tuple{}, u, v, scale) = (u, v), scale

function modify_image(model, t::Tuple, x, y, scale)
    xt, yt = transform_image(model, last(t), x, y)
    scalet = scale_image(model, last(t), x, y)
    return modify_image(model, Base.front(t), xt, yt, scalet*scale)
end
modify_image(model, ::Tuple{}, x, y, scale) = x, y, scale


# @inline function transform_uv(model, t::Tuple, u, v)
#     ut, vt = transform_uv(model, last(t), u, v)
#     return transform_uv(model, Base.front(t), ut, vt)
# end
# @inline transform_uv(model, ::Tuple{}, u, v) = u, v

# @inline function scale_uv(model, t::Tuple, u, v)
#     scale = scale_uv(model, last(t), u, v)
#     return scale*scale_uv(model, Base.front(t), u, v)
# end
# @inline scale_uv(::M, t::Tuple{}, u, v) where {M} = unitscale(eltype(u), M)

# @inline function transform_image(model, t::Tuple, x, y)
#     xt, yt = transform_image(model, last(t), x, y)
#     return transform_image(model, Base.front(t), xt, yt)
# end
# @inline transform_image(model, ::Tuple{}, x, y) = x, y

# @inline function scale_image(model, t::Tuple, x, y)
#     scale = scale_image(model, last(t), x, y)
#     return scale*scale_image(model, Base.front(t), x, y)
# end
# @inline scale_image(::M, t::Tuple{}, x, y) where {M} = unitscale(eltype(x), M)

@inline radialextent_modified(r::Real, t::Tuple) = radialextent_modified(radialextent_modified(r, last(t)), Base.front(t))
@inline radialextent_modified(r::Real, ::Tuple{}) = r


"""
    modify(m::AbstractModel, transforms...)

Modify a given `model` using the set of `transforms`. This is the most general
function that allows you to apply a sequence of model transformation for example

```julia-repl
modify(Gaussian(), Stretch(2.0, 1.0), Rotate(π/4), Shift(1.0, 2.0), Renorm(2.0))
```
will create a asymmetric Gaussian with position angle `π/4` shifted to the position
(1.0, 2.0) with a flux of 2 Jy. This is similar to Flux's chain.
"""
function modify(m::AbstractModel, transforms...)
    return ModifiedModel(m, transforms)
end




# @inline function apply_uv_transform(m::AbstractModifier, t::TransformState)
#     ut, vt = transform_uv(m, t.u, t.v)
#     scale = t.scale*scale_uv(m, t.u, t.v)
#     return apply_uv_transform(basemodel(m), TransformState(ut, vt, scale))
# end

# @inline function apply_uv_transform(::AbstractModel, t::TransformState)
#     return t
# end


# @inline function apply_uv_transform(m::AbstractModifier, u, v, scale)
#     ut, vt = transform_uv(m, u, v)
#     scale = scale*scale_uv(m, u, v)
#     return apply_uv_transform(basemodel(m), ut, vt, scale)
# end

# @inline function apply_uv_transform(::AbstractModel, u, v, scale)
#     return (u, v), scale
# end

# @inline function _visibilities(m::AbstractModifier, u, v, time, freq)
#     uv, scale = apply_uv_transform(m, u, v)
#     ut = first.(uv)
#     vt = last.(uv)
#     scale.*_visibilities(unmodified(m), ut, vt, time, freq)
# end


# # function visibilities(m, p::NamedTuple)
# #     m = Base.Fix1(m∘NamedTuple{keys(p)})
# #     return visibilities(m, NamedTuple{keys(p)}(p))
# # end

# function apply_uv_transform(m::AbstractModifier, u::AbstractVector, v::AbstractVector)
#     res = apply_uv_transform.(Ref(m), u, v, 1.0)
#     return first.(res), last.(res)
# end



# function apply_uv_transform(m::AbstractModifier, u::AbstractVector, v::AbstractVector)
#     res = apply_uv_transform.(Ref(m), u, v, 1.0)
#     return getindex.(res,1), getindex.(res,2), getindex.(res,3)
#     res = apply_uv_transform.(Ref(m), u, v, 1.0)
#     return getindex.(res,1), getindex.(res,2), getindex.(res,3)
# end
# @inline function _visibilities(m::M, p) {M<:AbstractModifier}

# @inline function _visibilities(m::M, p) {M<:AbstractModifier}

#     return _visibilities(visanalytic(M), m, u, v, args...)
# end

# @inline function _visibilities(m::AbstractModifier{M}, p) where {M}
#     return _visibilities(ispolarized(M), m, p)
# end



# @inline function _visibilities(m::AbstractModifier, p)
#     (;U, V) = p
#     st = StructArray{TransformState{eltype(U), Complex{eltype(U)}}}(u=U, v=V, scale=fill(one(Complex{eltype(U)}), length(U)))
#     auv = Base.Fix1(apply_uv_transform, m)
#     mst = map(auv, st)
#     mst.scale.*visibilities(unmodified(m), (U=mst.u, V=mst.v))
# end

function update_uv(p::NamedTuple, uv)
    p1 = @set p.U = uv.U
    p2 = @set p1.V = uv.V
    return p2
end

function update_xy(p::NamedTuple, xy)
    p1 = @set p.X = xy.X
    p2 = @set p1.Y = xy.Y
    return p2
end


# @inline function _visibilities(::IsPolarized, m::AbstractModifier, p)
#     (;U, V) = p

#     S = eltype(U)
#     unit = StokesParams(complex(one(S)), complex(one(S)), complex(one(S)),complex(one(S)))
#     st = StructArray{TransformState{eltype(U), typeof(unit)}}(u=U, v=V, scale=Fill(unit, length(U)))
#     mst = apply_uv_transform.(Ref(m), st)

#     pup = update_uv(p, (U=mst.u, V=mst.v))
#     mst.scale.*visibilities(unmodified(m), pup)
# end

# @inline function _visibilities(m::AbstractModifier, p)
#     (;U, V) = p
#     st = StructArray{TransformState{eltype(U), Complex{eltype(U)}}}(u=U, v=V, scale=fill(one(Complex{eltype(U)}), length(U)))
#     mst = apply_uv_transform.(Ref(m), st)
#     pup = update_uv(p, (U=mst.u, V=mst.v))
#     mst.scale.*visibilities(unmodified(m), pup)
# end



# I need some special pass-throughs for the non-analytic FFT transform
# since you evaluate the visibilities pointwise
function modelimage(::NotAnalytic,
    model::ModifiedModel,
    image::IntensityMap, alg::FFTAlg, pulse = DeltaPulse(), thread::StaticBool=False())

    @set model.model = modelimage(NotAnalytic(), model.model, image, alg, pulse, thread)
end

# I need some special pass-throughs for the non-analytic NUFT transform
# since you evaluate the visibilities as a vector
function modelimage(::NotAnalytic,
    model::ModifiedModel,
    image::IntensityMap, alg::NUFT,
    pulse = DeltaPulse(), thread::StaticBool=False())
    _modelimage(model, image, alg, pulse, thread)
end

@inline function modelimage(::NotAnalytic, model, cache::FFTCache, thread::StaticBool)
    img = cache.img
    intensitymap!(img, model, thread)
    newcache = update_cache(cache, img, cache.pulse)
    return ModelImage(model, img, newcache)
end



@inline function visibility_point(m::ModifiedModel, u, v, time, freq)
    mbase = m.model
    transform = m.transform
    (ut, vt), scale = modify_uv(mbase, transform, u, v, unitscale(Complex{eltype(u)}, typeof(mbase)))
    scale*visibility_point(mbase, ut, vt, time, freq)
end

@inline function intensity_point(m::ModifiedModel, p)
    mbase = m.model
    transform = m.transform
    xt, yt, scale = modify_image(mbase, transform, p.X, p.Y, unitscale(eltype(p.X), typeof(mbase)))
    scale*intensity_point(mbase, update_xy(p, (X=xt, Y=yt)))
end

function modify_uv(model, transform::Tuple, u::AbstractVector, v::AbstractVector, scale)
    uvscale = modify_uv.(Ref(model), Ref(transform), u, v, scale)
    return first.(uvscale), last.(uvscale)
end


function visibilities_analytic(m::ModifiedModel, u, v, time, freq)
    T = typeof(visibility_point(m, first(u), first(v), first(time), first(freq)))
    vis = similar(u, T)
    visibilities_analytic!(vis, m, u, v, time, freq)
    return vis
end

using NamedTupleTools
function ChainRulesCore.rrule(::typeof(visibilities_analytic), m::ModifiedModel, u, v, time, freq)
    vis = visibilities_analytic(m, u, v, time, freq)
    function _visibilities_analytic_mm_pullback(Δ)
        du = zero(u)
        dv = zero(v)
        df = zero(freq)
        dt = zero(time)

        dvis = zero(vis)
        dvis .= unthunk(Δ)
        rvis = zero(vis)
        d = autodiff(Reverse, visibilities_analytic!, Const, Duplicated(rvis, dvis), Active(m), Duplicated(u, du), Duplicated(v, dv), Duplicated(time, dt), Duplicated(freq, df))
        dm = d[1]
        tm = __extract_tangent(dm)
        return NoTangent(), tm, du, dv, df, dt
    end
    return vis, _visibilities_analytic_mm_pullback
end


function __extract_tangent(dm::ModifiedModel)
    tm = __extract_tangent(dm.model)
    dtm = dm.transform
    ttm = map(x->Tangent{typeof(x)}(;ntfromstruct(x)...), dtm)
    tm = Tangent{typeof(dm)}(model=tm, transform=ttm)
end


function __extract_tangent(dm::GeometricModel)
    ntm = ntfromstruct(dm)
    if ntm isa NamedTuple{(), Tuple{}}
        tbm = ZeroTangent()
    else
        tbm = Tangent{typeof(dm)}(;ntm...)
    end
    return tbm
end







"""
    $(TYPEDEF)

Shifts the model by `Δx` units in the x-direction and `Δy` units
in the y-direction.

An end user should not call this directly but instead
the [`shifted`](@ref) function instead.
"""
struct Shift{T} <: ModelModifier
    Δx::T
    Δy::T
end

#function ShiftedModel(model::AbstractModel, Δx::Number, Δy::Number)
#    T =
#    return ShiftedModel{(model, promote(Δx, Δy)...)
#end

"""
    $(SIGNATURES)

Shifts the model `m` in the image domain by an amount `Δx,Δy`
in the x and y directions respectively.
"""
shifted(model, Δx, Δy) = ModifiedModel(model, Shift(Δx, Δy))
# This is a simple overload to simplify the type system
@inline radialextent_modified(r::Real, t::Shift) = r + max(abs(t.Δx), abs(t.Δy))

@inline transform_image(model, transform::Shift, x, y) = (x-transform.Δx, y-transform.Δy)
@inline transform_uv(model, ::Shift, u, v) = (u, v)

@inline scale_image(::M, transform::Shift, x, y) where {M} = unitscale(typeof(transform.Δx), M)
# Curently we use exp here because enzyme has an issue with cispi that will be fixed soon.
@inline scale_uv(::M, transform::Shift, u, v) where {M}    = exp(2im*π*(u*transform.Δx + v*transform.Δy))*unitscale(typeof(transform.Δx), M)



"""
    $(TYPEDEF)

Renormalizes the flux of the model to the new value `scale*flux(model)`.
We have also overloaded the Base.:* operator as syntactic sugar
although I may get rid of this.

An end user should not call this directly but instead
the [`renormed`](@ref) function or Base.:* instead.

# Example
```julia-repl
julia> renormed(Gaussian(), 2.0) == 2.0*Gaussian()
true
```
"""
struct Renormalize{T} <: ModelModifier
    scale::T
end



"""
    $(SIGNATURES)

Renormalizes the model `m` to have total flux `f*flux(m)`.
This can also be done directly by calling `Base.:*` i.e.,

```julia-repl
julia> renormed(m, f) == f*M
true
```
"""
renormed(model::M, f) where {M<:AbstractModel} = ModifiedModel(model, Renormalize(f))
Base.:*(model::AbstractModel, f::Number) = renormed(model, f)
Base.:*(f::Number, model::AbstractModel) = renormed(model, f)
Base.:/(f::Number, model::AbstractModel) = renormed(model, inv(f))
Base.:/(model::AbstractModel, f::Number) = renormed(model, inv(f))
# Dispatch on RenormalizedModel so that I just make a new RenormalizedModel with a different f
# This will make it easier on the compiler.
# Base.:*(model::ModifiedModel, f::Number) = renormed(model.model, model.scale*f)
# Overload the unary negation operator to be the same model with negative flux
Base.:-(model::AbstractModel) = renormed(model, -1.0)
flux(t::Renormalize) = t.scale

@inline transform_image(m, ::Renormalize, x, y) = (x, y)
@inline transform_uv(m, ::Renormalize, u, v) = (u, v)

@inline scale_image(::M, transform::Renormalize, x, y) where {M} = transform.scale*unitscale(typeof(transform.scale), M)
@inline scale_uv(::M, transform::Renormalize, u, v) where {M}    = transform.scale*unitscale(typeof(transform.scale), M)

@inline radialextent_modified(r::Real, ::Renormalize) = r


"""
    $(TYPEDEF)

Stretched the model in the x and y directions, i.e. the new intensity is
    Iₛ(x,y) = 1/(αβ) I(x/α, y/β),
where were renormalize the intensity to preserve the models flux.

An end user should not call this directly but instead
the [`stretched`](@ref) function instead.
"""
struct Stretch{T} <: ModelModifier
    α::T
    β::T
end


"""
    $(SIGNATURES)

Stretches the model `m` according to the formula
    Iₛ(x,y) = 1/(αβ) I(x/α, y/β),
where were renormalize the intensity to preserve the models flux.
"""
stretched(model, α, β) = ModifiedModel(model, Stretch(α, β))

@inline transform_image(m, transform::Stretch, x, y) = (x/transform.α, y/transform.β)
@inline transform_uv(m, transform::Stretch, u, v) = (u*transform.α, v*transform.β)

@inline scale_image(::M, transform::Stretch{T}, x, y) where {M,T} = inv(transform.α*transform.β)*unitscale(T, M)
@inline scale_uv(::M, ::Stretch{T}, u, v) where {M,T} = unitscale(T, M)

@inline radialextent_modified(r::Real, t::Stretch) = r*max(t.α, t.β)


"""
    $(TYPEDEF)

Type for the rotated model. This is more fine grained constrol of
rotated model.

An end user should not call this directly but instead
the [`rotated`](@ref) function instead.
"""
struct Rotate{T} <: ModelModifier
    s::T
    c::T
end
function Rotate(ξ::F) where {F}
    s,c = sincos(ξ)
    return Rotate(s, c)
end


"""
    $(SIGNATURES)

Rotates the model by an amount `ξ` in radians in the clockwise direction.
"""
rotated(model, ξ) = ModifiedModel(model, Rotate(ξ))

"""
    $(SIGNATURES)

Returns the rotation angle of the rotated `model`
"""
posangle(model::Rotate) = atan(model.s, model.c)

@inline function transform_image(m, transform::Rotate, x, y)
    s,c = transform.s, transform.c
    return c*x - s*y, s*x + c*y
end

@inline function transform_uv(m, transform::Rotate, u, v)
    s,c = transform.s, transform.c
    return c*u - s*v, + s*u + c*v
end


@inline function scale_image(::M, transform::Rotate{T}, x, y) where {M,T}
    scale_image(ispolarized(M), transform, x, y)
end


@inline scale_image(::NotPolarized, model::Rotate{T}, x, y) where {T} = one(T)

@inline function spinor2_rotate(c, s)
    u = oneunit(c)
    z = zero(s)
    c2 = c^2 - s^2
    s2 = 2*c*s
    return SMatrix{4,4}(u,  z,  z,  z,
                        z,  c2, s2, z,
                        z, -s2, c2, z,
                        z,   z,  z, u)

end

@inline function scale_image(::IsPolarized, model::Rotate, x, y)
    return spinor2_rotate(model.c, model.s)
end

@inline function scale_uv(::M, model::Rotate, u, v) where {M}
    scale_uv(ispolarized(M), model, u, v)
end


@inline scale_uv(::NotPolarized, model::Rotate{T}, u, v) where {T} = one(T)

@inline function scale_uv(::IsPolarized, model::Rotate, x, y)
    return spinor2_rotate(model.c, model.s)
end
@inline radialextent_modified(r::Real, ::Rotate) = r
