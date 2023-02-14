# Model Interface


## Primitive Geometric Models

`Comrade` aims to be more modular and extensible than previous VLBI modeling packages. Namely, simple models are composed to construct complicated source morphologies instead of making many different models. This is accomplished with a type and trait-based hierarchy.

Additionally, [`ComradeBase`](https://github.com/ptiede/ComradeBase.jl) is a low-dependency version of this package that defines this type and trait hierarchy that users can more easily incorporate into their packages.

To see how this works, we will go through a simplified implementation of the Gaussian model in `Comrade`. The Gaussian model is a simple, compact emission structure model that can constrain the typical characteristic size
of an image feature from VLBI data. To construct a Gaussian model, we will first define a struct:

```julia
struct MyGaussian <: Comrade.AbstractModel end
```

Notice that we don't provide any more information about the model, e.g., *size, shape, flux*, etc. This is because we will use `Comrade`'s extensive set of modifiers to change the structure of the model.
Now a `Gaussian` is the simplest model structure we can consider. We can consider this Gaussian to be a **primitive** model. That means a Gaussian is not a combination or modification of
an existing model. To tell `Comrade` that this is the case, we define the following method:

```julia
# Tell Comrade Gaussian is a primitive model
ComradeBase.isprimitive(::Type{<:MyGaussian}) = IsPrimitive()
```

Note that if the Gaussian wasn't a primitive model, we could've used `NotPrimitive()` instead.
Now a Gaussian has an analytic expression in the image and Fourier domain. We can tell `Comrade` this by setting:

```julia
# Fourier and image domains are analytic
Comrade.visanalytic(::Type{<:MyGaussian}) = IsAnalytic()
Comrade.imanalytic(::Type{<:MyGaussian}) = IsAnalytic()
```

Finally, we can specify if the model is intrinsically polarized by using the `IsPolarized` and `NotPolarized()` trait
```julia
Comrade.ispolarized(::Type{<:MyGaussian}) = NotPolarized()
```
!!! note
    The actual implementation defines the Gaussian to be a subtype of `Comrade.GeometricModel`, which automatically defines these methods. However, for models that aren't a subtype of `GeometricModel`, we assume the image domain `IsAnalytic()` and the Fourier domain is `NotAnalytic()`.

Since both the image and visibility domain representation of the Gaussian are analytic, we need to define an `intensity_point` and `visibility_point` method. For a Gaussian, these are given by

```julia
function ComradeBase.intensity_point(::MyGaussian, p)
    (;X, Y) = p
    return exp(-(X^2+Y^2)/2)/2π
end

function ComradeBase.visibility_point(::MyGaussian, u, v, time, freq) where {T}
    return exp(-2π^2*(u^2 + v^2)) + 0im
end
```

Additionally, most models in `Comrade` has two additional functions one can implement if possible:

1. `flux(m::MyGaussian)`: This defines the flux of a model. If this isn't defined, the model won't have a flux until an image is created. For a Gaussian, the definition is `flux(::MyGaussian) = 1.0`.
2. `radialextent(::MyGaussian)`: This defines the model's default radial extent. For a Gaussian, we will consider the radial extent to be $5σ$, so `radialextent(::MyGaussian) = 5.0`.

This completely defines the model interface for `Comrade`. With this, you can call the usual user API to evaluate, fit, and plot the model. Additionally, we can start talking about
adding multiple Gaussians and modifying them. For instance, suppose you want an elliptical Gaussian with a flux of 2 Jy. This can be created by `Comrade` as follows:

```julia
gauss = MyGaussian()
ellgauss = 2.0*rotated(stretched(gauss, 1.0, 0.5), π/4)
fig = plot(gauss, layout=(1,2), size=(800,300))
plot!(fig[2], ellgauss, size=(800,350))
```
![Image](assets/gauss.png)


```julia
u = rand(100)*0.5; v=rand(100)*0.5
vg  = visibilities(gauss, u, v)
veg = visibilities(ellgauss, u, v)

scatter(hypot.(u, v), abs.(vg), label="Gaussian")
scatter!(hypot.(u, v), abs.(veg), label="Elliptical Gaussian")
```
![Image](assets/vis.png)

## Models without an Analytic Fourier Transform

Now suppose your model does not have an analytic Fourier transform. In this case, the procedure is very similar to the above, except you define `visanalytic(::Type{<:MyModel}) = NotAnalytic()`.
However, everything else is the same. To compute visibilities, you just then create a `ModelImage` type using the [`modelimage`](@ref) function. To see how this see [Modeling with non-analytic Fourier transforms](@ref).
