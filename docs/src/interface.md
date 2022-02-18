# Model Interface

`Comrade` aims to be more modular and extensible than previous VLBI modeling packages. Namely, instead of making many different models, simple models are composed to construct complicated source morphologies. This is accomplished with a type, and trait based hierarchy.

To see how this works we will implement a  model

```julia
struct Gamma{T} <: Comrade.AbstractModel
    α::T
end
```

This will a compact source that has more extent than a usual Gaussian blob. This is a **primtive** model in that it can't easily be constructed from other models. Additionally, this model has analytic expressions in the image and Fourier domain. To tell `Comrade` that `Gamma` has these properties we need to specify a couple of traits

```julia
# Tell Comrade Gamma is a primitive model
Comrade.isprimitive(::Type{<:Gamma}) = IsPrimitive()

# Fourier and image domain are analytic
Comrade.visanalytic(::Type{<:Gamma}) = IsAnalytic()
Comrade.imanalytic(::Type{<:Gamma}) = IsAnalytic()
```

Now since both the image and visibilities are analytic we need to specify how to calculate them. First we will specify the image domain function:

```julia
function Comrade.intensity_point(m::Gamma, x, y)
    r = hypot(x,y) + eps()
    return inv(2π*gamma(α))*r^(α-1)*exp(-r)
end
```

Similarly we specify the visibility function as follows:

```julia
function Comrade.visibility_point(m, u, v)
    return (1- )
end
```

