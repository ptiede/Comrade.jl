# ComradeBase API



## Contents

```@contents
Pages = ["base_api.md"]
```

## Index

```@index
Pages = ["base_api.md"]
```

```@meta
CurrentModule = ComradeBase
```

## Model API

```@docs
ComradeBase.AbstractModel
ComradeBase.flux
ComradeBase.isprimitive
ComradeBase.visanalytic
ComradeBase.imanalytic
ComradeBase.ispolarized
ComradeBase.radialextent
ComradeBase.visibility_point
ComradeBase.intensity_point
ComradeBase.PrimitiveTrait
ComradeBase.IsPrimitive
ComradeBase.NotPrimitive
ComradeBase.DensityAnalytic
ComradeBase.IsAnalytic
ComradeBase.NotAnalytic
ComradeBase.IntensityMap
```

### Image Types
```@docs
ComradeBase.IntensityMap(::AbstractArray, ::AbstractDims)
ComradeBase.StokesIntensityMap
ComradeBase.imagepixels
ComradeBase.GriddedKeys
ComradeBase.dims
ComradeBase.named_dims
ComradeBase.axisdims
ComradeBase.stokes
ComradeBase.imagegrid
ComradeBase.fieldofview
ComradeBase.pixelsizes
ComradeBase.phasecenter
ComradeBase.centroid
ComradeBase.second_moment
```


## Polarization

```@docs
ComradeBase.AbstractPolarizedModel
ComradeBase.StokesParams
ComradeBase.ElectricFieldBasis
ComradeBase.RPol
ComradeBase.LPol
ComradeBase.XPol
ComradeBase.YPol
ComradeBase.PolBasis
ComradeBase.CirBasis
ComradeBase.LinBasis
ComradeBase.CoherencyMatrix
ComradeBase.evpa
ComradeBase.m̆
ComradeBase.linearpol
ComradeBase.SingleStokes
ComradeBase.innerprod
ComradeBase.basis_components
ComradeBase.basis_transform
```
