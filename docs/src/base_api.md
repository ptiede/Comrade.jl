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
ComradeBase.flux
ComradeBase.visibility
ComradeBase.visibilitymap
ComradeBase.visibilitymap!
ComradeBase.intensitymap
ComradeBase.intensitymap!
ComradeBase.IntensityMap
ComradeBase.amplitude(::Any, ::Any)
ComradeBase.amplitudes
ComradeBase.bispectrum
ComradeBase.bispectra
ComradeBase.closure_phase
ComradeBase.closure_phases
ComradeBase.logclosure_amplitude
ComradeBase.logclosure_amplitudes
```

### Model Interface
```@docs
ComradeBase.AbstractModel
ComradeBase.isprimitive
ComradeBase.visanalytic
ComradeBase.imanalytic
ComradeBase.ispolarized
ComradeBase.radialextent
ComradeBase.PrimitiveTrait
ComradeBase.IsPrimitive
ComradeBase.NotPrimitive
ComradeBase.DensityAnalytic
ComradeBase.IsAnalytic
ComradeBase.NotAnalytic
ComradeBase.visibility_point
ComradeBase.visibilitymap_analytic
ComradeBase.visibilitymap_analytic!
ComradeBase.visibilitymap_numeric
ComradeBase.visibilitymap_numeric!
ComradeBase.intensity_point
ComradeBase.intensitymap_analytic
ComradeBase.intensitymap_analytic!
ComradeBase.intensitymap_numeric
ComradeBase.intensitymap_numeric!
```

### Image Types
```@docs
ComradeBase.IntensityMap(::AbstractArray, ::AbstractRectiGrid)
ComradeBase.StokesIntensityMap
ComradeBase.imagepixels
ComradeBase.RectiGrid
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
ComradeBase.header
ComradeBase.ComradeBase.NoHeader
ComradeBase.MinimalHeader
ComradeBase.load
ComradeBase.save
```


## Polarization

```@docs
ComradeBase.AbstractPolarizedModel
```
