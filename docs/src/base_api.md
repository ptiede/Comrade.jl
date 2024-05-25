```@meta
CurrentModule = ComradeBase
```

# Base API



## Index

```@index
Pages = ["base_api.md"]
```



## Model API

```@docs
ComradeBase.flux
ComradeBase.visibility
ComradeBase.visibilitymap
ComradeBase.visibilitymap!
ComradeBase.intensitymap
ComradeBase.intensitymap!
ComradeBase.allocate_vismap
ComradeBase.allocate_imgmap
ComradeBase.create_imgmap
ComradeBase.create_vismap
ComradeBase.amplitude(::Any, ::Any)
ComradeBase.amplitudemap
ComradeBase.bispectrum
ComradeBase.bispectrummap
ComradeBase.closure_phase
ComradeBase.closure_phasemap
ComradeBase.logclosure_amplitude
ComradeBase.logclosure_amplitudemap
```

### Model Interface
```@docs
ComradeBase.AbstractModel
ComradeBase.AbstractPolarizedModel
ComradeBase.visanalytic
ComradeBase.imanalytic
ComradeBase.ispolarized
ComradeBase.radialextent
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

### Image Domain
```@docs
ComradeBase.imagepixels
ComradeBase.RectiGrid
ComradeBase.UnstructuredDomain
ComradeBase.dims
ComradeBase.named_dims
ComradeBase.axisdims
ComradeBase.domainpoints
ComradeBase.fieldofview
ComradeBase.pixelsizes
ComradeBase.phasecenter
ComradeBase.executor
ComradeBase.Serial
ComradeBase.ThreadsEx
ComradeBase.header
ComradeBase.NoHeader
ComradeBase.MinimalHeader
```

### Image Types
```@docs
ComradeBase.IntensityMap
ComradeBase.IntensityMap(::AbstractArray, ::AbstractRectiGrid)
ComradeBase.UnstructuredMap
ComradeBase.baseimage
ComradeBase.centroid
ComradeBase.second_moment
ComradeBase.load
ComradeBase.save
ComradeBase.stokes
```

## Internal Methods not part of public API
```@docs
ComradeBase._visibilitymap
ComradeBase._visibilitymap!
ComradeBase.create_map
```
