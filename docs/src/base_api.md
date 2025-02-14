```@meta
CurrentModule = StokedBase
```

# Base API



## Index

```@index
Pages = ["base_api.md"]
```



## Model API

```@docs
StokedBase.flux
StokedBase.visibility
StokedBase.visibilitymap
StokedBase.visibilitymap!
StokedBase.intensitymap
StokedBase.intensitymap!
StokedBase.allocate_vismap
StokedBase.allocate_imgmap
StokedBase.create_imgmap
StokedBase.create_vismap
StokedBase.amplitude(::Any, ::Any)
StokedBase.amplitudemap
StokedBase.bispectrum
StokedBase.bispectrummap
StokedBase.closure_phase
StokedBase.closure_phasemap
StokedBase.logclosure_amplitude
StokedBase.logclosure_amplitudemap
```

### Model Interface
```@docs
StokedBase.AbstractModel
StokedBase.AbstractPolarizedModel
StokedBase.visanalytic
StokedBase.imanalytic
StokedBase.ispolarized
StokedBase.radialextent
StokedBase.DensityAnalytic
StokedBase.IsAnalytic
StokedBase.NotAnalytic
StokedBase.visibility_point
StokedBase.visibilitymap_analytic
StokedBase.visibilitymap_analytic!
StokedBase.visibilitymap_numeric
StokedBase.visibilitymap_numeric!
StokedBase.intensity_point
StokedBase.intensitymap_analytic
StokedBase.intensitymap_analytic!
StokedBase.intensitymap_numeric
StokedBase.intensitymap_numeric!
```

### Image Domain
```@docs
StokedBase.imagepixels
StokedBase.RectiGrid
StokedBase.UnstructuredDomain
StokedBase.dims
StokedBase.named_dims
StokedBase.axisdims
StokedBase.domainpoints
StokedBase.fieldofview
StokedBase.pixelsizes
StokedBase.phasecenter
StokedBase.executor
StokedBase.Serial
StokedBase.ThreadsEx
StokedBase.header
StokedBase.NoHeader
StokedBase.MinimalHeader
```

### Image Types
```@docs
StokedBase.IntensityMap
StokedBase.IntensityMap(::AbstractArray, ::AbstractRectiGrid)
StokedBase.UnstructuredMap
StokedBase.baseimage
StokedBase.centroid
StokedBase.second_moment
StokedBase.stokes
```

## Internal Methods not part of public API
```@docs
StokedBase._visibilitymap
StokedBase._visibilitymap!
StokedBase.create_map
```
