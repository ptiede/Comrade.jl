# Comrade API



## Contents

```@contents
Pages = ["api.md"]
```

## Index

```@index
Pages = ["api.md"]
```

## Model Definitions

### Calibration Models

```@docs
Comrade.corrupt
Comrade.DesignMatrix
Comrade.GainCache
Comrade.GainCache(::Comrade.ScanTable)
Comrade.GainModel
Comrade.GainPrior
Comrade.GainPrior(::Any, ::Comrade.ScanTable)
Comrade.RIMEModel
```

### Combinators

```@docs
Base.:+(::Comrade.AbstractModel, ::Comrade.AbstractModel)
Comrade.added
Comrade.convolved
Comrade.components
Comrade.smoothed
Comrade.CompositeModel
Comrade.AddModel
Comrade.ConvolvedModel
```

### GeometricModels

```@docs
Comrade.scale_uv
Comrade.scale_image
Comrade.transform_uv
Comrade.transform_image
Comrade.GeometricModel
Comrade.ConcordanceCrescent
Comrade.Crescent
Comrade.Disk
Comrade.ExtendedRing
Comrade.Gaussian
Comrade.MRing
Comrade.Ring
Comrade.ParabolicSegment
```

### Model Image (non analytic FFT)

```@docs
Comrade.create_cache
Comrade.update_cache
Comrade.modelimage
Comrade.uviterator
Comrade.fouriermap
Comrade.ModelImage
Comrade.DFTAlg
Comrade.DFTAlg(::Comrade.EHTObservation)
Comrade.DFTAlg(::AbstractArray. ::AbstractArray)
Comrade.DFTAlg(::Comrade.ArrayConfiguration)
Comrade.FFTAlg
Comrade.FFTCache
Comrade.NFFTAlg
Comrade.NUFTCache
Comrade.ObservedNUFT
```


### Modifiers

```@docs
Comrade.basemodel
Comrade.renormed
Comrade.rotated
Comrade.posangle
Comrade.shifted
Comrade.stretched
Comrade.AbstractModifier
Comrade.RenormalizedModel
Comrade.RotatedModel
Comrade.ShiftedModel
Comrade.StretchedModel
```

### Polarized Models

```@docs
Comrade.m̆(pimg::ComradeBase.AbstractPolarizedModel, u, v)
Comrade.mbreve
Comrade.evpa(pimg::ComradeBase.AbstractPolarizedModel, u, v)
Comrade.coherencymatrix(pimg::PolarizedModel, u, v)
Comrade.PolarizedModel
```


### Model Evaluation

For more docstrings on how to evaluate models see [ComradeBase API](@ref).

```@docs
Comrade.amplitude
Comrade.amplitudes
Comrade.bispectra
Comrade.bispectrum
Comrade.closure_phase
Comrade.closure_phases
Comrade.logclosure_amplitude
Comrade.logclosure_amplitudes
Comrade.visibilities
Comrade.visibility
Comrade.visibility_point
Comrade.intensitymap
Comrade.intensitymap!
Comrade.intensity_point
```

## Data Types


```@docs
Comrade.amplitude(::Comrade.EHTVisibilityDatum)
Comrade.amplitude(::Comrade.EHTVisibilityAmplitudeDatum)
Comrade.baselines
Comrade.arrayconfig
Comrade.closure_phase(::Comrade.EHTVisibilityDatum, ::Comrade.EHTVisibilityDatum, ::Comrade.EHTVisibilityDatum)
Comrade.getdata
Comrade.getuv
Comrade.getuvtimefreq
Comrade.rescaleuv!
Comrade.scantable
Comrade.stations
Comrade.uvpositions
Comrade.ClosureConfig
Comrade.ArrayBaselineDatum
Comrade.EHTObservation
Comrade.EHTArrayConfiguration
Comrade.EHTClosurePhaseDatum
Comrade.EHTLogClosureAmplitudeDatum
Comrade.EHTVisibilityDatum
Comrade.EHTVisibilityAmplitudeDatum
Comrade.Scan
Comrade.ScanTable
```

## eht-imaging interface

```@docs
Comrade.extract_amp
Comrade.extract_cphase
Comrade.extract_lcamp
Comrade.extract_vis
Comrade.load_ehtim
```

## Bayesian Tools

### Distributions

```@docs
Comrade.AmpNormal
Comrade.ComplexNormal
Comrade.CPVonMises
```

### Posterior Constructions

```@docs
Comrade.ascube
Comrade.asflat
Comrade.flatten
Comrade.inverse
Comrade.prior_sample
Comrade.sample(::Posterior)
Comrade.transform
Comrade.MultiRadioLikelihood
Comrade.Posterior
Comrade.TransformedPosterior
Comrade.RadioLikelihood
```

## Misc

```@docs
Comrade.μas2rad
Comrade.rad2μas
Comrade.fileio_load
Comrade.fileio_save
Comrade.make_pullback
```
