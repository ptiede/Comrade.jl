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
Comrade.CalTable
Comrade.caltable(::Comrade.JonesCache, ::AbstractVector)
Comrade.caltable(::Comrade.EHTObservation, ::AbstractVector)
Comrade.DesignMatrix
Comrade.JonesCache
Comrade.TransformCache
Comrade.CorruptionModel
Comrade.VLBIModel
Comrade.CalPrior
Comrade.CalPrior(::NamedTuple, ::JonesCache)
Comrade.CalPrior(::NamedTuple, ::NamedTuple, ::JonesCache)
Comrade.RIMEModel
Comrade.ObsSegmentation
Comrade.IntegSeg
Comrade.ScanSeg
Comrade.TrackSeg
Comrade.FixedSeg
Comrade.jonescache(::Comrade.EHTObservation, ::Comrade.ObsSegmentation)
Comrade.jonesStokes
Comrade.jonesG
Comrade.jonesD
Comrade.jonesT
Comrade.PoincareSphere2Map
Comrade.caltable
Comrade.JonesPairs
Comrade.GainSchema
```


### Polarized Models

```@docs
ComradeBase.mbreve
ComradeBase.evpa(pimg::ComradeBase.AbstractPolarizedModel, p)
```



## Data Types


```@docs
Comrade.extract_table
Comrade.ComplexVisibilities
Comrade.VisibilityAmplitudes
Comrade.ClosurePhases
Comrade.LogClosureAmplitudes
Comrade.Coherencies
Comrade.amplitude(::Comrade.EHTVisibilityDatum)
Comrade.amplitude(::Comrade.EHTVisibilityAmplitudeDatum)
Comrade.baselines
Comrade.arrayconfig
Comrade.closure_phase(::Comrade.EHTVisibilityDatum, ::Comrade.EHTVisibilityDatum, ::Comrade.EHTVisibilityDatum)
Comrade.getdata
Comrade.getuv
Comrade.getuvtimefreq
Comrade.scantable
Comrade.stations
Comrade.uvpositions
Comrade.ArrayConfiguration
Comrade.ClosureConfig
Comrade.AbstractInterferometryDatum
Comrade.ArrayBaselineDatum
Comrade.EHTObservation
Comrade.EHTArrayConfiguration
Comrade.EHTCoherencyDatum
Comrade.EHTClosurePhaseDatum
Comrade.EHTLogClosureAmplitudeDatum
Comrade.EHTVisibilityDatum
Comrade.EHTVisibilityAmplitudeDatum
Comrade.Scan
Comrade.ScanTable
```

## eht-imaging interface (Internal)

```@docs
Comrade.extract_amp
Comrade.extract_cphase
Comrade.extract_lcamp
Comrade.extract_vis
```

## Bayesian Tools

### Posterior Constructions

```@docs
Comrade.ascube
Comrade.asflat
Comrade.flatten
Comrade.inverse
Comrade.prior_sample
Comrade.likelihood
Comrade.simulate_observation
Comrade.dataproducts
Comrade.skymodel
Comrade.instrumentmodel
Comrade.vlbimodel
Comrade.sample(::Posterior)
Comrade.transform
Comrade.MultiRadioLikelihood
Comrade.Posterior
Comrade.TransformedPosterior
Comrade.RadioLikelihood
Comrade.IsFlat
Comrade.IsCube
```

## Misc

```@docs
Comrade.μas2rad
Comrade.rad2μas
Comrade.load
Comrade.save
Comrade.NonAnalyticTest
Comrade.station_tuple
Comrade.center_image
Comrade.convolve!
Comrade.convolve
Comrade.dirty_image
Comrade.dirty_beam
```

## Internal (Not Public API)
```@docs
VLBISkyModels.scale_uv
VLBISkyModels.scale_image
VLBISkyModels.transform_uv
VLBISkyModels.transform_image
VLBISkyModels.ThreadedModel
Comrade.extract_FRs
```

