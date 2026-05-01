# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Comrade.jl is a Julia package for Bayesian modeling of Very-Long-Baseline Interferometry (VLBI) radio data (e.g. EHT). It provides composable sky models, instrument/calibration models, likelihoods, priors, and a unified `VLBIPosterior` exposed via the `DensityInterface` and `LogDensityProblems` interfaces — sampler/optimizer backends are wired in through weak-dep extensions.

Supported Julia versions: **1.10 (LTS) and 1.11**. Newer versions are not officially supported due to ongoing Enzyme compatibility work.

## Common commands

All commands run from the repo root unless noted.

### Tests

```bash
# Full test suite (uses test/Project.toml; matches CI)
julia --project=. -e 'using Pkg; Pkg.test()'

# Faster iteration on a single test file — activate the test env once, then include the file repeatedly in a live REPL
julia --project=test
# in REPL:
julia> using Pkg; Pkg.develop(path="."); Pkg.instantiate()
julia> include("test/Core/models.jl")     # or test/ext/comradeahmc.jl, etc.
```

`test/runtests.jl` is the entry point. It loads `test/test_util.jl` (shared `load_data`, `test_model`, `test_prior`) and includes:
- `test/Core/core.jl` — pure-Comrade tests (observation, models, partially_fixed, imgnormal, bayes)
- `test/ext/*.jl` — one file per weak-dep extension (AHMC, Optimization, Pigeons, Dynesty, Nested)

### Sub-packages under `lib/`

`lib/Comrade{AHMC,AdaptMCMC,Dynesty,Nested,Optimization}` are separate registered packages with their own `Project.toml` and `test/runtests.jl`. Test one with:

```bash
julia --project=lib/ComradeAHMC -e 'using Pkg; Pkg.test()'
```

Note: these libs `Pkg.develop(path="../..")` against the parent Comrade in dev workflows — keep them in sync when changing public APIs.

### Formatting

CI enforces **Runic** formatting (`.github/workflows/formatting.yml`). Format before pushing:

```bash
julia -e 'using Pkg; Pkg.add("Runic")'      # one-time
julia -e 'using Runic; Runic.format_file("path/to/file.jl"; inplace=true)'
# or run on the whole repo via the runic CLI if installed
```

### Docs

```bash
julia --project=docs docs/make.jl
```

## Architecture

The top-level `Comrade.jl` module wires together five subsystems via `include`s. Reading `src/Comrade.jl` first gives the dependency order.

### Layered package stack (re-exported)

Comrade re-exports and builds on a stack of sister packages — many "Comrade" symbols actually live elsewhere:

- **ComradeBase** — abstract types (`AbstractModel`, `AbstractPolarizedModel`, `AbstractDomain`, `AbstractRectiGrid`, `UnstructuredDomain`), trait types (`IsAnalytic`/`NotAnalytic`, `IsPolarized`/`NotPolarized`), and the `visibilitymap`/`intensitymap` interface that subtypes implement.
- **VLBISkyModels** — concrete sky models, Fourier transforms, `FourierDualDomain`, `ContinuousImage`.
- **PolarizedTypes** — Stokes/coherency types.
- **VLBIImagePriors** — image-space priors (Markov random fields, Dirichlet, etc.).
- **VLBILikelihoods** — likelihood distributions for closure phases, log-closure amplitudes, complex visibilities.
- **HypercubeTransform** / **TransformVariables** — bijections between unconstrained ℝⁿ / unit cube and the constrained prior space; `ascube` and `asflat` produce the transformed posterior used by samplers.

When something looks "missing" from this repo, check those upstream packages.

### Subsystems (under `src/`)

1. **`observations/`** — observation data containers. `EHTObservationTable` and `AbstractArrayConfiguration` describe the array, baselines, and measured data products (`VisData`, `DualData`, closure quantities). Includes `obstable.jl`, `timetable.jl`, and operations like coherency↔Stokes conversion.
2. **`skymodels/`** — `AbstractSkyModel` and `SkyModel` (the canonical concrete type). A `SkyModel` bundles a function `f(params) -> AbstractModel`, a prior over `params`, and an image-domain grid. `set_array(skymodel, array)` produces an `ObservedSkyModel` with a `FourierDualDomain` ready for fast visibility evaluation. See `src/skymodels/abstract.jl` for the trait-based dispatch on `idealmaps`/`forward_dualmap`.
3. **`instrument/`** — calibration / instrument model. `JonesMatrix`-based corruption model with site-level priors, integration-time and frequency-channel binning, `SiteArray` storage, and `caltable` output. The instrument model has its own parameter transform (`instrument_transforms.jl`) that composes with the sky transform.
4. **`posterior/`** — `VLBIPosterior(skymodel, instrumentmodel, dataproducts...)` ties it all together. `AbstractVLBIPosterior` implements `DensityInterface` and `LogDensityProblems`. `transformed.jl` wraps a posterior with an `ascube`/`asflat` transform so samplers see an unconstrained/uniform domain. Key entry points: `logdensityof`, `logprior`, `loglikelihood`, `forward_model`, `prior_sample`, `simulate_observation`, `chi2`, `residuals`.
5. **`inference/`** — sampler-agnostic `AbstractMCMC` glue, `MemoryStore`/`DiskStore` checkpointing for long chains, `PosteriorSamples` container, optimization helpers, and `load_samples` for reloading on-disk results.

`dirty_image.jl` and `mrf_image.jl` are utilities for initial-guess imaging. `rules.jl` holds ChainRules / Enzyme rules for AD correctness.

### Extensions (`ext/`)

Sampler and visualization backends are weak deps, loaded on demand:

| Extension | Trigger | Purpose |
|---|---|---|
| `ComradeAdvancedHMCExt` | `using AdvancedHMC` | NUTS / HMC sampling |
| `ComradeOptimizationExt` | `using Optimization` | MAP / MLE via Optimization.jl |
| `ComradeDynestyExt` | `using Dynesty` | Nested sampling |
| `ComradePigeonsExt` | `using Pigeons` | Parallel tempering |
| `ComradeEnzymeExt` | `using Enzyme` | AD rules / shadows for Enzyme |
| `ComradeMakieExt` | `using Makie` | Plotting recipes |
| `ComradePyehtimExt` | `using Pyehtim` | Read EHT-format `.uvfits` via Python's `ehtim` |

If you change a public API consumed inside an `ext/` file, also update the corresponding `lib/Comrade<X>` sub-package — they target the same backends and tests run separately in CI.

## Conventions

- **AD**: Enzyme is the primary AD; the `Enzyme` compat bound in `Project.toml` is intentionally narrow because of known regressions — don't widen it without testing. Some functions carry explicit `EnzymeRules.inactive` or ChainRules definitions in `src/rules.jl`; respect these when refactoring.
- **Parameter NamedTuples**: posterior parameters are `NamedTuple`s with `:sky` and (optionally) `:instrument` fields. Sky-only models pass `θ.sky`; calibrated models add `θ.instrument`.
- **Domain types**: prefer `AbstractRectiGrid` for image-space and `UnstructuredDomain` for visibility-space rather than raw arrays — Fourier transforms dispatch on these.
- **Re-exports**: prefer `using Comrade` in user-facing code; internal code may import directly from `ComradeBase` / `VLBISkyModels` for clarity.
