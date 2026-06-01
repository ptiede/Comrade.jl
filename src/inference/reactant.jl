"""
    ReactantNUTS(; kwargs...)

Experimental NUTS sampler that runs through Reactant's ProbProg `mcmc_logpdf`.
Like AdvancedHMC's `NUTS`, all sampler configuration lives in the struct; the
warmup adaptation **schedule** (Stan-style windowed/doubling) is part of it too.

# Keyword arguments
  n_adapts = 1000          total warmup steps; the windowed schedule is DERIVED from
                           this (Stan-style), not hardcoded.
  init_step_size = 0.01    initial leapfrog step size (adapted during warmup)
  max_tree_depth = 10
  max_delta_energy = 1000.0
  strong_zero = true       turn 0*Inf / 0*NaN in the gradient into 0; REQUIRED for
                           stiff image models or every proposal NaN-rejects and the
                           chain freezes -- keep `true` unless you know better.

  # --- warmup schedule (Stan's windowed adaptation; all derived from n_adapts) ---
  init_buffer = 75         fixed-metric step-size warmup at the start
  term_buffer = 50         fixed-metric step-size warmup at the end
  base_window = 25         size of the first metric-estimation window; each subsequent
                           window doubles, and the last absorbs the remainder.
  windows = nothing        advanced override: pass an explicit collection of window
                           sizes to bypass the derived schedule (then init/term are
                           used verbatim).

The derived schedule exactly matches Stan: if `init_buffer + base_window +
term_buffer > n_adapts` the buffers fall back to 15%/75%/10% of `n_adapts` with a
single window.

!!! note "Not yet exposed"
    `mcmc_logpdf` does not currently take a target acceptance rate, so unlike
    `NUTS(0.8)` there is no `target_accept` knob here.
"""
Base.@kwdef struct ReactantNUTS{W}
    n_adapts::Int = 1000
    init_step_size::Float64 = 0.01
    max_tree_depth::Int = 10
    max_delta_energy::Float64 = 1000.0
    strong_zero::Bool = true
    init_buffer::Int = 75
    term_buffer::Int = 50
    base_window::Int = 25
    windows::W = nothing
end
