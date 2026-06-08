"""
    ReactantNUTS(; kwargs...)

Experimental NUTS sampler that runs through Reactant's ProbProg `mcmc_logpdf`.
Like AdvancedHMC's `NUTS`, all sampler configuration lives in the struct.

Warmup is a **single fused** `mcmc_logpdf` call over `n_adapts` steps. The ProbProg
NUTS pass runs Stan's windowed adaptation *internally* (fixed `init_buffer = 75`,
`term_buffer = 50`, base window `25`, doubling, with the same 15%/10% fallback when
the buffers don't fit) — one continuous Nesterov dual-averaging step-size process
(γ=0.05, t₀=10, κ=0.75, target accept 0.8) plus a Welford diagonal mass matrix
updated at window ends. This is algorithmically identical to AdvancedHMC's NUTS
with its `StanHMCAdaptor`.

# Keyword arguments
  n_adapts = 1000          total warmup steps; the internal Stan windowed schedule
                           is derived from this by the ProbProg backend.
  init_step_size = 0.01    initial leapfrog step size; seeds dual averaging (see note)
  max_tree_depth = 10
  max_delta_energy = 1000.0
  strong_zero = true       turn 0*Inf / 0*NaN in the gradient into 0; REQUIRED for
                           stiff image models or every proposal NaN-rejects and the
                           chain freezes -- keep `true` unless you know better.

!!! note "Differences from AdvancedHMC's `NUTS`"
    Two pieces of AdvancedHMC's setup are not (yet) reachable through ProbProg:
    - **No `find_reasonable_step_size`.** AdvancedHMC seeds dual averaging from
      `find_good_stepsize`; ProbProg has no such heuristic and instead anchors dual
      averaging at `init_step_size` (its prox-center is `log(10*init_step_size)`).
      Pick `init_step_size` in the right ballpark for your problem.
    - **Target acceptance is fixed at 0.8** in the backend (matching `NUTS(0.8)`'s
      default) and is not currently configurable.
"""
Base.@kwdef struct ReactantNUTS
    n_adapts::Int = 1000
    init_step_size::Float64 = 0.01
    max_tree_depth::Int = 10
    max_delta_energy::Float64 = 1000.0
    strong_zero::Bool = true
end
