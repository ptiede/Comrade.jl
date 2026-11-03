```@meta
CurrentModule = Comrade
```

# Writing and using moves

A move is a Metropolis–Hastings transition on the latent point of a posterior. A correct
move leaves the posterior invariant, so it can be interleaved with the transitions of any
MCMC sampler that also leaves the posterior invariant. Moves target directions that
gradient samplers explore slowly, for example a degeneracy between the total flux and the
gain amplitudes, the 2π ambiguity of a gain phase, or the hyperparameters of a gain chain
against the chain values.

## What a move does

A move acts on the base latent point `x` of the posterior, in one of two spaces:

  - the flat space of [`asflat`](@ref), or
  - the standard normal space of [`transport_to`](@ref)`(post, StdNormal())`.

Given `x`, a move draws a step from its step distribution and proposes a point
`x′ = propose(move, x, step, ctx)`. The proposal is accepted with probability `min(1, α)`,

```math
\log \alpha = \ell(x') - \ell(x) + \log\left|\det \frac{\partial x'}{\partial x}\right|,
```

where ``\ell`` is the log density in the base space. This is a valid Metropolis–Hastings
step when the reverse step returns `x′` to `x` and is drawn with the same probability as the
forward step.

In both spaces the log density is ``\ell(x) = \log p(y \mid T(x)) + \log \pi(x)``, where
``T`` maps the latent point to the model parameters and ``\pi`` is the prior density of the
latent point. A move that leaves the likelihood unchanged ([`likelihood_invariant`](@ref)) is
accepted using only the change in ``\log \pi`` and the log-determinant, without evaluating
the likelihood. In the standard normal space ``\log \pi(x)`` is
``-\lVert x \rVert^2/2`` up to a constant, which needs no evaluation of ``T``; in the flat space it is
the prior density of ``T(x)`` times the Jacobian of ``T``. After a call that accepted such a
move, [`MoveSet`](@ref) evaluates ``\ell`` once and throws an error if it differs from the
value implied by the accepted moves. Moves with a [`ComponentwiseNormalStep`](@ref) must
leave the likelihood unchanged.

## Using moves

### Built-in moves

Parameters are addressed by their path, a tuple of field names such as `(:sky, :f1)` or
`(:instrument, :lg)`. A [`CoordinateView`](@ref) maps paths to latent coordinates in a
chosen space. Comrade provides these moves:

  - [`flux_gain_move`](@ref): trades a scalar flux parameter against the gain
    log-amplitudes, keeping `F · exp(power · g)` fixed;
  - [`PhaseSheetMove`](@ref): adds `±2π` to one site's path of a real-line Gauss–Markov
    phase chain, from a free point up to the next point the model fixes;
  - [`chain_hyper_moves`](@ref): one [`ChainHyperMove`](@ref) per fitted hyperparameter of a
    Gauss–Markov chain term, which updates the hyperparameter of each site while keeping the
    chain values fixed;
  - [`CompensatedMove`](@ref): a general random-walk move that shifts one latent coordinate
    and recomputes a block of parameters to compensate (see [Writing a move](@ref)).

### Building a `MoveSet`

A [`MoveSet`](@ref) holds the moves, their step scales and their statistics:

```julia
using Comrade

sp = StdNormal()
view = CoordinateView(post, sp)
fg = flux_gain_move(view; flux = (:sky, :f1), gains = (:instrument, :lg))
sheet = PhaseSheetMove(post, (:gp,); space = sp)
hyper = chain_hyper_moves(view, (:instrument, :lg))

moves = (fg, sheet, hyper...)
ms = MoveSet(
    post, moves; space = sp, θ0 = prior_sample(post),
    rounds = [5; 3; fill(2, length(hyper))],
    initial_scale = [0.05; nothing; fill(0.5, length(hyper))],
)
```

The options are per move. A number applies to every move and a vector gives one value per
move.

  - `rounds`: proposals of each move per call.
  - `initial_scale`: the starting step scale `τ` of moves with a [`NormalStep`](@ref) or
    [`ComponentwiseNormalStep`](@ref). It is ignored for other moves, so their entries may be
    `nothing`.
  - `target_accept`: the acceptance rate each step scale is adapted toward during warmup
    (default `0.45`).
  - `θ0`: a constrained parameter point. If given, construction checks there that every move
    flagged [`likelihood_invariant`](@ref) leaves the likelihood unchanged, and that the
    components of every componentwise move are independent.
  - `rng`: the random number generator for the steps of the `θ0` checks.
  - `output`: a file path. If given, [`move_summary`](@ref) is serialized there after every
    call.

The moves, the `MoveSet`, and the sampler must all use the same latent space. A call
throws an error otherwise.

### With any sampler

Call the `MoveSet` on the latent point between the sampler's steps:

```julia
using Random

function run_chain(ms, tpost, x, rng; nwarmup, nsamples)
    for i in 1:nwarmup
        x = my_sampler_step(x)
        x = ms(x, tpost, :warmup, rng; step = i, total = nwarmup, report = i % 100 == 0)
    end
    for i in 1:nsamples
        x = my_sampler_step(x)
        x = ms(x, tpost, :sampling, rng)
    end
    return x
end

tpost = transport_to(post, sp)
x = run_chain(ms, tpost, prior_sample(tpost), Random.default_rng(); nwarmup = 1000, nsamples = 1000)
```

`tpost` must sample in the space of the `MoveSet`. The call returns `x` itself when no
proposal was accepted, and a new array otherwise; a sampler that caches the log density or
gradient at `x` must recompute them when `x′ !== x`. Each call makes `rounds[j]` proposals of
move `j`. Traceable moves run first, then the others, each group interleaved round by
round.

The phase controls adaptation. During `:warmup` the log step scale of each scaled move
(and of each component of a componentwise move) follows a Robbins–Monro recursion toward
`target_accept`. During `:sampling` the scales are fixed. A warmup call with `report = true`
logs a table of acceptance rates and step scales over the calls since the last table,
labeled with `step` out of `total`.

### With `ReactantNUTS`

Pass the `MoveSet` to `sample` as `between_chunks`, and request the same latent space with
`transport_method`:

```julia
using Reactant

out = sample(
    post_device, ReactantNUTS(; n_adapts = 2000), 2000;
    transport_method = sp, between_chunks = ms,
    saveto = DiskStore(name = "Results", stride = 50),
)
```

Build the `MoveSet` from the host posterior `post`; `post_device` may be its device copy,
e.g. `todevice(post, Reactant.Sharding.NoSharding())`.
The moves run after every warmup chunk (`warmup_chunk`, by default the `DiskStore` stride or
`chunk_size`) and after every sampling chunk except the last, so the chunk length sets how
often they run. The chain stores only the NUTS draws. When the moves change the position,
the next chunk recomputes the gradient and keeps the step size, metric, and adaptation
state.

With a device position, moves that are [`traceable`](@ref) propose inside one compiled
Reactant program per call, and the position stays on the device. Other moves propose on
the host from a copy of the position, with `ctx.view` over the host posterior the `MoveSet`
was built from rather than the device posterior; their log density is evaluated on the
device.

The method `ms(state, tpost, info, rng)` is the `between_chunks` form used by
`ReactantNUTS`. It updates `state.position` and reads `phase`, `step`, `total`, and
`report` from the `NamedTuple` `info`.

### Reading the statistics

```julia
for s in move_summary(ms)
    rate = s.sampling.accepted / s.sampling.proposed
    println(s.name, ": τ = ", s.τ, ", sampling acceptance = ", rate)
end
move_seconds(ms)   # total wall time spent in the moves
```

Each entry of [`move_summary`](@ref) has `name`, the current step scale `τ` (`nothing` for
a move without one, a vector for a componentwise move), `rounds`, the `proposed` and
`accepted` counts in `warmup` and `sampling`, and `bystep`, a dictionary of
`label => [proposed, accepted]` keyed by [`step_label`](@ref) (or by
[`component_label`](@ref) for a componentwise move).

## Writing a move

A move is a subtype of [`AbstractMove`](@ref) that implements three methods:

  - [`propose`](@ref)`(move, x, step, ctx)`: the proposal from the latent point `x`;
  - [`step_distribution`](@ref)`(move)`: one of the three step distributions below;
  - [`move_name`](@ref)`(move)`: a name for logs and statistics, unique within a `MoveSet`.

These functions are not exported, so extend them as `Comrade.propose` and so on.

The step distribution determines the step that `propose` receives and what it must return:

| Step distribution                          | `step` passed to `propose`                     | `propose` returns                      | Reverse step    |
|:------------------------------------------ |:---------------------------------------------- |:-------------------------------------- |:--------------- |
| [`NormalStep`](@ref)`()`                   | scalar `u ~ N(0, τ²)`                          | `(x′, logdet)`                         | `-u`            |
| [`ComponentwiseNormalStep`](@ref)`(n)`     | vector `u` with `uᵢ ~ N(0, τᵢ²)`               | `(x′, logdets)`, one per component     | `-u`            |
| [`SymmetricDiscreteStep`](@ref)`(draw, reverse)` | `draw(rng)`                              | `(x′, logdet)`                         | `reverse(step)` |

Here `logdet = log|det ∂x′/∂x|`. Proposing with the reverse step from `x′` must return `x`
with log-determinant `-logdet`. `propose` must not modify `x`.

A componentwise move also implements [`component_coords`](@ref)`(move)`, the latent
coordinates each component changes. These sets must be disjoint, and component `i` of the
proposal must depend only on `uᵢ` and the coordinates of component `i`, so that each component
is accepted on its own. A componentwise move must be likelihood invariant.

The argument `ctx` is a `NamedTuple`. `ctx.view` is a [`CoordinateView`](@ref) over the
posterior being sampled; use it to read and write parameters by path:

  - [`coords`](@ref)`(ctx.view, path)`: the latent coordinates of a parameter;
  - [`value`](@ref)`(ctx.view, x, path)`: its value at `x`;
  - [`latent`](@ref)`(ctx.view, path, v)`: the latent coordinates that give it the value `v`.

The remaining entries of `ctx` are the arrays returned by the moves'
[`context_data`](@ref).

Optional methods:

  - [`likelihood_invariant`](@ref)`(move)` (default `false`): declares that the move leaves
    the likelihood unchanged, which enables acceptance without the likelihood.
  - [`step_label`](@ref)`(move, step)` (default `nothing`): a key under which proposals and
    acceptances of `step` are counted in `bystep`.
  - [`component_label`](@ref)`(move, i)` (default `i`): the label of component `i` of a
    componentwise move.
  - [`context_data`](@ref)`(move, view)` (default empty): arrays the move reads from `ctx`.
    [`update_context!`](@ref) overwrites them in place, including the copies used by a
    compiled device program.
  - [`traceable`](@ref)`(move)` (default `false`): declares that `propose` compiles with
    Reactant. A traceable move has a `NormalStep` and reads every array it needs
    from `ctx` instead of capturing host objects.

### A compensated move

Suppose the sky model depends on two parameters `a` and `b` only through their product
`a * b`. A [`CompensatedMove`](@ref) can random-walk the latent coordinate of `a` and set
`b′ = b a / a′`, which leaves the likelihood unchanged:

```julia
view = CoordinateView(post, sp)
ab = CompensatedMove(
    "a_b", view, (:sky, :a), (:sky, :b),
    (vb, x, x′, ctx) -> vb * Comrade.value(ctx.view, x, (:sky, :a)) /
        Comrade.value(ctx.view, x′, (:sky, :a));
    likelihood_invariant = true,
)
```

The function `compensate(vb, x, x′, ctx)` receives the value `vb` of the block at `x`, and
`x′`, which is `x` with only the shift applied. It returns the new value of the block. With
log-normal priors on `a` and `b`, the latent coordinate of each is an affine function of its
logarithm in both spaces, so the latent map is a translation and the default
`logdet = 0` is correct. Otherwise pass `logdet = (vb, x, x′, ctx) -> ...`, which receives
the complete proposal `x′`.

### A custom move

If the likelihood depends on a phase `φ` only through `cos φ` and `sin φ`, and `φ` has an
unbounded prior such as a normal distribution, then shifts by `±2π` leave the likelihood
unchanged:

```julia
struct PhaseWrap{P <: Tuple} <: Comrade.AbstractMove
    path::P
end

Comrade.step_distribution(::PhaseWrap) = SymmetricDiscreteStep(rng -> rand(rng, (-1, 1)), -)
Comrade.move_name(::PhaseWrap) = "phase_wrap"
Comrade.likelihood_invariant(::PhaseWrap) = true

function Comrade.propose(m::PhaseWrap, x, s, ctx)
    x′ = copy(x)
    φ = Comrade.value(ctx.view, x, m.path)
    x′[Comrade.coords(ctx.view, m.path)] = Comrade.latent(ctx.view, m.path, φ + 2π * s)
    return x′, 0.0
end

move = PhaseWrap((:sky, :φ))
```

The step `s = ±1` is reversed by `-s`. For a normal prior the latent coordinate is an affine
function of `φ` in both spaces, so the log-determinant is zero.

## Checking a move

Run [`check_move`](@ref) on every new move before using it in a sampler. It takes the move,
the posterior, and a collection of constrained parameter points, and throws an error that
names the first failed property:

```julia
θs = [prior_sample(post) for _ in 1:3]
r = check_move(move, post, θs; space = sp, τ = 0.1, nprior = 400, ngeweke = 1000, geweke_steps = 20)
```

At each point in `θs` it checks:

  - reversal: the reverse step from the proposal returns to the start, with the opposite
    log-determinant;
  - log-determinant: `logdet` agrees with a central finite-difference Jacobian (step `h`)
    over the coordinates the proposal changes;
  - invariance: a move flagged [`likelihood_invariant`](@ref) leaves the log-likelihood
    unchanged;
  - components: for a componentwise move, the components are independent.

Two statistical tests are optional:

  - `nprior > 0`: a prior-invariance test. `nprior` prior draws each take `nsteps`
    Metropolis–Hastings steps of the move targeting the prior, and are compared with fresh
    prior draws on up to `ncoords` changed coordinates by two-sample Kolmogorov–Smirnov tests
    at level `ks_alpha`.
  - `ngeweke > 0`: Geweke's (2004) joint-distribution test of the kernel that `MoveSet`
    runs. Each replicate draws `θ` from the prior and data from `p(y | θ)`
    ([`simulate_observation`](@ref)), then alternates `geweke_steps` move steps targeting
    `p(θ | y)` with fresh draws of the data. The kernel leaves the posterior invariant if and
    only if the results remain draws of the joint distribution `p(θ, y)`. The replicates are
    compared with fresh joint draws on the changed coordinates, on the projection onto the
    move direction, and on the log-likelihood. This test requires a posterior whose data can
    be simulated, which excludes image data.

Steps are drawn at scale `τ`.

The result reports the largest reversal and log-determinant errors, the largest
log-likelihood change, the KS statistics (`ks`, `geweke`) with their critical values, and
the acceptance rate in the joint test (`geweke_accept`). Check that `geweke_accept` is not
close to zero: a move that rarely moves passes the joint test without being tested.

[`check_moveset`](@ref) runs the joint-distribution test on the kernel of a whole
`MoveSet`, at its current step scales and rounds:

```julia
r = check_moveset(ms; ngeweke = 1000, geweke_steps = 20)
all(<=(r.geweke_critical), r.geweke)
```

It returns `geweke`, `geweke_critical`, and `geweke_accept` (one rate per move), and throws
an error if a statistic exceeds its critical value.
