# Writing the PDE function

The PDE function is where the model lives: it encodes the equation at one grid point, and
`pdesolve` assembles it into a global nonlinear system. This page is the reference for
writing it. Everything in `EconPDEs.jl` is keyed by names: the names of the state
variables and the names of the unknown functions. You choose them once, in the grid and
the guess, and they determine the field names of `u`, the required return names, and the
keys of the result. [Getting started](getting_started.md) shows these conventions in
action on a first model.

## Grids and guesses

The state grid is a `NamedTuple` mapping each state variable to an `AbstractVector`; the
guess is a `NamedTuple` mapping each unknown function to an array of starting values with
the same shape as the grid:

```julia
stategrid = (; y = range(0.5, 1.5, length = 10), a = range(0.0, 100.0, length = 500))
guess = (; v = [initial_value(y, a) for y in stategrid[:y], a in stategrid[:a]])
```

Grid entries can be different vector containers and element types: for example, one
dimension can be a `range` and another can be a `Vector`.

Use the `(; name = value)` form; an `OrderedDict` is also accepted. A plain `Dict` is
rejected, because its iteration order is arbitrary and the order of names determines how
the solution arrays are laid out (the first state is the first array dimension, and so
on).

With `N` state variables, each guess array has the same `N`-dimensional shape as the
state grid, and so do the solved unknowns.

## The function signature

The PDE function must have the form

```julia
pde(state::NamedTuple, u::NamedTuple) -> NamedTuple
```

(or `pde(state, u, t)` for a time-dependent equation, see
[Time-dependent problems](time_dependent.md)). `state` gives the current
grid point, one entry per state variable. `u` is the local finite-difference bundle: the
value of each unknown and all the derivative stencils the equation may need at that point.
The model then selects the economically correct stencil, usually by the sign of a drift
(see [Upwinding](#Upwinding)).

## Derivative naming

Field names concatenate the unknown name, the state name or names, and an optional
direction suffix. If the unknown is `v` and the states are `k` and `l`, then `u` contains:

| Field | Meaning |
|---|---|
| `v` | value of the unknown at this grid point |
| `vk_up`, `vk_down` | forward and backward first derivatives in `k` |
| `vl_up`, `vl_down` | forward and backward first derivatives in `l` |
| `vkk`, `vll` | own second derivatives (central) |
| `vkl_up`, `vkl_down` | directional cross derivatives, for positive and negative cross-diffusion terms |

Note there is *no* plain `vk`: the package never chooses a first-derivative direction for
you. Cross-derivative names use the states in grid order (`vkl_up`, not `vlk_up`).

With several unknowns, the same rule applies to each one: unknowns `v` and `w` on states
`k` and `l` give fields such as `v`, `vk_up`, `vkl_up`, `w`, `wk_up`, and `wkl_up`. The argument
name `u` is only a convention; the fields are what matter.

The examples typically unpack the fields they need at the top of the function with Julia's
`NamedTuple` destructuring syntax:

```julia
(; v, vk_up, vk_down, vkk) = u
```

which is shorthand for `v = u.v; vk_up = u.vk_up; …`. The same syntax unpacks the current
grid point (`(; k, l) = state`) and, when the model is a struct, its parameters
(`(; ρ, γ) = m`). Writing `u.vk_up` directly, as in
[Getting started](getting_started.md), is equivalent.

Because names are formed by concatenation, some combinations are ambiguous — states
`(k, a)` with unknowns `(v, vk)` would both generate a field `vka_up`. `pdesolve` checks
for this upfront and throws an error asking you to rename.

This naming step is the one piece of magic in the package, so it is worth knowing what
happens behind it. Before solving, `pdesolve` reads the state names off the grid and the
unknown names off the guess, and generates a function specialized to those names that, at
every grid point, fills each field in the table with the corresponding finite-difference
stencil — one-sided differences from the neighboring grid points for `_up`/`_down`,
central differences for the second derivatives, ghost-node values at the boundaries (see
[Boundary conditions](boundary_conditions.md)). By the time your equation runs, `u` is an
ordinary `NamedTuple` of numbers: every derivative the equation might ask for has already
been computed, and the function simply picks which ones to use.

## The return value: one time derivative per unknown

The PDE function must return one *time derivative* per unknown, named
`Symbol(unknown, :t)` — e.g. `(; vt)` for one unknown, `(; pAt, pBt)` for two. Write the
stationary equation as `0 = RHS - ρ v` and return `vt = -(RHS - ρ v)`:

```julia
vt = -(c^(1 - γ) / (1 - γ) + μk * vk - ρ * u.v)
return (; vt)
```

The minus sign is the package's core convention. `pdesolve` finds the stationary solution
by integrating the *time-dependent* equation ``\partial_t v = \text{vt}`` (a false
transient) until it stops moving, and it is the backward-in-time HJB — with exactly this
sign — whose transient is stable and converges to the stationary solution. If the sign is
flipped, the iteration diverges — `pdesolve` recognizes the typical signature (a residual
Jacobian whose diagonal is negative at every grid point) and warns on the first iteration.
The same convention is what lets a single PDE function serve both stationary solves and
genuinely [time-dependent problems](time_dependent.md).

## Saving intermediate outputs

To store objects computed inside the PDE function (policies, drifts, prices of risk),
return a second `NamedTuple`:

```julia
return (; vt), (; c, μa, r)
```

Each object in the second `NamedTuple` is evaluated at every grid point and stored as an
array with the same shape as the grid, available in `result.saved`.

## The result

`pdesolve` returns an [`EconPDEResult`](api.md) with three fields:

- `solution`: the solved unknowns, a `NamedTuple` of arrays (`result.solution.v`);
- `residual_norm`: the norm of the residual at the solution — check it is small before
  using the output;
- `saved`: the saved objects, together with the solved unknowns for convenience, or an
  empty `NamedTuple` if the PDE saves nothing.

For a time-dependent problem each array gains a trailing time dimension
(`result.solution.v[.., i]` is the solution at time `τs[i]`), and `residual_norm` is a
vector of residual norms.

## Upwinding

`u` gives you two first derivatives per unknown and state — `vk_up` (forward) and
`vk_down` (backward) — and never chooses between them. This section explains why, and gives
the three patterns that cover essentially every model.

### Why upwind at all

A finite-difference scheme for an HJB converges to the correct (viscosity) solution when it
is *monotone* (Barles–Souganidis): loosely, the value at a grid point must be a weighted
average of its neighbors with nonnegative weights, so the discrete solution inherits the
maximum principle of the continuous equation. For the first-derivative (drift) term, this
requires taking the derivative *in the direction the state moves*:

- drift ``\mu > 0`` → use the forward difference (`_up`);
- drift ``\mu < 0`` → use the backward difference (`_down`).

Choosing the wrong direction puts a negative weight on a neighbor. The practical symptom is
not subtle: oscillating iterates, a residual that stalls, or a pseudo-time step that
collapses (see [Solver and troubleshooting](solver.md)). Central
differences for the *second* derivative are always fine (their weights are automatically
positive), which is why `vkk` comes in only one flavor.

Monotonicity also has a probabilistic reading: a monotone scheme is exactly one whose
discretized operator is a valid Markov generator, with non-negative jump rates to
neighboring grid points. That is the condition under which
[InfinitesimalGenerators.jl](https://matthieugomez.github.io/InfinitesimalGenerators.jl/dev/univariate/)
turns the same discretization into a Markov chain — and why an upwinded solution hands off
exactly to its [distributions and expectations](infinitesimal_generators.md).

### Pattern 1: exogenous drift

When the drift of a state does not depend on the derivative you are choosing, upwinding is
one line. For example, income mean-reverting at rate ``\mu_y = \kappa(\bar y - y)``:

```julia
μy = κ * (ybar - y)
vy = (μy >= 0) ? u.vy_up : u.vy_down
```

### Pattern 2: endogenous drift (consumption-saving)

In a consumption-saving problem the asset drift ``\mu_a = y + ra - c`` depends on
consumption, which comes from the first-order condition ``c = (\partial_a v)^{-1/\gamma}``
— which depends on the derivative you are choosing. The standard resolution tries each
direction and keeps the one consistent with itself:

1. compute `c` from `va_up`; if the implied drift is positive, forward was the right
   direction — keep it;
2. otherwise compute `c` from `va_down`; if the implied drift is negative, backward was
   right — keep it;
3. otherwise (the two candidates straddle zero, or the borrowing constraint binds at the
   bottom of the grid), the state does not move: set ``\mu_a = 0``, i.e. consume exactly
   ``c = y + ra``, and read the marginal value off the FOC, ``v_a = c^{-\gamma}``.

```julia
c_up = va_up > 0 ? min(va_up^(-1 / γ), cmax) : cmax
μa_up = y + r * a - c_up
if μa_up >= 0
    va, c, μa = va_up, c_up, μa_up
else
    c_down = va_down > 0 ? min(va_down^(-1 / γ), cmax) : cmax
    μa_down = y + r * a - c_down
    if (μa_down <= 0) && (a > amin)
        va, c, μa = va_down, c_down, μa_down
    else
        c = y + r * a
        va, μa = c^(-γ), 0.0
    end
end
```

Two remarks on the guards:

- `cmax` caps consumption when the marginal value is nonpositive. Newton iterates can
  visit economically meaningless regions (negative `va`) on the way to the solution; the
  cap keeps the FOC well-defined there without affecting the converged answer.
- Case 3 with `a == amin` is precisely how a **borrowing constraint** is imposed: the
  household would like to dissave past the limit, but the constraint pins ``\mu_a = 0``.
  No separate boundary condition is needed — see
  [Boundary conditions](boundary_conditions.md).

The [consumption-saving examples](examples/consumption_saving/consumption_saving_diffusion_income.md)
use exactly this block.

### Pattern 3: general equilibrium (recompute the direction)

In general-equilibrium models the drift of the state is pinned down by equilibrium objects
— prices, volatilities — that themselves depend on the derivatives being solved for. The
logic is the same as Pattern 2, but the derivative now enters through a whole block of
equilibrium computations rather than a single first-order condition, so switching direction
means redoing the block. The pattern: compute everything with `_up` derivatives; if the
resulting drift is negative, recompute once with `_down`:

```julia
pAx, pBx = u.pAx_up, u.pBx_up
iter = 0
@label start
# ... compute σx, prices, and finally the drift μx from pAx, pBx ...
if (iter == 0) && (μx <= 0)
    iter += 1
    pAx, pBx = u.pAx_down, u.pBx_down
    @goto start
end
```

See [Gârleanu-Panageas](examples/asset_pricing/garleanu_panageas.md) for this pattern in
context.

### Cross derivatives

With two or more states and a nonzero *cross*-diffusion term (correlated shocks), the same
monotonicity logic applies to the cross derivative. `u` exposes two directional versions,
`vkl_up` / `vkl_down`. Choose by the sign of whatever multiplies the cross derivative — the
instantaneous covariance of the two states:

```julia
vkl = (σk * σl >= 0) ? u.vkl_up : u.vkl_down
```

### Checking your upwinding

Upwinding on the sign of the raw drift is a shortcut that is guaranteed to produce a
monotone scheme only for terms that are linear in the derivative. As an opt-in diagnostic,
`pdesolve(...; check_monotonicity = true)` inspects the fully assembled scheme and warns
about neighbor weights with the wrong sign — see
[the monotonicity diagnostic](solver.md#The-monotonicity-diagnostic).
