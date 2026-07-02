# Getting started

`EconPDEs.jl` solves nonlinear ODEs and PDEs that arise in economic models, especially
Hamilton–Jacobi–Bellman equations. You write the local equation; the package supplies the
finite-difference derivatives, upwinding, sparse Jacobians, and pseudo-transient Newton
iteration.

## Installation

The package is registered in the Julia `General` registry:

```julia
] add EconPDEs
```

Current versions require Julia 1.10 or later.

## A first model

The main function is [`pdesolve`](api.md). A stationary problem needs three ingredients: a
**state grid**, an **initial guess**, and a **PDE function** encoding the equation at one
grid point. Here is the deterministic neoclassical growth model,

```math
\rho v(k) = \max_c \left\{ \frac{c^{1-\gamma}}{1-\gamma}
    + v'(k) \left(A k^\alpha - \delta k - c\right) \right\},
```

whose first-order condition is ``c = v'(k)^{-1/\gamma}``.

### The state grid and the guess

We build the grid and the guess first, because their names fix the names used everywhere
else. The grid is a `NamedTuple` whose keys are the state variables (here just `k`); the
guess is a `NamedTuple` whose keys are the unknown functions (here just `v`), holding one
starting value per grid point:

```julia
using EconPDEs

const A = 0.5
const α = 0.3
const δ = 0.05
const ρ = 0.05
const γ = 2.0

# Steady-state capital satisfies α A k^(α - 1) = ρ + δ.
kbar = (α * A / (ρ + δ))^(1 / (1 - α))
stategrid = (; k = range(0.1 * kbar, 5.0 * kbar, length = 200))

guess = (; v = [(A * k^α)^(1 - γ) / (1 - γ) / ρ for k in stategrid[:k]])
```

### The PDE function

The PDE function receives the current grid point `state` and a `NamedTuple` `u` holding
each unknown together with its finite-difference derivatives at that point. Derivative
fields are named by concatenation: `u.vk_up` and `u.vk_down` are the forward and backward
first derivatives of `v` in `k` (see [Model conventions](#Model-conventions) for the full
naming table). The model chooses between them by the sign of the drift — this is
*upwinding*, explained [below](#Upwinding); here capital can drift either way, so
we try both directions:

```julia
function hjb(state::NamedTuple, u::NamedTuple)
    k = state.k

    c_up = max(u.vk_up, eps())^(-1 / γ)
    μ_up = A * k^α - δ * k - c_up

    c_down = max(u.vk_down, eps())^(-1 / γ)
    μ_down = A * k^α - δ * k - c_down

    if μ_up > 0
        c, vk, μk = c_up, u.vk_up, μ_up
    elseif μ_down < 0
        c, vk, μk = c_down, u.vk_down, μ_down
    else
        c, vk, μk = A * k^α - δ * k, u.vk_up, 0.0
    end

    vt = -(c^(1 - γ) / (1 - γ) + μk * vk - ρ * u.v)
    return (; vt)
end
```

The function returns one *time derivative* per unknown, named `Symbol(unknown, :t)` — here
`vt`. The sign convention matters:

!!! note "Why `vt = -(...)`?"
    `pdesolve` treats the stationary equation as the long-run limit of the time-dependent
    equation ``\partial_t v = \text{vt}`` and integrates this *false transient* until it
    stops moving. For an HJB written as
    ``\rho v = u(c) + v'(k)\mu_k``, the time-dependent version is the backward equation
    ``\rho v = u(c) + v'(k)\mu_k + \partial_t v``, so the PDE function must return
    ``\text{vt} = -\left(u(c) + v'(k)\mu_k - \rho v\right)``.
    With this sign, the false transient is stable — each pseudo-time step is a contraction
    and the iteration converges to the stationary solution. Flip the sign and the iteration
    runs *away* from the solution: if a model diverges immediately, this sign is the first
    thing to check. See [Solving models](solving.md#Solver-and-troubleshooting) for how
    the pseudo-transient scheme works.

### Solving

```julia
result = pdesolve(hjb, stategrid, guess)
result.zero[:v]            # the solved value function, on the grid
result.residual_norm       # should be ≈ 0
```

`pdesolve` returns an [`EconPDEResult`](api.md) with the solved unknowns in `result.zero`,
indexed by name.

## Saving policies and other outputs

Objects computed inside the PDE function — the optimal policy, a drift, an interest rate —
can be stored on the grid by returning a *second* `NamedTuple`:

```julia
function hjb(state::NamedTuple, u::NamedTuple)
    # ... as above ...
    vt = -(c^(1 - γ) / (1 - γ) + μk * vk - ρ * u.v)
    return (; vt), (; c, μk)
end

result = pdesolve(hjb, stategrid, guess)
consumption = result.optional[:c]     # same shape as the grid
```

Each saved object has the same shape as the state grid, so it plots directly against it:

```julia
using Plots
plot(stategrid[:k], result.optional[:c]; xlabel = "k", ylabel = "consumption")
```

## Model conventions

Everything in `EconPDEs.jl` is keyed by names: the names of the state variables and the
names of the unknown functions. You choose them once, in the grid and the guess, and they
determine the field names of `u`, the required return names, and the keys of the result.

### Grids and guesses

The state grid is a `NamedTuple` mapping each state variable to an `AbstractVector`; the
guess is a `NamedTuple` mapping each unknown function to an array of starting values with
the same shape as the grid:

```julia
stategrid = (; y = range(0.5, 1.5, length = 10), a = range(0.0, 100.0, length = 500))
yend = (; v = [initial_value(y, a) for y in stategrid[:y], a in stategrid[:a]])
```

Grid entries can be different vector containers and element types: for example, one
dimension can be a `range` and another can be a `Vector`.

Use the `(; name = value)` form; an `OrderedDict` is also accepted. A plain `Dict` is
rejected, because its iteration order is arbitrary and the order of names determines how
the solution arrays are laid out (the first state is the first array dimension, and so
on).

With one, two, or three state variables, the guess is a `Vector`, `Matrix`, or
3-dimensional `Array`, and so are the solved unknowns.

### The PDE function

The PDE function encodes the equation at one grid point. It must have the form

```julia
f(state::NamedTuple, u::NamedTuple) -> NamedTuple
```

(or `f(state, u, t)` for a time-dependent equation, see
[Time-dependent problems](solving.md#Time-dependent-problems)). `state` gives the current
grid point, one entry per state variable. `u` is the local finite-difference bundle: the
value of each unknown and all the derivative stencils the equation may need at that point.
The model then selects the economically correct stencil, usually by the sign of a drift
(see [Upwinding](#Upwinding)).

### Derivative naming

Field names concatenate the unknown name, the state name or names, and an optional
direction suffix. If the unknown is `v` and the states are `k` and `l`, then `u` contains:

| Field | Meaning |
|---|---|
| `v` | value of the unknown at this grid point |
| `vk_up`, `vk_down` | forward and backward first derivatives in `k` |
| `vl_up`, `vl_down` | forward and backward first derivatives in `l` |
| `vkk`, `vll` | own second derivatives (central) |
| `vkl` | central cross derivative |
| `vkl_up`, `vkl_down` | directional cross derivatives, for positive and negative cross-diffusion terms |

Note there is *no* plain `vk`: the package never chooses a first-derivative direction for
you. Cross-derivative names use the states in grid order (`vkl`, not `vlk`).

With several unknowns, the same rule applies to each one: unknowns `v` and `w` on states
`k` and `l` give fields such as `v`, `vk_up`, `vkl`, `w`, `wk_up`, and `wkl`. The argument
name `u` is only a convention; the fields are what matter.

Because names are formed by concatenation, some combinations are ambiguous — states
`(k, a)` with unknowns `(v, vk)` would both generate a field `vka_up`. `pdesolve` checks
for this upfront and throws an error asking you to rename.

### The return value: one time derivative per unknown

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
flipped, the iteration diverges. The same convention is what lets a single PDE function
serve both stationary solves and genuinely
[time-dependent problems](solving.md#Time-dependent-problems).

### Saving intermediate outputs

To store objects computed inside the PDE function (policies, drifts, prices of risk),
return a second `NamedTuple`:

```julia
return (; vt), (; c, μa, r)
```

Each object in the second `NamedTuple` is evaluated at every grid point and stored as an
array with the same shape as the grid, available in `result.optional`.

### The result

`pdesolve` returns an [`EconPDEResult`](api.md) with three fields:

- `zero`: the solved unknowns, indexed by name (`result.zero[:v]`);
- `residual_norm`: the norm of the residual at the solution — check it is small before
  using the output;
- `optional`: the saved objects, together with the solved unknowns for convenience.

For a time-dependent problem, `zero` is instead a vector of solutions, one per time point
(`result.zero[i][:v]`), and `residual_norm` a vector of residual norms.

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
collapses (see [Solving models](solving.md#Solver-and-troubleshooting)). Central
differences for the *second* derivative are always fine (their weights are automatically
positive), which is why `vkk` comes in only one flavor.

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
  [Solving models](solving.md#Boundary-conditions).

The [consumption-saving examples](examples/consumption_saving/consumption_saving_diffusion_income.md)
use exactly this block.

### Pattern 3: endogenous volatility (recompute the direction)

In general-equilibrium models the *drift of the state itself* can depend on the derivatives
being solved for (through endogenous volatility or prices). Then the upwind direction is
not known until the drift is computed. The pattern: compute everything with `_up`
derivatives; if the resulting drift is negative, recompute once with `_down`:

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
monotonicity logic applies to the cross derivative. `u` exposes three versions: the central
`vkl` (second-order accurate but not monotone) and the directional `vkl_up` / `vkl_down`.
Choose by the sign of whatever multiplies the cross derivative — the instantaneous
covariance of the two states:

```julia
vkl = (σk * σl >= 0) ? u.vkl_up : u.vkl_down
```

### Checking your upwinding

Upwinding on the sign of the raw drift is a shortcut that is guaranteed to produce a
monotone scheme only for terms that are linear in the derivative. As an opt-in diagnostic,
`pdesolve(...; check_monotonicity = true)` inspects the fully assembled residual Jacobian —
after endogenous policies, risk adjustments, and nonlinear terms — and warns about
neighbor weights with the wrong sign, naming the state, the neighbor, and the stencil that
most likely caused it. For transformed equations, treat the warning as a debugging aid
rather than proof that the economic upwind rule is wrong. See
[Solving models](solving.md#Solver-and-troubleshooting) for the related tolerances.

## Where to go next

- [Solving models](solving.md) — boundary conditions, time-dependent problems, and solver
  troubleshooting.
- The **Examples** in the sidebar, starting with
  [Neoclassical growth](examples/neoclassical_growth.md) — each solves a published model,
  plots the solution, and explains it economically.
