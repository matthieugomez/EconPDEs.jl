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
first derivatives of `v` in `k` (see [Model conventions](conventions.md) for the full
naming table). The model chooses between them by the sign of the drift — this is
*upwinding*, explained in [Upwinding](upwinding.md); here capital can drift either way, so
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
    thing to check. See [Solver](solver.md) for how the pseudo-transient scheme works.

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

## Where to go next

- [Model conventions](conventions.md) — the full naming scheme for `u`, grids, and guesses.
- [Upwinding](upwinding.md) — how to choose between `_up` and `_down` derivatives.
- [Boundary conditions](boundary_conditions.md) — reflecting boundaries, the `bc` keyword,
  state constraints, and optimal stopping.
- [Time-dependent problems](time_dependent.md) — solving backward from a terminal value.
- [Solver](solver.md) — how the pseudo-transient scheme works and what to do when a solve
  fails.
- The **Examples** in the sidebar, starting with
  [Neoclassical growth](examples/neoclassical_growth.md) — each solves a published model,
  plots the solution, and explains it economically.
