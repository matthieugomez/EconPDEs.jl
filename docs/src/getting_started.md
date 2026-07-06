# Getting started

This page solves a first model end to end. The main function is [`pdesolve`](api.md). A
stationary problem needs three ingredients: a **state grid**, an **initial guess**, and a
**PDE function** encoding the equation at one grid point.

## A first model

Here is the deterministic neoclassical growth model,

```math
\rho v(k) = \max_c \left\{ \frac{c^{1-\gamma}}{1-\gamma}
    + v'(k) \left(A k^\alpha - \delta k - c\right) \right\},
```

whose first-order condition is ``c = v'(k)^{-1/\gamma}``.

### The state grid and the guess

We build the grid and the guess first, because their names fix the names used everywhere
else. The grid is a `NamedTuple` whose keys are the state variables (here just `k`); the
guess is a `NamedTuple` whose keys are the unknown functions (here just `v`), holding one
starting value per grid point. The full rules — accepted containers, name ordering, array
shapes — are in [Grids and guesses](pde_function.md#Grids-and-guesses).

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
first derivatives of `v` in `k` ([Derivative naming](pde_function.md#Derivative-naming)
gives the full table). The model chooses between them by the sign of the drift — this is
*upwinding*. Here the drift depends on consumption, which comes from a first-order
condition on the very derivative being chosen, so the code tries each direction and keeps
the one consistent with itself, pinning the drift at zero in between; the `cmax` cap keeps
the first-order condition well-defined when Newton iterates visit regions with nonpositive
marginal value. [Upwinding](pde_function.md#Upwinding) explains why the direction matters
and gives the patterns that cover most models.

```julia
function pde(state::NamedTuple, u::NamedTuple)
    k = state.k
    cmax = 10 * A * k^α

    c_up = u.vk_up > 0 ? min(u.vk_up^(-1 / γ), cmax) : cmax
    μk_up = A * k^α - δ * k - c_up

    if μk_up >= 0
        c, vk, μk = c_up, u.vk_up, μk_up
    else
        c_down = u.vk_down > 0 ? min(u.vk_down^(-1 / γ), cmax) : cmax
        μk_down = A * k^α - δ * k - c_down
        if μk_down <= 0
            c, vk, μk = c_down, u.vk_down, μk_down
        else
            c = A * k^α - δ * k
            vk, μk = c^(-γ), 0.0
        end
    end

    vt = -(c^(1 - γ) / (1 - γ) + μk * vk - ρ * u.v)
    return (; vt)
end
```

The function returns one *time derivative* per unknown, named `Symbol(unknown, :t)` — here
`vt`. For a stationary equation written as `0 = RHS - ρv`, return the negative residual:
`vt = -(RHS - ρv)`. The minus sign is the package's core convention — with the opposite
sign the iteration diverges;
[The return value](pde_function.md#The-return-value:-one-time-derivative-per-unknown)
explains why.

### Solving

```julia
result = pdesolve(pde, stategrid, guess)
result.zero[:v]            # the solved value function, on the grid
result.residual_norm       # should be ≈ 0
```

`pdesolve` returns an [`EconPDEResult`](api.md) with the solved unknowns in `result.zero`,
indexed by name. Check `result.residual_norm` is small before using the output; the full
set of fields is described in [The result](pde_function.md#The-result).

## Exploring the solution

`result.zero` holds the solved unknowns. Objects computed inside the PDE function — the
optimal policy, a drift, an interest rate — can also be stored on the grid by returning a
*second* `NamedTuple`:

```julia
function pde(state::NamedTuple, u::NamedTuple)
    # ... as above ...
    vt = -(c^(1 - γ) / (1 - γ) + μk * vk - ρ * u.v)
    return (; vt), (; c, μk)
end

result = pdesolve(pde, stategrid, guess)
consumption = result.saved[:c]     # same shape as the grid
```

Each saved object has the same shape as the state grid, so it plots directly against it.
The saved drift is often the most useful diagnostic: it shows where the state moves and
where it stops.

```julia
using Plots
plot(stategrid[:k], result.saved[:c]; xlabel = "k", ylabel = "consumption")
plot(stategrid[:k], result.saved[:μk]; xlabel = "k", ylabel = "capital drift")
```

Saved drifts are also the bridge to stationary distributions: once the HJB is solved, the
states follow a known, policy-implied law of motion, and the companion package
[InfinitesimalGenerators.jl](https://github.com/matthieugomez/InfinitesimalGenerators.jl)
turns that law of motion into stationary distributions and expectations. See
[InfinitesimalGenerators](infinitesimal_generators.md) for the workflow.

## Where to go next

- [Writing the PDE function](pde_function.md) — the reference for the naming conventions,
  the sign convention of the return value, and the upwinding patterns.
- [Boundary conditions](boundary_conditions.md),
  [Time-dependent problems](time_dependent.md), and
  [Solver and troubleshooting](solver.md) — what usually matters after the first model
  works.
- The **Examples** in the sidebar, starting with
  [Neoclassical growth](examples/neoclassical_growth.md) — the same model as this page,
  packaged the way the other examples are. Each example solves a published model, plots
  the solution, and explains it economically.
