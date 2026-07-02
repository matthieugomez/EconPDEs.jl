[![Build status](https://github.com/matthieugomez/EconPDEs.jl/workflows/CI/badge.svg)](https://github.com/matthieugomez/EconPDEs.jl/actions)
[![Coverage Status](http://codecov.io/github/matthieugomez/EconPDEs.jl/coverage.svg?branch=main)](http://codecov.io/github/matthieugomez/EconPDEs.jl/?branch=main)

# EconPDEs.jl

`EconPDEs.jl` solves nonlinear ODEs and PDEs that arise in economic models, especially Hamilton-Jacobi-Bellman equations. You write the local equation; the package supplies finite-difference derivatives, upwinding, sparse Jacobians, and pseudo-transient Newton iteration.

Use it when an economic model gives a stationary or time-dependent HJB on a grid, possibly with several value functions and several state variables.

## Installation

The package is registered in the Julia `General` registry:

```julia
] add EconPDEs
```

Current versions require Julia 1.10 or later.

## Quickstart

The main function is `pdesolve`. A stationary problem needs a PDE function, a state grid, and an initial guess. Here is a simple deterministic neoclassical growth model:

```math
\rho v(k) = \max_c \left\{ \frac{c^{1-\gamma}}{1-\gamma}
    + v'(k) \left(A k^\alpha - \delta k - c\right) \right\}.
```

The first-order condition is `c = v'(k)^(-1 / γ)`. 

```julia
using EconPDEs

# Model parameters.
const A = 0.5
const α = 0.3
const δ = 0.05
const ρ = 0.05
const γ = 2.0

# Declare the state grid: one state variable named `k`.
# Steady-state capital satisfies α A k̄^(α - 1) = ρ + δ.
k̄ = (α * A / (ρ + δ))^(1 / (1 - α))
stategrid = OrderedDict(:k => range(0.1 * k̄, 5.0 * k̄, length = 200))

# Declare the unknown value function: one function named `v`,
# with one initial value at each point of the `k` grid.
guess = OrderedDict(:v => [(A * k^α)^(1 - γ) / (1 - γ) / ρ for k in stategrid[:k]])

# `hjb` encodes the HJB equation at one grid point.
# `u.vk_up` and `u.vk_down` are one-sided finite-difference derivatives.
function hjb(state::NamedTuple, u::NamedTuple)
    k = state.k

    # Since capital can drift up or down, use the upwinded derivative.
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

# Solve the stationary HJB.
result = pdesolve(hjb, stategrid, guess)
value = result.zero[:v]
@assert result.residual_norm <= 1e-4
```

## Model Conventions

The PDE function (named `hjb` in example above) must have the form

```julia
f(state::NamedTuple, u::NamedTuple) -> NamedTuple
```

`state` gives the current grid point. `u` gives each unknown function and its finite-difference derivatives at that point. Think of `u` as the local finite-difference bundle: it contains the value and all stencil candidates the HJB may need at this grid node. The model then chooses the economically correct derivative, usually by the sign of a drift.

Field names concatenate the unknown name, the state name or names, and an optional direction suffix. If the unknown is `v` and the states are `k` and `l`, then `u` contains:

| Field | Meaning |
|---|---|
| `v` | local value of the unknown |
| `vk_up`, `vk_down` | forward and backward first derivatives in `k` |
| `vl_up`, `vl_down` | forward and backward first derivatives in `l` |
| `vkk`, `vll` | own second derivatives |
| `vkl` | central cross derivative |
| `vkl_up`, `vkl_down` | directional cross derivatives for positive and negative cross terms |

With several unknowns, the same naming rule applies to each one. For example, if the unknowns are `v` and `w` on states `k` and `l`, then `u` contains fields such as `v`, `vk_up`, `vkl`, `w`, `wk_up`, and `wkl`. The argument name `u` is only a convention; the fields are what matter. The PDE function must return one time derivative per unknown, e.g. `(; vt, wt)`.

## Saving Intermediate Outputs

Sometimes you compute useful objects inside the PDE function: a selected first derivative, an optimal policy, a drift, an interest rate, or a market price of risk. To store these objects on the grid, return two `NamedTuple`s from the PDE function. The first one is the PDE residual as before. The second one contains the extra outputs to save:

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
    return (; vt), (; c, μk)
end

result = pdesolve(hjb, stategrid, guess)
value = result.zero[:v]
consumption = result.optional[:c]
capital_drift = result.optional[:μk]
```

Each object in the second `NamedTuple` is evaluated at every state-grid point and stored as an array with the same shape as the grid. `result.zero` contains the solved unknowns. `result.optional` contains the saved outputs, and also includes the solved unknowns for convenience.

## Time-Dependent Problems

For a time-dependent problem, pass an increasing time grid as the fourth argument. If the HJB itself depends on time, define the PDE function with a third argument. For example, suppose productivity follows a deterministic path:

```julia
function hjb(state::NamedTuple, u::NamedTuple, t)
    k = state.k
    A_t = A * (1 + 0.1 * exp(-0.05 * t))

    c_up = max(u.vk_up, eps())^(-1 / γ)
    μ_up = A_t * k^α - δ * k - c_up

    c_down = max(u.vk_down, eps())^(-1 / γ)
    μ_down = A_t * k^α - δ * k - c_down

    if μ_up > 0
        c, vk, μk = c_up, u.vk_up, μ_up
    elseif μ_down < 0
        c, vk, μk = c_down, u.vk_down, μ_down
    else
        c, vk, μk = A_t * k^α - δ * k, u.vk_up, 0.0
    end

    vt = -(c^(1 - γ) / (1 - γ) + μk * vk - ρ * u.v)
    return (; vt)
end

τs = range(0, 100, length = 50)
result = pdesolve(hjb, stategrid, guess, τs)
```

`pdesolve` solves backward over this grid. In this case, `guess` is the terminal value at `τs[end]`, and `result.zero[i]` is the solution dictionary at time `τs[i]`; for example, `result.zero[i][:v]` is the value function at that time.

If the PDE function has only two arguments, `pdesolve` uses the same equation at each time on the grid.


## Boundary Conditions

Finite-difference schemes typicall need ghost-node values outside the grid to construct
derivatives at the boundaries. More precisely, ghost-node values are required to construct second derivatives when the state volatility is nonzero, or for first derivatives when the state drift points toward the boundary. By default, `EconPDEs.jl` uses reflecting boundaries: the ghost-node value
equals the boundary-node value, so the first derivative is zero beyond the boundary.

Use the `bc` keyword to impose a different boundary derivative:

```julia
pdesolve(f, grid, guess; bc = OrderedDict(:vx => (0.0, 1.0)))
```



Borrowing constraints and other endogenous state constraints are usually part of the PDE. Put assets, wealth, or the constrained object on the state grid. Then compute the policy inside the PDE function, form the state drift, and choose the upwind derivative from the sign of that drift.

For example, in a consumption-saving problem with asset state `a`, compute consumption from the marginal value of assets, form the asset drift `μa`, and then use `va_up` or `va_down` depending on the sign of `μa`. At the borrowing limit, impose feasibility directly: prevent the drift from moving below the lower asset bound, or impose the boundary derivative implied by the constraint.

See [WangWangYang.jl](examples/ConsumptionProblem/WangWangYang.jl), [AchdouHanLasryLionsMoll_Diffusion.jl](examples/ConsumptionProblem/AchdouHanLasryLionsMoll_Diffusion.jl), and [BoltonChenWang.jl](examples/InvestmentProblem/BoltonChenWang.jl) for examples.

## Optimal Stopping

Optimal stopping problems are supported through HJB variational inequalities. For a payoff `S(x)`, the HJBVI is:

```math
\min\left\{
    \rho v(x) - f(x) - \mu(x) v'(x) - \frac{1}{2}\sigma^2(x) v''(x),
    v(x) - S(x)
\right\} = 0.
```

Pass the exercise payoff as a lower bound with `y̲`; use `ȳ` for upper bounds in minimization problems. This formulation implies the usual value-matching and smooth-pasting conditions. See [Leland.jl](examples/OptimalStoppingTime/Leland.jl) for a working example and Ben Moll's [notes on stopping time problems](https://benjaminmoll.com/codes/) for background.

Internally, these bounded problems are mixed complementarity problems and use `NLsolve.mcpsolve`.

## Diagnostics and Solver Backend

For ordinary PDE solves, `EconPDEs.jl` constructs a sparse colored finite-difference Jacobian and passes it to `NonlinearSolve.jl`. The default nonlinear method is Newton's method:

```julia
pdesolve(f, grid, guess; method = :newton)
```

The trust-region solver is also available:

```julia
pdesolve(f, grid, guess; method = :trust_region)
```

As an opt-in diagnostic, `pdesolve` can check the sparse residual Jacobian for monotonicity violations:

```julia
pdesolve(f, grid, guess; check_monotonicity = true)
```

This checks the effective neighbor weights in the fully assembled equation, after endogenous policies, risk adjustments, and nonlinear terms are included. Upwinding on the raw state drift is only a shortcut that is guaranteed to work for a linear drift term. Under the package's `vt = -(...)` convention, a positive same-variable spatial off-diagonal entry usually means that an `_up` or `_down` stencil was chosen with the wrong sign. For transformed equations or endogenous-control problems, treat the warning as a debugging diagnostic rather than a proof that the economic upwind rule is wrong. Adjust `monotonicity_tol` (default `1e-6`) and `monotonicity_max_warnings` to tune it.

## Examples

The [examples folder](https://github.com/matthieugomez/EconPDEs.jl/tree/master/examples) works through a wide range of models, grouped by type. Together they exercise the full range of what `pdesolve` handles — one to several state variables, single and coupled value functions, stationary and time-dependent HJBs, and reflecting, state-constraint, degenerate, and free (optimal-stopping) boundaries.

**Asset pricing** ([examples/AssetPricing](examples/AssetPricing))
- *Habit*: Campbell–Cochrane (1999), Wachter (2005)
- *Long-run risk*: Bansal–Yaron (2004) — expected growth + stochastic variance
- *Disaster risk*: Wachter (2013)
- *Arbitrage with holding costs*: Tuckman–Vila (1992) — finite-horizon, time-dependent
- *Heterogeneous agents*: He–Krishnamurthy (2013), Brunnermeier–Sannikov (2013), Gârleanu–Panageas (2015), Di Tella (2017), Haddad — endogenous volatility, occasionally-binding constraints, up to several coupled value functions

**Consumption-saving with a borrowing constraint** ([examples/ConsumptionProblem](examples/ConsumptionProblem))
- Achdou–Han–Lasry–Lions–Moll (2018) — diffusion, two-asset, and two-income-state variants
- Wang–Wang–Yang (2016)

**Investment with a financing constraint** ([examples/InvestmentProblem](examples/InvestmentProblem))
- Bolton–Chen–Wang (2009)

**Optimal stopping** ([examples/OptimalStoppingTime](examples/OptimalStoppingTime)) — free-boundary problems as HJB variational inequalities (value-matching + smooth-pasting)
- Leland (1994) — endogenous default

**Growth** ([examples/GrowthModel](examples/GrowthModel))
- Deterministic neoclassical (Ramsey) growth, on a grid centered on its closed-form steady state

## Citation
If you use EconPDEs.jl in your work, please cite it. You can use the "Cite this repository" button on the [GitHub page](https://github.com/matthieugomez/EconPDEs.jl) (generated from [`CITATION.cff`](CITATION.cff)), or the following BibTeX entry:

```bibtex
@software{EconPDEs.jl,
  author = {Gomez, Matthieu},
  title  = {{EconPDEs.jl}: {Solving PDEs in economics}},
  url    = {https://github.com/matthieugomez/EconPDEs.jl},
}
```
