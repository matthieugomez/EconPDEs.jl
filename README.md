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

The main function is `pdesolve`. A stationary problem needs a PDE function, a state grid, and an initial guess. Here is an example for the PDE `∂v_t = x^2 - 0.04 * v + μx * vx + 0.05 * vxx`. 

```julia
using EconPDEs

# Declare the state grid: one state variable named `x`.
grid = OrderedDict(:x => range(-1.0, 1.0, length = 100))

# Declare the unknown value function: one function named `v`,
# with one initial value at each point of the `x` grid.
guess = OrderedDict(:v => ones(length(grid[:x])))

# `hjb` encodes the PDE at one point in the state space
# `state.x` is one value from `grid[:x]`.
# `y.v` is the local value of `v`; `y.vx_up`, `y.vx_down`, and `y.vxx`
# are finite-difference derivatives inferred from the names `v` and `x`.
function hjb(state::NamedTuple, y::NamedTuple)
    x = state.x
    v, vx_up, vx_down, vxx = y.v, y.vx_up, y.vx_down, y.vxx

    μx = -x
    vx = μx >= 0 ? vx_up : vx_down

    vt = -(x^2 + μx * vx + 0.05 * vxx - 0.04 * v)
    return (; vt)
end

# Solve the stationary HJB.
result = pdesolve(hjb, grid, guess; verbose = false)
value = result.zero[:v]
residual_norm = result.residual_norm
```

## Model Conventions

The PDE function has the form

```julia
f(state::NamedTuple, y::NamedTuple) -> NamedTuple
```

`state` gives the current grid point. `y` gives each unknown function and its finite-difference derivatives at that point. If the unknown is `v` and the state is `x`, then `y` contains:

- `v`: value of the function.
- `vx_up` and `vx_down`: forward and backward first derivatives.
- `vxx`: second derivative.

With several states, cross derivatives are also available, such as `vxy`, `vxy_up`, and `vxy_down`. The PDE function must return one time derivative for each unknown; for `v`, return `vt`.

For higher-dimensional systems, the same rule applies:

```julia
agrid = range(0.0, 10.0, length = 100)
zgrid = range(-0.5, 0.5, length = 50)

# Two state variables, named `a` and `z`.
grid2 = OrderedDict(:a => agrid, :z => zgrid)

# Two unknown value functions, named `v` and `w`.
# Each array has dimensions `(length(agrid), length(zgrid))`,
# following the order of the state variables in `grid2`.
guess2 = OrderedDict(
    :v => ones(length(agrid), length(zgrid)),
    :w => ones(length(agrid), length(zgrid)),
)

# A PDE function called with this `grid2` and `guess2` receives
# `state.a` and `state.z`.
#
# For `v`, `y` contains `y.v`, first derivatives such as
# `y.va_up`, `y.va_down`, `y.vz_up`, and `y.vz_down`,
# second derivatives such as `y.vaa` and `y.vzz`,
# and cross derivatives such as `y.vaz`, `y.vaz_up`, and `y.vaz_down`.
#
# For `w`, `y` contains the analogous fields: `y.w`, `y.wa_up`,
# `y.wz_down`, `y.waa`, `y.wzz`, `y.waz`, and so on.
# The PDE function must return one time derivative per unknown: `(; vt, wt)`.
```

## Saving Intermediate Outputs

Sometimes you compute useful objects inside the PDE function: a selected first derivative, an optimal policy, a drift, an interest rate, or a market price of risk. To store these objects on the grid, return two `NamedTuple`s from the PDE function. The first one is the PDE residual as before. The second one contains the extra outputs to save:

```julia
function hjb(state::NamedTuple, y::NamedTuple)
    x = state.x
    v, vx_up, vx_down, vxx = y.v, y.vx_up, y.vx_down, y.vxx

    μx = -x
    vx = μx >= 0 ? vx_up : vx_down

    vt = -(x^2 + μx * vx + 0.05 * vxx - 0.04 * v)
    return (; vt), (; vx, μx)
end

result = pdesolve(hjb, grid, guess; verbose = false)
value = result.zero[:v]
vx = result.optional[:vx]
drift = result.optional[:μx]
```

Each object in the second `NamedTuple` is evaluated at every state-grid point and stored as an array with the same shape as the grid. `result.zero` contains the solved unknowns. `result.optional` contains the saved outputs, and also includes the solved unknowns for convenience.

## Time-Dependent Problems

For a time-dependent problem, pass an increasing time grid as the fourth argument:

```julia
τs = range(0, 100, length = 50)
result = pdesolve(hjb, grid, guess, τs; verbose = false)
```

`pdesolve` solves backward over this grid. In this case, `guess` is the terminal value at `τs[end]`, and `result.zero[i]` is the solution dictionary at time `τs[i]`; for example, `result.zero[i][:v]` is the value function at that time.

If the HJB itself depends on time, define the PDE function with a third argument:

```julia
function hjb(state::NamedTuple, y::NamedTuple, t)
    x = state.x
    v, vx_up, vx_down, vxx = y.v, y.vx_up, y.vx_down, y.vxx

    μx = -x
    vx = μx >= 0 ? vx_up : vx_down

    flow_payoff = exp(-0.01 * t) * x^2
    vt = -(flow_payoff + μx * vx + 0.05 * vxx - 0.04 * v)
    return (; vt)
end

result = pdesolve(hjb, grid, guess, τs; verbose = false)
```

If the PDE function has only two arguments, `pdesolve` uses the same equation at each time on the grid.

## Endogenous Choices and Borrowing Constraints

Borrowing constraints and other endogenous state constraints are usually part of the PDE. Put assets, wealth, or the constrained object on the state grid. Then compute the policy inside the PDE function, form the state drift, and choose the upwind derivative from the sign of that drift.

For example, in a consumption-saving problem with asset state `a`, compute consumption from the marginal value of assets, form the asset drift `μa`, and then use `va_up` or `va_down` depending on the sign of `μa`. At the borrowing limit, impose feasibility directly: prevent the drift from moving below the lower asset bound, or impose the boundary derivative implied by the constraint.

The `y̲` and `ȳ` keywords have a different role. They impose bounds on the unknown function itself in a variational inequality, as in optimal stopping. They are not the right tool for an ordinary borrowing constraint on a state variable or policy choice.

See [WangWangYang.jl](examples/ConsumptionProblem/WangWangYang.jl), [AchdouHanLasryLionsMoll_Diffusion.jl](examples/ConsumptionProblem/AchdouHanLasryLionsMoll_Diffusion.jl), and [BoltonChenWang.jl](examples/InvestmentProblem/BoltonChenWang.jl) for examples.

## Boundary Conditions

Finite-difference schemes need ghost-node values outside the grid to construct derivatives at the boundary. By default, `EconPDEs.jl` uses reflecting boundaries: the ghost-node value equals the boundary-node value, so the first derivative is zero beyond the boundary.

Use the `bc` keyword to impose a different boundary derivative:

```julia
pdesolve(f, grid, guess; bc = OrderedDict(:vx => (0.0, 1.0)))
```

## Optimal Stopping

Optimal stopping problems are supported through HJB variational inequalities. For a payoff `S(x)`, the HJBVI is:

<img src="img/hjbvi.png">

Pass the exercise payoff as a lower bound with `y̲`; use `ȳ` for upper bounds in minimization problems. This formulation implies the usual value-matching and smooth-pasting conditions. See [Leland.jl](examples/OptimalStoppingTime/Leland.jl) for a working example and Ben Moll's [notes on stopping time problems](https://benjaminmoll.com/codes/) for background.

Internally, these bounded problems are mixed complementarity problems and use `NLsolve.mcpsolve`. Ordinary nonlinear PDE solves use `NonlinearSolve.jl`.

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
