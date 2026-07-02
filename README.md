[![Build status](https://github.com/matthieugomez/EconPDEs.jl/workflows/CI/badge.svg)](https://github.com/matthieugomez/EconPDEs.jl/actions)
[![Coverage Status](http://codecov.io/github/matthieugomez/EconPDEs.jl/coverage.svg?branch=main)](http://codecov.io/github/matthieugomez/EconPDEs.jl/?branch=main)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://www.matthieugomez.com/EconPDEs.jl/dev/)

# EconPDEs.jl

`EconPDEs.jl` solves nonlinear ODEs and PDEs that arise in economic models, especially Hamilton-Jacobi-Bellman equations.

You write the equation at one grid point. The package builds the finite-difference derivatives, exposes upwind stencils, assembles sparse Jacobians, and solves the system by pseudo-transient Newton iteration.

Use it for stationary or time-dependent HJBs on one or more state variables, with one or more value functions.

## Installation

The package is registered in the Julia `General` registry:

```julia
] add EconPDEs
```

Current versions require Julia 1.10 or later.

## Quickstart

The main function is `pdesolve`. It takes three objects:

- a local equation `f(state, u)`;
- a named state grid, such as `(; k = range(...))`;
- a named initial guess, such as `(; v = [...])`.

Here is a small linear equation on one state, with a mean-reverting drift:

```math
\rho v(x) = x + \mu(x) v_x(x),
\qquad
\mu(x) = -\kappa (x - 0.5).
```

`pdesolve` finds the stationary solution by embedding the equation in a time-dependent
problem and driving `v_t` to zero. For this example, the stationary residual is

```math
x + \mu(x) v_x(x) - \rho v(x),
```

so the PDE function returns its negative:

```julia
using EconPDEs

const rho = 0.05
const kappa = 0.2

grid = (; x = range(0.0, 1.0, length = 50))
guess = (; v = zeros(length(grid.x)))

function equation(state, u)
    mu = -kappa * (state.x - 0.5)
    vx = mu >= 0 ? u.vx_up : u.vx_down
    vt = -(state.x + mu * vx - rho * u.v)
    return (; vt)
end

result = pdesolve(equation, grid, guess; verbose = false)

value = result.zero[:v]
@assert result.residual_norm <= 1e-4
```

The drift is positive below `0.5` and negative above `0.5`, so the equation chooses the
forward or backward derivative at each grid point. Full HJB examples, including policy
choice and endogenous upwinding, are in the documentation.

## Interface

The state grid and the initial guess are `NamedTuple`s keyed by state names and unknown names. Grid entries can be ranges or vectors, and different dimensions can use different vector containers:

```julia
grid = (; a = range(0.0, 100.0, length = 500),
          y = collect(range(0.5, 1.5, length = 10)))
guess = (; v = zeros(length(grid.a), length(grid.y)))
```

Inside your equation, `state` is the current grid point and `u` contains the value and derivative stencils. If the unknown is `v` and the state is `k`, then `u.v` is the local value and `u.vk_up`, `u.vk_down`, and `u.vkk` are finite-difference derivatives.

Return one time derivative per unknown, named by appending `t`: `(; vt)` for unknown `v`, or `(; vt, wt)` for unknowns `v` and `w`. Return a second `NamedTuple` to save policies, drifts, prices, or other objects on the grid.

## Documentation

The manual gives the details:

- [Getting started](https://www.matthieugomez.com/EconPDEs.jl/dev/getting_started/)
- [Solving models](https://www.matthieugomez.com/EconPDEs.jl/dev/solving/)
- [Why EconPDEs](https://www.matthieugomez.com/EconPDEs.jl/dev/design/)
- [InfinitesimalGenerators](https://www.matthieugomez.com/EconPDEs.jl/dev/infinitesimal_generators/)

The documentation also includes a gallery of runnable examples:

- consumption-saving models, including borrowing constraints, income risk, and risky assets;
- asset-pricing models, including habit, long-run risk, rare disasters, and intermediary finance;
- corporate-finance models, including optimal default and investment with financing constraints.

The source scripts live in [`examples/`](examples).

## Numerical Details

For the finite-difference discretization and the pseudo-transient Newton scheme behind
`pdesolve`, see [Solving models](https://www.matthieugomez.com/EconPDEs.jl/dev/solving/#Solver-and-troubleshooting).

## Citation

If you use EconPDEs.jl in your work, please cite it. You can use the "Cite this repository" button on the [GitHub page](https://github.com/matthieugomez/EconPDEs.jl), or the following BibTeX entry:

```bibtex
@software{EconPDEs.jl,
  author = {Gomez, Matthieu},
  title  = {{EconPDEs.jl}: {Solving PDEs in economics}},
  url    = {https://github.com/matthieugomez/EconPDEs.jl},
}
```
