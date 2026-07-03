[![Build status](https://github.com/matthieugomez/EconPDEs.jl/workflows/CI/badge.svg)](https://github.com/matthieugomez/EconPDEs.jl/actions)
[![Coverage Status](http://codecov.io/github/matthieugomez/EconPDEs.jl/coverage.svg?branch=main)](http://codecov.io/github/matthieugomez/EconPDEs.jl/?branch=main)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://www.matthieugomez.com/EconPDEs.jl/stable/)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://www.matthieugomez.com/EconPDEs.jl/dev/)

# EconPDEs.jl

`EconPDEs.jl` solves nonlinear ODEs and PDEs that arise in economic models, especially Hamilton-Jacobi-Bellman equations.

You write the equation at one grid point. The package builds the finite-difference derivatives, exposes upwind stencils, assembles sparse Jacobians, and solves the system by pseudo-transient Newton iteration.

Use it for stationary or time-dependent HJBs on one or more state variables, with one or more value functions.

`EconPDEs.jl` is robust, fast, and solves virtually all textbook continuous-time models; see the Examples for a list.

## Installation

The package is registered in the Julia `General` registry:

```julia
] add EconPDEs
```

Current versions require Julia 1.10 or later.

## Quickstart

Here is a small linear equation on one state, with a mean-reverting drift:

```math
\rho v(x) = x + \mu(x) v_x(x),
\qquad
\mu(x) = -\kappa (x - 0.5).
```

`pdesolve` finds the stationary solution by embedding the equation in a time-dependent
problem, solved backward in time from an initial guess until the time derivative is zero.
```math
\rho v(x, t) = x + \mu(x) v_x(x, t) + v_t(x, t),
```

The main function `pdesolve` takes three objects:

- a local equation `f(state, u)` encoding the PDE;
- a `NamedTuple` encoding the state grid, such as `(; k = range(...))`;
- a `NamedTuple` encoding the initial guess, such as `(; v = [...])`.

Here is how to solve this equation using this package:
```julia
using EconPDEs

const rho = 0.05
const kappa = 0.2

grid = (; x = range(0.0, 1.0, length = 50))
guess = (; v = zeros(length(grid.x)))

function equation(state, u)
    mu = -kappa * (state.x - 0.5)
    # The grid and guess names determine the fields: state.x, u.v, u.vx_up, and u.vx_down.
    vx = mu >= 0 ? u.vx_up : u.vx_down
    # Return the time derivative implied by rho*v = x + mu*v_x + v_t.
    vt = -(state.x + mu * vx - rho * u.v)
    return (; vt)
end

result = pdesolve(equation, grid, guess; verbose = false)

value = result.zero[:v]
@assert result.residual_norm <= 1e-4
```

The drift is positive below `0.5` and negative above `0.5`, so the equation chooses the
forward or backward derivative at each grid point.

The same interface extends to multiple states or unknowns. Grid entries can be ranges or
vectors, and different dimensions can use different vector containers:

```julia
grid = (; a = range(0.0, 100.0, length = 500),
          y = collect(range(0.5, 1.5, length = 10)))
guess = (; v = zeros(length(grid.a), length(grid.y)))
```

For two unknowns, return `(; vt, wt)`. Return a second `NamedTuple` to save policies,
drifts, prices, or other objects on the grid; they are available as `result.saved`.
Full HJB examples, including policy choice and endogenous upwinding, are in the
documentation.

## Documentation

The manual gives the details:

- [Getting started](https://www.matthieugomez.com/EconPDEs.jl/dev/getting_started/)
- [Writing the PDE function](https://www.matthieugomez.com/EconPDEs.jl/dev/pde_function/)
- [Boundary conditions](https://www.matthieugomez.com/EconPDEs.jl/dev/boundary_conditions/)
- [Time-dependent problems](https://www.matthieugomez.com/EconPDEs.jl/dev/time_dependent/)
- [Solver and troubleshooting](https://www.matthieugomez.com/EconPDEs.jl/dev/solver/)
- [Why EconPDEs](https://www.matthieugomez.com/EconPDEs.jl/dev/#Why-EconPDEs)
- [InfinitesimalGenerators](https://www.matthieugomez.com/EconPDEs.jl/dev/infinitesimal_generators/) — using the companion package [InfinitesimalGenerators.jl](https://github.com/matthieugomez/InfinitesimalGenerators.jl) for stationary distributions and expectations after solving

The documentation also includes a gallery of runnable examples:

- consumption-saving models, including borrowing constraints, income risk, and risky assets;
- asset-pricing models, including habit, long-run risk, rare disasters, and intermediary finance;
- corporate-finance models, including optimal default and investment with financing constraints.

The source scripts live in [`examples/`](examples).

## Citation

If you use EconPDEs.jl in your work, please cite it. You can use the "Cite this repository" button on the [GitHub page](https://github.com/matthieugomez/EconPDEs.jl), or the following BibTeX entry:

```bibtex
@software{EconPDEs.jl,
  author = {Gomez, Matthieu},
  title  = {{EconPDEs.jl}: {Solving PDEs in economics}},
  url    = {https://github.com/matthieugomez/EconPDEs.jl},
}
```
