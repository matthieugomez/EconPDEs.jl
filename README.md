[![Build status](https://github.com/matthieugomez/EconPDEs.jl/workflows/CI/badge.svg)](https://github.com/matthieugomez/EconPDEs.jl/actions)
[![Coverage Status](http://codecov.io/github/matthieugomez/EconPDEs.jl/coverage.svg?branch=main)](http://codecov.io/github/matthieugomez/EconPDEs.jl/?branch=main)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://www.matthieugomez.com/EconPDEs.jl/stable/)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://www.matthieugomez.com/EconPDEs.jl/dev/)

# EconPDEs.jl

`EconPDEs.jl` solves nonlinear ODEs and PDEs that arise in economic models, especially Hamilton-Jacobi-Bellman equations.

You specify three things: (i) the state space, (ii) an initial guess for the solution, and (iii) a function that returns the PDE residual at each grid point. The package handles the rest — it builds the finite-difference derivatives, applies upwind stencils, assembles the sparse Jacobian, and solves the system by pseudo-transient Newton iteration.

Use it for stationary or time-dependent HJBs on one or more state variables, with one or more value functions.

`EconPDEs.jl` is robust, fast, and solves virtually all textbook continuous-time models; see the [Examples](https://www.matthieugomez.com/EconPDEs.jl/dev/examples_overview/) for a list.

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



Here is how to solve this equation using this package:
```julia
using EconPDEs

const rho = 0.05
const kappa = 0.2

# 1. Define a `NamedTuple` representing the state space as a grid
grid = (; x = range(0.0, 1.0, length = 50))

# 2. Define a `NamedTuple` representing an initial guess for the solution
guess = (; v = zeros(length(grid.x)))

# 3. Define a function encoding the PDE at each grid point
function pde(state, u)
    # The grid names determine the fields of state: state.x.
    mu = -kappa * (state.x - 0.5)
    # The guess names determine the fields of u: u.v, u.vx_up, and u.vx_down.
    # Choose first-derivative direction based on the sign of the drift (upwinding)
    vx = mu >= 0 ? u.vx_up : u.vx_down
    # Return the time derivative implied by rho*v = x + mu*v_x + v_t.
    vt = -(state.x + mu * vx - rho * u.v)
    return (; vt)
end

# 4. Pass these three objects to `pdesolve`
result = pdesolve(pde, grid, guess; verbose = false)
# EconPDEResult
#   solution:      v (50)
#   residual_norm: 1.22e-10
#   converged:     true (tolerance 1.49e-08)
```

The same interface extends to multiple states...

```julia
grid = (; x = range(0.0, 1.0, length = 50),
          y = collect(range(0.5, 1.5, length = 10)))
guess = (; v = zeros(length(grid.x), length(grid.y)))

function pde(state, u)
    (; x, y) = state
    (; v, vx_up, vx_down, vy_up, vy_down, vxx, vyy, vxy_up, vxy_down) = u
    ...
    return (; vt)
end
```

...as well as multiple value functions (i.e., coupled PDES). 
```julia
grid = (; x = range(0.0, 1.0, length = 50))
guess = (; v = zeros(length(grid.x)), w = zeros(length(grid.x)))

function pde(state, u)
    (; x) = state
    (; v, vx_up, vx_down, vxx, w, wx_up, wx_down, wxx) = u
    ...
    return (; vt, wt)
end
```


## Documentation

Read the [documentation](https://www.matthieugomez.com/EconPDEs.jl/dev/) for further details. 
The documentation also includes a gallery of runnable examples (see [`examples/`](examples) for the original scripts): 

- consumption-saving models, including borrowing constraints, income risk, and risky assets;
- asset-pricing models, including habit, long-run risk, rare disasters, and intermediary finance;
- corporate-finance models, including optimal default and investment with financing constraints.

## Citation

Please cite EconPDEs.jl if it is useful in your work. You can use the "Cite this repository" button on the [GitHub page](https://github.com/matthieugomez/EconPDEs.jl), or the following BibTeX entry:

```bibtex
@software{EconPDEs.jl,
  author = {Gomez, Matthieu},
  title  = {{EconPDEs.jl}: {Solving PDEs in economics}},
  url    = {https://github.com/matthieugomez/EconPDEs.jl},
}
```
