# EconPDEs.jl

`EconPDEs.jl` solves nonlinear ODEs and PDEs that arise in economic models, especially
Hamilton–Jacobi–Bellman equations. You write the local equation; the package supplies the
finite-difference derivatives, upwinding, sparse Jacobians, and pseudo-transient Newton
iteration.

Repository: [github.com/matthieugomez/EconPDEs.jl](https://github.com/matthieugomez/EconPDEs.jl)

## Installation

The package is registered in the Julia `General` registry:

```julia
] add EconPDEs
```

Current versions require Julia 1.10 or later.

## Where to start

Start with [Getting started](getting_started.md), which solves a first model end to end.
[Writing the PDE function](pde_function.md) is the reference for the package's naming
conventions and the upwinding patterns used in most models.
[Boundary conditions](boundary_conditions.md),
[Time-dependent problems](time_dependent.md), and
[Solver and troubleshooting](solver.md) cover what usually matters after the first model
works. The [`pdesolve`](api.md) reference lists every solver option.

The **Examples** in the sidebar each solve a model, plot the solution, and explain what it
means economically. The [overview](examples_overview.md) lists them all, flagging what each
demonstrates beyond the baseline workflow. Start with
[Neoclassical growth](examples/neoclassical_growth.md). The
[InfinitesimalGenerators](infinitesimal_generators.md) page explains when to combine
EconPDEs with `InfinitesimalGenerators.jl` for stationary distributions or lower-level
finite-difference residuals.

## Why EconPDEs

Julia already has good PDE packages. `EconPDEs.jl` has a narrower goal: solve the
continuous-time HJBs that come up in economic models.

These equations have a structure that is easy to lose in a generic PDE interface. The
state drift and volatility are often chosen by the agent, so the differential operator
depends on the current guess for the value function. At the same time, the equation is
strongly forward-looking: today's value depends on future values, future policies, and the
way the state distribution moves under those policies.

`EconPDEs.jl` keeps that feedback loop explicit. You write the local economic equation,
including the policy rules and the implied drift or volatility. The package then builds the
finite-difference derivatives, boundary treatment, sparse Jacobian, and nonlinear solve
needed to turn that equation into a stationary or time-dependent solution.
