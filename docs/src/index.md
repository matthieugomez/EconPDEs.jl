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

In these models, the PDE is usually tied to an optimal policy. The policy depends on the
value function, and the value function depends on the policy. That feedback is the part
that general PDE interfaces often do not make convenient.

| Economic HJBs need | What this means in practice |
| --- | --- |
| Infinitesimal generators | The operator describes how the state moves locally: drift times the first derivative plus volatility squared times the second derivative. The companion package [InfinitesimalGenerators.jl](https://github.com/matthieugomez/InfinitesimalGenerators.jl) works with these operators directly; see [InfinitesimalGenerators](infinitesimal_generators.md). |
| Policies inside the equation | Drift, diffusion, controls, and nonlinear terms can depend on the current guess for the value function. You write these terms directly, in the same place as the policy rule. |
| Upwinding by drift sign | The first derivative must be taken from the right side. If the drift changes during the solve, the stencil changes too. See [Writing the PDE function](pde_function.md#Upwinding). |
| Economic boundaries | Borrowing limits, reflecting states, and degenerate diffusions are not just Dirichlet or periodic boundary conditions. See [Boundary conditions](boundary_conditions.md). |
| Stationary nonlinear solves | Sparse finite-difference Jacobians and continuation make fine grids practical. See [Solver and troubleshooting](solver.md). |

`EconPDEs.jl` handles those details so the model code can stay close to the economics:
write the local equation, and the package builds the derivatives, boundary treatment,
upwinding, and nonlinear solve.
