# EconPDEs.jl

`EconPDEs.jl` solves nonlinear ODEs and PDEs that arise in economic models, especially
Hamilton–Jacobi–Bellman equations. You write the local equation; the package supplies the
finite-difference derivatives, upwinding, sparse Jacobians, and pseudo-transient Newton
iteration.

Repository: [github.com/matthieugomez/EconPDEs.jl](https://github.com/matthieugomez/EconPDEs.jl)

## Installation

```julia
] add EconPDEs
```

## Getting started

Start with [Getting started](getting_started.md), which solves a first model end to end.
It also covers the package's naming conventions and the upwinding patterns used in most
models. [Solving models](solving.md) collects boundary conditions, time-dependent
problems, and solver troubleshooting. [Why EconPDEs](design.md) explains what this package
does that general PDE software does not. The [`pdesolve`](api.md) reference lists every
solver option.

The **Examples** in the sidebar each solve a model, plot the solution, and explain what it
means economically. They are grouped by field — consumption–saving, asset pricing,
corporate finance — and ordered roughly from the simplest to the most involved. Start with
[Neoclassical growth](examples/neoclassical_growth.md). The
[InfinitesimalGenerators](infinitesimal_generators.md) page explains when to combine
EconPDEs with `InfinitesimalGenerators.jl` for stationary distributions or lower-level
finite-difference residuals.
