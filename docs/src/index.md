# EconPDEs.jl

`EconPDEs.jl` solves nonlinear ODEs and PDEs that arise in economic models, especially
Hamilton–Jacobi–Bellman equations. You write the local equation; the package supplies the
finite-difference derivatives, upwinding, sparse Jacobians, and pseudo-transient Newton
iteration.

## Installation

```julia
] add EconPDEs
```

## Getting started

The [README](https://github.com/matthieugomez/EconPDEs.jl) documents the modeling
conventions, boundary conditions, and optimal-stopping support; the [`pdesolve`](api.md)
reference lists every solver option.

The **Examples** in the sidebar each solve a model, plot the solution, and explain what it
means economically. They are grouped into macro and finance models and ordered roughly from the simplest to the
most involved — start with [Neoclassical growth](examples/neoclassical_growth.md).
