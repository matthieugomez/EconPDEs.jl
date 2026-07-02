# Why EconPDEs

Julia already has good PDE packages. `EconPDEs.jl` has a narrower goal: solve the
continuous-time HJBs that come up in economic models.

In these models, the PDE is usually tied to an optimal policy. The policy depends on the
value function, and the value function depends on the policy. That feedback is the part
that general PDE interfaces often do not make convenient.

## The Difference

| Economic HJBs need | What this means in practice |
| --- | --- |
| Infinitesimal generators | The operator describes how the state moves locally: drift times the first derivative plus volatility squared times the second derivative. The companion package [InfinitesimalGenerators.jl](https://github.com/matthieugomez/InfinitesimalGenerators.jl) works with these operators directly; see [InfinitesimalGenerators](infinitesimal_generators.md). |
| Policies inside the equation | Drift, diffusion, controls, and nonlinear terms can depend on the current guess for the value function. You write these terms directly, in the same place as the policy rule. |
| Upwinding by drift sign | The first derivative must be taken from the right side. If the drift changes during the solve, the stencil changes too. See [Getting started](getting_started.md#Upwinding). |
| Economic boundaries | Borrowing limits, reflecting states, and degenerate diffusions are not just Dirichlet or periodic boundary conditions. See [Solving models](solving.md#Boundary-conditions). |
| Stationary nonlinear solves | Sparse finite-difference Jacobians and continuation make fine grids practical. See [Solving models](solving.md#Solver-and-troubleshooting). |

`EconPDEs.jl` handles those details so the model code can stay close to the economics:
write the local equation, and the package builds the derivatives, boundary treatment,
upwinding, and nonlinear solve.
