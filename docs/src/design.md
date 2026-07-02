# Why EconPDEs

Julia already has good PDE packages. `EconPDEs.jl` has a narrower goal: solve the
continuous-time HJBs that come up in economic models.

In these models, the PDE is usually tied to an optimal policy. The policy depends on the
value function, and the value function depends on the policy. That feedback is the part
that general PDE interfaces often do not make convenient.

## The Difference

| Economic HJBs need | What this means in practice |
| --- | --- |
| Infinitesimal generators | The operator looks like ``b(x) v_x + \frac12 \sigma(x)^2 v_{xx}``, not a conservation-law flux. |
| Endogenous policies | Drift, diffusion, and other coefficients can depend on the current guess for the value function. |
| Upwinding by drift sign | The first derivative must be taken from the right side. If the drift changes during the solve, the stencil changes too. See [Getting started](getting_started.md#Upwinding). |
| Economic boundaries | Borrowing limits, reflecting states, and degenerate diffusions are not just Dirichlet or periodic boundary conditions. See [Solving models](solving.md#Boundary-conditions). |
| Direct nonlinear terms | Terms like ``\sigma^2 p_x^2 / p`` are easier to write directly than to force into a flux/source split. |
| Stationary nonlinear solves | Sparse finite-difference Jacobians and continuation make fine grids practical. See [Solving models](solving.md#Solver-and-troubleshooting). |

`EconPDEs.jl` handles those details so the model code can stay close to the economics:
write the local equation, and the package builds the derivatives, boundary treatment,
upwinding, and nonlinear solve.
