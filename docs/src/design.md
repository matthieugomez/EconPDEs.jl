# Why EconPDEs

Julia has excellent general PDE software. Why a package specific to economics? Because a
continuous-time HJB has structure that general PDE frameworks are not organized around —
and that structure, not the nonlinear solve, is the hard part.

## What an economics HJB demands

1. **Non-conservative, non-divergence form.** The spatial operator is the infinitesimal
   generator ``\mathcal{A}v = b(x)\cdot\nabla v + \tfrac12 \mathrm{tr}(\Sigma(x)\nabla^2 v)``,
   not a flux divergence ``\nabla \cdot F``. Finite-volume and discontinuous-Galerkin
   machinery is organized entirely around fluxes; the generator does not factor through
   one.
2. **Gradient-nonlinear, value-divided terms.** Recursive utility and market prices of risk
   inject terms like ``\sigma^2 (\nabla p)^2 / p`` — nonlinear in the gradient *and*
   divided by the value. This fits no slot in a conservation-law API (it is neither a flux
   nor a node-local reaction).
3. **State-dependent, endogenous coefficients.** Drift and diffusion vary over the domain
   and, in controlled problems, depend on the unknown itself.
4. **Monotone upwinding by drift sign.** Convergence to the viscosity solution
   (Barles–Souganidis) requires a monotone scheme, with the first-derivative stencil chosen
   by the sign of a drift that — for endogenous drift — changes every Newton iteration.
   See [Upwinding](upwinding.md).
5. **Reflecting, state-constraint, and degenerate boundaries.** Not Dirichlet or periodic.
   The economically important boundaries are state constraints (the control keeps the state
   in the domain) and degenerate edges (diffusion vanishes, e.g. a square-root process at
   zero, where no condition should be imposed). See
   [Boundary conditions](boundary_conditions.md).
6. **Endogenous state constraints tied to the control.** A borrowing limit is enforced by
   overriding the marginal value so the implied drift is zero — a model-specific operation
   entangled with the optimal policy, not a declarative boundary condition.
7. **Stationary problems solved by Newton with continuation**, using a sparse
   finite-difference Jacobian with matrix coloring so fine grids stay cheap. See
   [Solver](solver.md).

`EconPDEs.jl` is a thin layer that encodes exactly these pieces, on top of standard
nonlinear solvers. You write the local equation; the scaffolding is the package.

## A case study against general packages

As a controlled experiment, we re-solved one EconPDEs example — the Bansal–Yaron (2004)
long-run-risk model, a stationary nonlinear 2-state HJB — with two general Julia PDE
packages, checking each against the converged EconPDEs reference (residual `1.6e-9`) on the
same 30×30 grid.

- **Trixi.jl** (high-order DG/FV for conservation laws): the HJB cannot be expressed
  through its equation interface — fluxes take constant coefficients and never see the
  position; sources see no gradients; the ``(\nabla p)^2/p`` term fits nowhere. The port
  ended up bypassing Trixi's PDE stack entirely, hand-rolling a residual around the one
  usable piece, an upwind summation-by-parts operator that Trixi re-exports from
  `SummationByPartsOperators.jl`. Best relative error: `6e-3` — on a *different* discrete
  solution, since high-order SBP is not monotone.
- **MethodOfLines.jl** (symbolic method of lines): structurally the natural fit — the full
  nonlinear HJB expresses symbolically, and its `UpwindScheme()` upwinds correctly. It
  converged, but to a solution about `1.25e-2` away on the coarse grid. The cause is its
  boundary paradigm: MOL *replaces the PDE at the boundary node* with a first-order Neumann
  equation, while EconPDEs applies the full HJB at the boundary using a reflecting ghost
  node. Both are consistent, but the ``O(h)`` boundary discrepancy is advected across the
  whole domain by mean-reverting drift, producing a ~1% level offset that only refines away
  at finer grids. Multi-Neumann corners were also left unassigned, and the default
  steady-state route (`DynamicSS` with a dense Jacobian) did not scale past small grids.

The point is fit, not solver speed: the value-adding part of solving these models is the
monotone upwinding, the boundary treatment, the sparse generator Jacobian, and the
continuation that makes a flat guess converge. General packages each provide one useful
piece and leave the econ-specific pieces to the user.

*(Caveats: one model, one grid, ports written by us rather than by those packages'
maintainers. The structural arguments follow from the packages' source and stand
regardless; the accuracy numbers are from honest but non-expert ports.)*
