# Boundary conditions

Finite-difference schemes need *ghost-node* values just outside the grid to construct some
derivatives at the boundary: second derivatives whenever the state volatility is nonzero
there, and first derivatives when the drift points toward the boundary. This page covers
the default, the `bc` keyword, and the two boundary types that are really part of the
economics — state constraints and free (optimal-stopping) boundaries.

## The default: reflecting boundaries

By default, the ghost-node value equals the boundary-node value, so the first derivative is
zero beyond the boundary — a reflecting boundary. This is the natural choice when the
boundary is far enough out that little probability mass reaches it (e.g. the top of an
asset grid), and it is exact when the state genuinely reflects.

Note that `pdesolve` still applies the *full PDE* at the boundary node, using the ghost
value only to build the stencil. It does not replace the equation at the boundary with the
boundary condition — on coarse grids this is noticeably more accurate (see
[Why EconPDEs](design.md)).

## The `bc` keyword

To impose a different boundary derivative, pass `bc`, a `NamedTuple` mapping
`Symbol(unknown, state)` to a `(lower, upper)` tuple of outward first derivatives:

```julia
pdesolve(f, grid, guess; bc = (; va = (0.0, 1.0)))
```

Scalars apply uniformly along the boundary; arrays matching the boundary slice are also
accepted (e.g. in 2D, a vector giving one derivative per grid point of the other state).
Entries may be given for any subset of (unknown, state) pairs; omitted pairs keep the
reflecting default.

A degenerate boundary — where the diffusion vanishes and the drift points inward, like a
square-root process at zero — needs no condition at all: the PDE at the boundary node only
uses one-sided derivatives there, which is exactly what the upwind scheme provides.

## State constraints (borrowing limits)

Borrowing constraints and other endogenous state constraints are usually *part of the PDE*,
not a `bc` entry. Put the constrained quantity on the state grid, compute the policy inside
the PDE function, and let the upwind logic impose feasibility at the constraint: at the
lower bound, if the household would like to dissave past the limit, the drift is pinned at
zero and consumption equals income — case 3 of the endogenous-drift pattern in
[Upwinding](upwinding.md). Alternatively, some models imply a known boundary derivative at
the constraint (e.g. from the FOC at zero saving), which can be imposed directly through
`bc`.

See the [consumption–saving examples](examples/consumption_saving/consumption_saving_diffusion_income.md)
and [Wang–Wang–Yang](examples/consumption_saving/wang_wang_yang.md).

## Optimal stopping (free boundaries)

Optimal stopping problems — default, exercise, exit — are HJB *variational inequalities*.
For a stopping payoff ``S(x)``:

```math
\min\left\{
    \rho v(x) - f(x) - \mu(x) v'(x) - \tfrac{1}{2}\sigma^2(x) v''(x),\;
    v(x) - S(x)
\right\} = 0.
```

Pass the stopping payoff as a *lower bound* on the unknown with the `y̲` keyword (`ȳ` for
upper bounds in minimization problems); at the REPL, type `y\underbar<tab>` and
`y\bar<tab>`:

```julia
result = pdesolve(f, grid, guess; y̲ = vec(payoff_on_grid))
```

The free boundary then comes out of the solution — the region where the bound binds is the
stopping region, and value matching and smooth pasting hold automatically at its edge.
Internally these bounded problems are solved as mixed complementarity problems with
`NLsolve.mcpsolve`.

See [Leland](examples/corporate_finance/leland.md) for a worked example (optimal default),
and Ben Moll's [notes on stopping-time problems](https://benjaminmoll.com/codes/) for
background.
