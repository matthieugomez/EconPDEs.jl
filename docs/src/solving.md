# Solving models

This page collects the pieces that usually matter after the first model works: boundary
conditions, time-dependent problems, and solver troubleshooting.

## Boundary conditions

Finite-difference schemes need *ghost-node* values just outside the grid to construct some
derivatives at the boundary: second derivatives whenever the state volatility is nonzero
there, and first derivatives when the drift points toward the boundary. This section
covers the default, the `bc` keyword, and the two boundary types that are really part of
the economics — state constraints and free (optimal-stopping) boundaries.

### The default: reflecting boundaries

By default, the ghost-node value equals the boundary-node value, so the first derivative is
zero beyond the boundary — a reflecting boundary. This is the natural choice when the
boundary is far enough out that little probability mass reaches it (e.g. the top of an
asset grid), and it is exact when the state genuinely reflects.

Note that `pdesolve` still applies the *full PDE* at the boundary node, using the ghost
value only to build the stencil. It does not replace the equation at the boundary with the
boundary condition — on coarse grids this is noticeably more accurate (see
[Why EconPDEs](design.md)).

### The `bc` keyword

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

### State constraints (borrowing limits)

Borrowing constraints and other endogenous state constraints are usually *part of the PDE*,
not a `bc` entry. Put the constrained quantity on the state grid, compute the policy inside
the PDE function, and let the upwind logic impose feasibility at the constraint: at the
lower bound, if the household would like to dissave past the limit, the drift is pinned at
zero and consumption equals income — case 3 of the endogenous-drift pattern in
[Getting started](getting_started.md#Pattern-2:-endogenous-drift-(consumption-saving)).
Alternatively, some models imply a known boundary derivative at the constraint (e.g. from
the FOC at zero saving), which can be imposed directly through `bc`.

See the [consumption-saving examples](examples/consumption_saving/consumption_saving_diffusion_income.md)
and [Wang-Wang-Yang](examples/consumption_saving/wang_wang_yang.md).

### Optimal stopping (free boundaries)

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

## Time-dependent problems

For a time-dependent problem — a transition path, a finite horizon, a deterministic change
in parameters — pass an increasing time grid as the fourth positional argument:

```julia
τs = range(0, 100, length = 50)
result = pdesolve(f, stategrid, guess, τs)
```

`pdesolve` solves *backward* over this grid: `guess` is the terminal value at `τs[end]`,
and each step solves one implicit time step of ``\partial_t v = \text{vt}`` — the same
equation, and the same sign convention, as the stationary false transient (see
[Getting started](getting_started.md#The-return-value:-one-time-derivative-per-unknown)).
The result is one solution per time point: `result.zero[i][:v]` is the value function at
time `τs[i]`, and `result.residual_norm[i]` the corresponding residual.

### Time-varying equations

If the equation itself depends on time, define the PDE function with a third argument. For
example, productivity following a deterministic path:

```julia
function hjb(state::NamedTuple, u::NamedTuple, t)
    A_t = A * (1 + 0.1 * exp(-0.05 * t))
    # ... same as the stationary function, with A_t in place of A ...
    return (; vt)
end
```

If the PDE function has only two arguments, the same equation is used at every time on the
grid. The third argument must be left untyped or typed loosely (e.g. `t` or `t::Number`) —
`pdesolve` detects it with `hasmethod(f, Tuple{NamedTuple, NamedTuple, Number})`.

### Typical uses

- **Finite-horizon problems**: the terminal value is known (e.g. zero, or a bequest value)
  and the object of interest is the whole path. See
  [Tuckman-Vila](examples/asset_pricing/tuckman_vila.md).
- **Transition dynamics**: solve the stationary problem first, then use it as the terminal
  value of a long time grid to trace convergence from an initial condition.
- **A robust path to a hard stationary solution**: when the stationary solve fails from a
  poor guess, integrating the time-dependent equation over a long horizon is a natural
  continuation method — it follows the physically meaningful transient.

## Solver and troubleshooting

### How `pdesolve` finds the solution

Discretizing the PDE on the grid turns it into a large nonlinear system: one equation
(the returned time derivative) per unknown per grid point, coupled to its neighbors through
the finite-difference stencils. `pdesolve` drives this residual to zero with **Newton's
method plus pseudo-transient continuation**:

1. Starting from the guess, take an *implicit time step* of size `Δ` of the false
   transient ``\partial_t v = \text{vt}`` — itself a Newton solve, with a sparse Jacobian
   computed by colored finite differences.
2. If the step succeeds and reduces the residual, grow `Δ` (by
   `scale × old/new residual`) and repeat; if the inner Newton solve fails, shrink `Δ`
   by 10 and retry.
3. Stop when the residual norm falls below `maxdist` — or report failure when the
   iteration count exceeds `iterations` or `Δ` collapses below `minΔ`.

Small `Δ` makes each step easy (the identity term dominates the Jacobian) but slow; large
`Δ` approaches a pure Newton step. The scheme automatically transitions from cautious to
aggressive as it approaches the solution, which is what lets crude guesses (even flat
functions) converge. For the finite-difference discretization and the implicit scheme in
detail, see the
[numerical appendix (PDF)](https://github.com/matthieugomez/EconPDEs.jl/blob/main/examples/details.pdf).

With `verbose = true` (the default), each outer iteration prints

```
Iter   TimeStep   Residual
---- ---------- ----------
   1 1.0000e+00 4.4708e-03
   2 1.4390e+01 6.4159e-05
```

`TimeStep` is the current `Δ` and `Residual` the norm of the stationary residual. A `NaN`
residual line means the inner Newton solve failed at that `Δ` and the step is being
retried with `Δ/10`. Healthy solves show `Δ` growing geometrically and the residual
collapsing in the last few iterations.

### Solver options

All of these are keyword arguments of `pdesolve` (forwarded to
[`finiteschemesolve`](api.md)):

| Keyword | Default | Meaning |
|---|---|---|
| `Δ` | `1.0` | initial pseudo-time step; `Inf` = single Newton solve, no continuation |
| `scale` | `10.0` | growth factor of `Δ` after a successful step |
| `minΔ`, `maxΔ` | `1e-9`, `Inf` | bounds on `Δ`; the solve aborts when `Δ < minΔ` |
| `iterations` | `100` | maximum outer iterations |
| `maxdist` | `sqrt(eps())` | convergence tolerance on the residual norm |
| `inner_iterations` | `10` | Newton iterations per implicit time step |
| `innerdist` | `sqrt(eps())` | tolerance of the inner Newton solve |
| `method` | `:newton` | `:newton` or `:trust_region` |
| `autodiff` | `:forward` | `:forward`/`:finite` (forward differences) or `:central`; with 1-3 states the Jacobian always uses colored sparse finite differences |
| `verbose` | `true` | print the iteration table |
| `y̲`, `ȳ` | `-Inf`, `Inf` | bounds for variational inequalities (see [Boundary conditions](#Boundary-conditions)) |
| `is_algebraic` | all `false` | mark equations with no time derivative (static/algebraic conditions) |

### When a solve fails

Roughly in order of likelihood:

1. **It diverges immediately (residual grows from iteration 1).** Check the sign
   convention: the PDE function must return `vt = -(RHS - ρv)`, not the raw residual. See
   [Getting started](getting_started.md#The-return-value:-one-time-derivative-per-unknown).
2. **`NaN` at the initial value.** `pdesolve` throws
   `ArgumentError: G! returns NaN with the initial value`: the PDE function hit an invalid
   operation at the *guess* (negative marginal value raised to a power, division by zero).
   Guard the FOC against nonpositive derivatives — e.g. cap the implied policy, as the
   examples do with `cmax` — and check the guess itself is in the function's domain.
3. **The residual stalls and `Δ` keeps shrinking until `minΔ`.** The classic symptom of a
   non-monotone scheme: an upwind direction chosen with the wrong sign, or a missing
   directional cross derivative. Run with `check_monotonicity = true` (below). If the
   scheme is clean, improve the initial guess — a closed-form limit of the model (no risk,
   log utility, unconstrained) is usually enough.
4. **The inner Newton solve keeps failing (many `NaN` lines).** Try `method =
   :trust_region`, a smaller initial `Δ` (e.g. `1e-2` for stiff systems with several
   coupled unknowns), or loosen `inner_iterations`.
5. **It converges, but to something economically wrong.** Check `residual_norm` is actually
   small; refine the grid (the solution should be stable under refinement); and inspect the
   solved policies at the boundaries — a misplaced boundary condition shows up there first.
6. **Slow.** The cost is dominated by the sparse Jacobian. Keep the guess good (fewer outer
   iterations), start from a coarse grid and interpolate the solution onto the fine one,
   and pass `Δ = Inf` when the guess is already close (a single Newton solve).

### The monotonicity diagnostic

As an opt-in diagnostic, `pdesolve` can check the assembled residual Jacobian for
monotonicity violations:

```julia
pdesolve(f, grid, guess; check_monotonicity = true)
```

This checks the effective neighbor weights in the fully assembled equation, after
endogenous policies, risk adjustments, and nonlinear terms are included. Under the
`vt = -(...)` convention, a positive same-variable spatial off-diagonal entry usually means
an `_up`/`_down` stencil was chosen with the wrong sign; the warning names the variable,
the grid point, the offending neighbor, and the stencil that most likely caused it. For
transformed equations or endogenous-control problems, treat it as a debugging aid rather
than proof that the economic upwind rule is wrong. Tune with `monotonicity_tol` (default
`1e-6`) and `monotonicity_max_warnings` (default `5`).
