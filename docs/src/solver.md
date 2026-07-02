# Solver and troubleshooting

## How `pdesolve` finds the solution

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

## Solver options

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
| `autodiff` | `:forward` | `:forward`/`:finite` (forward differences) or `:central`; with 1–3 states the Jacobian always uses colored sparse finite differences |
| `verbose` | `true` | print the iteration table |
| `y̲`, `ȳ` | `-Inf`, `Inf` | bounds for variational inequalities (see [Boundary conditions](boundary_conditions.md)) |
| `is_algebraic` | all `false` | mark equations with no time derivative (static/algebraic conditions) |

## When a solve fails

Roughly in order of likelihood:

1. **It diverges immediately (residual grows from iteration 1).** Check the sign
   convention: the PDE function must return `vt = -(RHS - ρv)`, not the raw residual. See
   [Model conventions](conventions.md).
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

## The monotonicity diagnostic

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
