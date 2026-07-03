# Solver and troubleshooting

## How `pdesolve` finds the solution

`pdesolve` first turns the local PDE into a finite nonlinear system. Start from a
time-dependent equation such as

```math
\partial_t V(x, t) = f(x, V, V_x, V_{xx}).
```

After choosing a state grid, the values of `V` on the grid are stacked into a vector `y`.
Finite differences approximate `V_x` and `V_{xx}`. Evaluating the user's equation at every
grid point gives a vector field

```math
\dot y = F(y).
```

For a stationary problem, the goal is simply

```math
F(y) = 0.
```

`pdesolve` finds this zero with **Newton's method plus pseudo-transient continuation**. The
basic step is a backward implicit time step:

```math
\frac{y_{n+1} - y_n}{\Delta} = F(y_{n+1}).
```

This is itself a nonlinear system, solved by Newton's method. With algebraic equations
(`is_algebraic = true`), the corresponding rows skip the time-step term and are solved as
static equations.

The size of `Δ` controls how aggressive the step is. Small `Δ` makes the equation close to
``y_{n+1} = y_n``, so Newton is stable but progress is slow. Large `Δ` makes the time-step
term negligible, so the step approaches a direct Newton solve of ``F(y) = 0``. The solver
starts cautiously and adapts:

1. Starting from the guess, solve one implicit step of size `Δ`.
2. If the step succeeds and reduces the residual, grow `Δ` and continue.
3. If the inner Newton solve fails, shrink `Δ` by 10 and retry.
4. Stop when the residual norm falls below `maxdist`.

This adaptive scheme is *pseudo-transient continuation*, imported from computational fluid
dynamics; formal convergence conditions are given in
[Kelley and Keyes (1998)](https://doi.org/10.1137/S0036142996304796).

The Jacobian is sparse because each finite-difference equation only depends on nearby grid
points. For one-, two-, and three-state problems, `pdesolve` builds this sparsity pattern
and computes the Jacobian by colored finite differences. This is why fine grids remain
practical: each Newton step uses the local stencil structure rather than a dense Jacobian.

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
| `autodiff` | `:forward` | `:forward`/`:finite` (forward differences) or `:central`; with 1-3 states the Jacobian always uses colored sparse finite differences |
| `verbose` | `true` | print the iteration table |
| `lower_bound`, `upper_bound` | `-Inf`, `Inf` | bounds for variational inequalities (see [Optimal stopping](boundary_conditions.md#Optimal-stopping-(free-boundaries))) |
| `is_algebraic` | all `false` | mark equations with no time derivative (static/algebraic conditions) |

## When a solve fails

Roughly in order of likelihood:

1. **It diverges immediately (residual grows from iteration 1).** Check the sign
   convention: the PDE function must return `vt = -(RHS - ρv)`, not the raw residual. See
   [Writing the PDE function](pde_function.md#The-return-value:-one-time-derivative-per-unknown).
   `pdesolve` detects the typical signature automatically and warns (see
   [the monotonicity diagnostic](#The-monotonicity-diagnostic) below).
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

A related check runs *by default*, without any option: under the `vt = -(RHS - ρv)`
convention the diagonal of the residual Jacobian is positive (the discount rate plus the
outflow rates of the discretized generator), so a diagonal that is negative at **every**
grid point is the signature of a flipped sign convention — `RHS - ρv` returned instead of
`-(RHS - ρv)`. `pdesolve` then warns on the first iteration, since the false transient is
unstable and the solve would otherwise just diverge. The check is skipped when `Δ = Inf`,
where the solve is a single Newton step and either sign is legitimate.
