# Solver and troubleshooting

## How `pdesolve` finds the solution

`pdesolve` first turns the local PDE into a finite system. Start from a
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

`pdesolve` finds this zero in three steps:

1. Build the finite system with upwind finite differences.
2. Solve backward implicit false-transient steps. In a linear valuation problem this is a
   linear M-matrix solve; when policies or equilibrium objects depend on the new value, the
   package solves the nonlinear equation with a nonlinear solver, Newton by default.
3. Adapt `Δ`: grow it when the stationary residual falls, and shrink it when a step fails or
   the residual rises.

The basic equation in step 2 is a backward implicit time step:

```math
\frac{y_{n+1} - y_n}{\Delta} = F(y_{n+1}).
```

When this equation is nonlinear, `pdesolve` can use different nonlinear solvers. The
default is Newton's method. With algebraic equations (`is_algebraic = true`), the
corresponding rows skip the time-step term and are solved as static equations.

The size of `Δ` controls how aggressive the step is. Small `Δ` makes the equation close to
``y_{n+1} = y_n``, so the step is cautious but progress is slow. Large `Δ` makes the
time-step term negligible; in nonlinear problems, the step approaches a direct Newton solve
of ``F(y) = 0``. The outer continuation loop in step 3 starts cautiously and adapts:

1. Starting from the guess, solve one implicit step of size `Δ`.
2. If the inner solve converges, accept the step; grow `Δ` when the residual falls, and shrink it when the residual rises.
3. If the inner solve fails, reject the step, shrink `Δ` by 10, and retry. Regrowth is
   then capped below the failed `Δ` until the residual has clearly improved, so a `Δ`
   that just failed is not immediately retried.
4. Stop when the residual norm falls below `maxdist`.

This adaptive scheme is *pseudo-transient continuation*, imported from computational fluid
dynamics; formal convergence conditions are given in
[Kelley and Keyes (1998)](https://doi.org/10.1137/S0036142996304796).

With `verbose = true` (the default), `pdesolve` prints what it is solving, one line per
outer iteration, and a convergence summary (this is the habit-formation model of the
example gallery):

```
Solving for 1 unknown (p) on a 998-point grid
Iter TimeStep Residual
---- -------- --------
   1 1.00e+00 2.51e-02
   2 1.12e+01 1.21e-02
   3 2.33e+03 8.01e-05
   4 3.50e+08 3.45e-12
Converged after 4 time steps (4ms)
```

`TimeStep` is the current `Δ` and `Residual` the norm of the stationary residual at the
accepted step. A step is kept even when its residual rises; the next `Δ` then shrinks
(it grows when the residual falls). A `✗` in place of the residual means the inner
nonlinear solve **failed** at that `Δ`: no step is taken, so there is no residual to
show, and the step is retried with `Δ/10` — the unusual causes, the model returned `NaN`
or a singular Jacobian, are flagged in parentheses:

```
   2 1.04e+01        ✗
```

Healthy solves show `Δ` growing geometrically and the residual collapsing in the last
few time steps.

With `verbose = false` a successful solve prints nothing, and convergence failures are
still reported with `@warn` (with a hint at the likely cause) — so `pdesolve` can run
inside a loop, e.g. an estimation, without flooding the log while failures still surface.
Either way, the returned result records whether the solve converged:

```
julia> result
EconPDEResult
  zero:          p (998)
  residual_norm: 3.45e-12
  converged:     true (tolerance 1.49e-08)
```

## Why the method works

**Upwinding makes the discretized problem an economy.** With first derivatives taken from
the side the drift points to, the discretized equation at a grid point reads

```math
\rho v_i = u_i + \lambda_i^+ (v_{i+1} - v_i) + \lambda_i^- (v_{i-1} - v_i),
\qquad \lambda_i^\pm = \frac{\mu_i^\pm}{\Delta x} + \frac{\sigma_i^2}{2\Delta x^2} \ge 0,
```

which is the valuation equation of a continuous-time Markov chain on the grid: the drift
becomes a directed jump intensity, the variance a symmetric one. A wrong-direction stencil
makes some ``\lambda`` negative — the discrete model assigns negative probabilities and the
false transient explodes; that is what
[the monotonicity diagnostic](#The-monotonicity-diagnostic) detects. In the numerical-PDE
literature this nonnegativity is called *monotonicity*, and monotone consistent stable
schemes converge to the (viscosity) solution of the PDE as the grid is refined
([Barles and Souganidis, 1991](https://doi.org/10.3233/ASY-1991-4305)); the chain
construction itself is
[Kushner and Dupuis's](https://doi.org/10.1007/978-1-4613-0007-6) Markov chain
approximation method.

**Case 1: Linear valuation problem.** Take the law of motion, payoff, and discount rate as
given while solving for the value function. Equivalently, write
``F(y) = u + Ay - \rho y`` with ``A`` the transition-intensity matrix of the discretized
state process. The implicit step is then a linear system, not a nonlinear problem. It has
an exact probabilistic meaning:

```math
y_{n+1}(x) = \mathbb{E}\left[\int_0^T e^{-\rho t}\, u(X_t)\, dt \;+\; e^{-\rho T} y_n(X_T)
\;\middle|\; X_0 = x\right], \qquad T \sim \text{Exponential}(1/\Delta).
```

One step values the model *exactly* out to a random horizon with mean `Δ`, then hands over
to the previous iterate as continuation value. The error therefore contracts by the average
discount factor ``\mathbb{E}[e^{-\rho T}] = 1/(1+\rho\Delta)`` per step — Blackwell's
value-function-iteration argument, with ``1/(1+\rho\Delta)`` as the discount factor. The
step is stable for *any* `Δ` because its matrix is an M-matrix. An explicit step
``y_n + \Delta F(y_n)`` is instead one round of value function iteration with period `Δ`,
and the period is capped by the grid (``\Delta \lesssim \Delta x^2/\sigma^2``, the CFL
condition) — on fine grids, that iteration crawls.

**Solving nonlinear steps.** When ``F`` depends on the new candidate value through policies,
prices, or equilibrium objects, the backward implicit step is a nonlinear equation.
`pdesolve` can use different nonlinear solvers; the default is Newton's method, with a
trust-region variant available for harder systems. Let the residual of one implicit step be

```math
R_\Delta(W;V_n) = \frac{W-V_n}{\Delta} - F(W).
```

Given an inner iterate ``W_k``, Newton computes the residual and Jacobian, solves

```math
D R_\Delta(W_k)\, \delta_k = -R_\Delta(W_k),
```

and updates ``W_{k+1}=W_k+\delta_k``, with safeguards when needed. The inner loop stops
when the residual of the implicit step is small. The outer loop then decides whether to
accept the step and how to change `Δ`.

**Case 2: Pure HJBs.** A pure HJB adds an optimizer to the linear valuation problem. Once a
policy ``c`` is fixed, the equation is again a linear valuation equation; the nonlinearity
comes from the fact that the optimal policy itself depends on the candidate value. Let
``A(c)`` denote the finite-difference operator obtained after fixing policy ``c``. Given the
previous outer iterate ``V_n``, the fully implicit step asks for a new value ``W`` satisfying

```math
\frac{W - V_n}{\Delta} =
\max_c \left\{u(c) + A(c)W - \rho W\right\}.
```

This is a nonlinear problem in ``W`` because a new candidate value changes marginal values
and therefore the optimal policy ``c``.

For this pure HJB, Newton has a familiar interpretation. At Newton iterate ``W_k``,
compute the policy ``c_k`` that is optimal for ``W_k``. By the envelope theorem, the terms
involving the derivative of ``c_k`` with respect to ``W_k`` drop from the Jacobian. The
Newton step therefore holds the policy fixed and solves the linear policy-evaluation
problem

```math
\frac{W_{k+1} - V_n}{\Delta} =
u(c_k) + A(c_k) W_{k+1} - \rho W_{k+1}.
```

The policy is then recomputed at ``W_{k+1}``, and the process repeats until the nonlinear
implicit equation is solved. This is the policy-iteration logic inside the implicit step.

The comparison to [Achdou et al. (2022)](https://doi.org/10.1093/restud/rdab002) is now
direct. Their update is semi-implicit for the nonlinear HJB: it computes the policy from
the old value ``V_n``, freezes it, and solves one linear system. Equivalently, it is the
first frozen-policy Newton step on the fully implicit nonlinear equation, started from
``W_0 = V_n``. `pdesolve` iterates the inner Newton solve to convergence, so the final
policy is optimal for the new value function, not just for the old guess.

In discrete-time dynamic-programming terms, Achdou's semi-implicit step is one
modified-policy-iteration step: choose the policy from ``V_n``, then partially or fully
evaluate it. Small `Δ` looks like value function iteration; intermediate `Δ` is modified
policy iteration; and ``\Delta=\infty`` is Howard policy iteration.

This interpretation is useful for judging the comparison. Achdou's method is a strong
semi-implicit policy-iteration baseline for clean pure HJBs. `pdesolve` can be more
forgiving with poor guesses or aggressive `Δ` because it reoptimizes inside the implicit
step and accepts a large step only after the policy is consistent with the new value.
Still, the larger robustness gain appears in Case 3, where freezing equilibrium objects is
no longer policy iteration.

**Case 3: Equilibrium systems.** In asset-pricing and general-equilibrium models, the
residual includes prices, risk premia, volatilities, market-clearing equations, and
sometimes algebraic constraints. These objects are not maximized-away controls, so the
envelope theorem does not remove their derivatives from the Jacobian. Freezing them gives a
semi-implicit fixed-point iteration, not Newton's method. `pdesolve` instead treats the
whole discretized equilibrium as a sparse nonlinear system and solves that system with
Newton by default, with a trust-region variant available for harder cases.

**The time step controls the regime.** The role of `Δ` is separate from the three cases
above. For finite `Δ`, each Newton or policy-evaluation step is anchored to the previous
outer iterate ``V_n``. As ``\Delta \to \infty`` the anchor disappears. In a linear valuation
problem, this is just the stationary linear system. In a pure stationary HJB, the Newton
loop becomes Howard policy iteration: improve the policy, evaluate it, repeat. In an
equilibrium system, the same limit is a direct Newton solve of the full stationary residual.
Growing `Δ` adaptively therefore moves the solver from stabilized iteration toward a fast,
local stationary solve.

The same ideas travel under different names across fields:

| ingredient here | numerical analysis / CFD | discrete-time dynamic programming |
|---|---|---|
| upwind differences | monotone scheme | Markov chain approximation |
| explicit step | forward Euler (CFL-capped) | value function iteration |
| implicit step | backward Euler (stable for all `Δ`) | valuation to a random horizon of mean `Δ` |
| `Δ = Inf` | stationary nonlinear solve (Newton by default) | Howard policy iteration for HJBs; direct Newton for equilibrium systems |
| adaptive `Δ` | pseudo-transient continuation | value iteration morphing into policy iteration or Newton |

In the nonlinear cases, the extra Newton work is cheap enough to be useful because the
Jacobian is sparse. For one-, two-, and three-state problems, `pdesolve` builds the
sparsity pattern from the local stencil and computes the Jacobian by colored finite
differences. Recomputing this Jacobian lets the solver take much larger accepted `Δ` steps,
so the benefit is speed as well as robustness.

A fuller writeup with references is the method note in the repository:
[`examples/details.pdf`](https://github.com/matthieugomez/EconPDEs.jl/blob/main/examples/details.pdf).

## Solver options

All of these are keyword arguments of `pdesolve` (forwarded to
[`finiteschemesolve`](api.md)):

| Keyword | Default | Meaning |
|---|---|---|
| `Δ` | `1.0` | initial pseudo-time step; `Inf` = stationary nonlinear solve, no continuation |
| `scale` | `10.0` | growth factor of `Δ` after a successful step |
| `minΔ`, `maxΔ` | `1e-9`, `Inf` | bounds on `Δ`; the solve aborts when `Δ < minΔ` |
| `iterations` | `100` | maximum outer iterations |
| `maxdist` | `sqrt(eps())` | convergence tolerance on the residual norm |
| `inner_iterations` | `10` | nonlinear-solver iterations per implicit time step |
| `innerdist` | `sqrt(eps())` | tolerance of the inner nonlinear solve |
| `method` | `:newton` | `:newton` or `:trust_region` |
| `autodiff` | `:forward` | `:forward`/`:finite` (forward differences) or `:central`; with 1-3 states the Jacobian always uses colored sparse finite differences |
| `verbose` | `true` | print the problem summary, the iteration table, and a convergence summary; failures warn even when `false` |
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
4. **The inner nonlinear solve keeps failing (many `rejected` lines).** Try `method =
   :trust_region`, a smaller initial `Δ` (e.g. `1e-2` for stiff systems with several
   coupled unknowns), or loosen `inner_iterations`.
5. **It converges, but to something economically wrong.** Check `residual_norm` is actually
   small; refine the grid (the solution should be stable under refinement); and inspect the
   solved policies at the boundaries — a misplaced boundary condition shows up there first.
6. **Slow.** The cost is dominated by the sparse Jacobian. Keep the guess good (fewer outer
   iterations), start from a coarse grid and interpolate the solution onto the fine one,
   and pass `Δ = Inf` when the guess is already close (a stationary solve with no continuation).

## The monotonicity diagnostic

As an opt-in diagnostic, `pdesolve` can check the assembled residual Jacobian for
monotonicity violations:

```julia
pdesolve(pde, grid, guess; check_monotonicity = true)
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
where the solve is a stationary nonlinear solve and either sign is legitimate.
