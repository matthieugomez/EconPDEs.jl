# Changelog

All notable changes to EconPDEs.jl are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/), and the project follows
[semantic versioning](https://semver.org/).

## [3.0.0] — 2026-07-19

### Breaking Changes

- Revert back to using NLsolve as the default nonlinear solver, selected with
  `alg = :newton` (the default) or `alg = :trust_region`. NonlinearSolve.jl is now
  optional and is no longer exported by EconPDEs; load it explicitly before passing one
  of its algorithm objects. Loading NonlinearSolve activates a package extension, so
  calls such as `alg = NonlinearSolve.TrustRegion()` remain supported.

### Improvements

- `Δ` may now be a `NamedTuple`/`OrderedDict` in `pdesolve` (or a vector in
  `finiteschemesolve`), giving each unknown its own pseudo-time step. The sign of an
  entry sets the direction of that unknown's false transient: positive marches its
  equation backward in time (the value-function convention), negative forward — so
  market-clearing equations written in their natural tâtonnement direction
  (`rt = r_implied - r`) can be relaxed toward equilibrium jointly with the value
  functions, e.g. `Δ = (; pA = Inf, r = -1.0)`. `±Inf` entries disable damping for an
  unknown. The entries fix a relative profile that the continuation scales by one common
  adaptive factor.
- Add a `max_residual_growth = 10.0` guard to the pseudo-transient continuation: an
  accepted implicit step whose stationary residual grows by more than this factor (a
  blowup, or a `NaN`, which counts as infinite growth) is now reverted and retried with
  a smaller time step instead of accepted. The default is deliberately loose — false
  transients may legitimately increase the residual, and rejecting every worsening step
  can deadlock the continuation — so well-behaved solves are unaffected.
- The sign-convention diagnostic now runs per unknown and understands relaxation
  directions: an unknown whose residual-Jacobian diagonal is opposite in sign to its
  pseudo-time step at every grid point is reported by name, with the suggestion to flip
  either its `Δ` entry or its returned time derivative.
- The Gârleanu–Panageas long-run-risk example now solves values and prices jointly with
  a signed per-unknown `Δ` (values backward on a fast clock, prices forward on a damped
  one), replacing the hand-rolled alternation between valuation solves and damped price
  updates. The solve has two stages — a `maxΔ`-capped continuation to a loose tolerance
  followed by one undamped Newton polish — and converges across grid sizes (tested from
  10×5×5 to 100×5×5 and 50×9×9) where both the uncapped continuation and the old
  alternation's damped price iteration degrade.

### Documentation

- Standardize section headings across the example gallery.

## [2.1.0] — 2026-07-10

### Improvements

- Add a three-state Gârleanu–Panageas long-run-risk asset-pricing example. The model
  keeps the interest rate and three prices of risk as algebraic unknown functions and
  alternates between valuation and market-clearing price updates.
- Solver summaries now report the names and sizes of the state-grid dimensions, such as
  `grid x×y×z = 2×3×4`.
- Add `warn = true` to `pdesolve` and `finiteschemesolve`, separating convergence-failure
  warnings from the existing `verbose = true` progress output. Pass `warn = false` when
  failures are expected and the caller checks convergence explicitly.

### Documentation

- Expand the README tutorial to separate the state grid, initial guess, and pointwise PDE
  residual, and add multidimensional and coupled-PDE interface examples.
- Document when to substitute equilibrium prices into an HJB and when to keep them as
  algebraic unknown functions.
- Clarify the numerical role of `Δ`: state the conditions under which `Δ = Inf` gives
  globally convergent Howard policy iteration for a discrete HJB, explain why finite `Δ`
  can help Newton, and qualify the convergence result for related linearized
  pseudo-transient methods.
- Refresh the example explanations and comments throughout the gallery.

### Compatibility

- Fix deprecated result-alias tests on Julia nightly by matching their deprecation
  messages explicitly.

## [2.0.0] — 2026-07-07

### Breaking Changes

#### Migration Map

| 1.x API | 2.0 API | Applies to |
|---|---|---|
| `method = :newton` | `alg = NonlinearSolve.NewtonRaphson()` | `pdesolve`, `finiteschemesolve` |
| `method = :trust_region` | `alg = NonlinearSolve.TrustRegion()` | `pdesolve`, `finiteschemesolve` |
| `iterations` | `maxiters` | `pdesolve`, `finiteschemesolve` |
| `inner_iterations` | `inner_maxiters` | `pdesolve`, `finiteschemesolve` |
| `maxdist` | `abstol` | `pdesolve`, `finiteschemesolve` |
| `innerdist` | `inner_abstol` | `pdesolve`, `finiteschemesolve` |
| `y̲`, `ȳ` | `lower_bound`, `upper_bound` | `pdesolve`, `finiteschemesolve` |
| `J0` | `jac_prototype` | `finiteschemesolve` |
| `result.zero` | `result.solution` (`result.zero` remains as a deprecated alias) | `pdesolve` output |
| `result.optional` | `result.saved` (`result.optional` remains as a deprecated alias) | `pdesolve` output |

Additionally, `result.solution` is now a `NamedTuple` rather than an `OrderedDict`.

#### `pdesolve`

- Keyword change: `method` was removed; use `alg = NonlinearSolve.NewtonRaphson()` or
  another compatible NonlinearSolve algorithm object such as
  `NonlinearSolve.TrustRegion()`.
- Keyword change: `iterations` was renamed to `maxiters`.
- Keyword change: `inner_iterations` was renamed to `inner_maxiters`.
- Keyword change: `maxdist` was renamed to `abstol`.
- Keyword change: `innerdist` was renamed to `inner_abstol` for time-dependent inner
  solves.
- Keyword change: `y̲` was renamed to `lower_bound`.
- Keyword change: `ȳ` was renamed to `upper_bound`.
- Removed keyword: `autodiff`; `pdesolve` uses its sparse stencil Jacobian and colored
  finite differences directly.
- Removed keyword: `reformulation`; bounded solves always use the minmax
  mixed-complementarity residual.
- Removed keyword: `autoscale`.
- Output (`EconPDEResult`): `result.zero` was renamed to `result.solution`.
- Output (`EconPDEResult`): `result.optional` was renamed to `result.saved`.
- Output (`EconPDEResult`): `result.zero` and `result.optional` remain as deprecated
  aliases.
- Output (`EconPDEResult`): `result.solution` and `result.saved` are `NamedTuple`s of arrays.
- Output (`EconPDEResult`): when the PDE saves nothing, `result.saved === NamedTuple()`.
- Output (`EconPDEResult`): time-dependent solutions and saved objects use a trailing time
  dimension, so `result.solution.v[.., i]` is the solution at `τs[i]`.
- Output (`EconPDEResult`): `converged` is a computed property using the stored solver
  tolerance. Tuple destructuring remains `solution, residual_norm, saved = result`.
- Output (`EconPDEResult`): the legacy three-argument constructor was removed.

#### `finiteschemesolve`

- Keyword change: `method` was removed; use `alg = NonlinearSolve.NewtonRaphson()` or
  another compatible NonlinearSolve algorithm object such as
  `NonlinearSolve.TrustRegion()`.
- Keyword change: `J0` was renamed to `jac_prototype`.
- Keyword change: `iterations` was renamed to `maxiters`.
- Keyword change: `inner_iterations` was renamed to `inner_maxiters`.
- Keyword change: `maxdist` was renamed to `abstol`.
- Keyword change: `innerdist` was renamed to `inner_abstol`.
- Keyword change: `y̲` was renamed to `lower_bound`.
- Keyword change: `ȳ` was renamed to `upper_bound`.
- New keyword: `jac` can provide an in-place Jacobian for the stationary residual.
- New keyword: `colorvec` can provide the column coloring of `jac_prototype`; otherwise a
  coloring is computed from the supplied sparsity pattern.
- Removed keyword: `autodiff`; choose AD/linear-solver behavior through the
  NonlinearSolve algorithm object where applicable.
- Removed keyword: `reformulation`; bounded solves always use the minmax
  mixed-complementarity residual.
- Removed keyword: `autoscale`.

#### `EconPDEs`

- `EconPDEs` no longer re-exports `OrderedDict` and no longer depends on
  OrderedCollections.

### Improvements

- `EconPDEs` now exports the `NonlinearSolve` module, so users can write
  `using EconPDEs` and then pass `NonlinearSolve.NewtonRaphson()` or any other compatible
  NonlinearSolve algorithm object.
- Sparse coloring now works for any number of state dimensions.
- `pdesolve` time grids must be strictly increasing, time-dependent PDE methods are
  detected by applicability, stationary-only pseudo-transient step controls are rejected
  when a time grid is supplied, bounds may be `NamedTuple`s keyed like `guess`, and
  convergence failures warn even with `verbose = false`.
- Sparse Jacobian coloring is computed in closed form from the stencil structure
  (coordinate mod 3 in each state dimension, crossed with the unknown index) instead of by
  greedy graph coloring of the assembled pattern.
- Residual evaluation uses a compile-once barrier when a sparsity pattern is available,
  reducing repeated model-solve latency.
- A `PrecompileTools` workload now solves a tiny model during package precompilation,
  caching the model-independent solver machinery.
- Solver output now includes a one-line problem summary, a final convergence summary, and
  clearer rejected-step markers. After a rejected step, `Δ` regrowth is capped at half the
  failed value.
- The Di Tella example no longer needs a hand-tuned initial `Δ`.

## [1.6.0] — 2026-07-03

- Add a documentation site built with Documenter and Literate, with a getting-started guide,
  a reorganized example gallery, solver-diagnostics documentation, and cross-links to
  InfinitesimalGenerators.jl.
- Rename the saved-outputs result field from `optional` to `saved`; `result.optional` remains
  as a backward-compatible alias.
- Add a risky-asset consumption-saving example.
- Allow state-grid dimensions to mix vector containers (e.g. a `range` and a `Vector`).
- Harden the monotonicity diagnostics.

## [1.5.0] — 2026-07-01

- Switch the nonlinear solver backend to `NonlinearSolve.jl` (Newton and trust-region methods).
- Add an opt-in monotonicity diagnostic (`check_monotonicity`) that flags finite-difference
  stencils with the wrong upwind sign in the assembled residual Jacobian.
- Replace the block-banded Jacobian prototypes with explicit sparse local-stencil patterns for
  1D, 2D, and 3D state grids, using FiniteDiff's colored sparse Jacobian path.
- Support bounded mixed-complementarity (optimal-stopping) solves even when no sparse Jacobian
  prototype is available.
- Remove the unused `method = :linearization` path and the direct `BlockArrays` /
  `BlockBandedMatrices` dependencies.
- Add a deterministic neoclassical growth model example and citation metadata (`CITATION.cff`).

## [1.4.0] — 2026-06-22

- Add monotone upwind stencils for second-order cross derivatives.
- Add 3D sparse Jacobian support and 3D boundary conditions.
- Refactor `differentiate` into a single N-dimensional implementation; harden tests and examples.

## [1.3.0] — 2026-03-29

- Fix bugs in the solver, 3D differentiation, and the Wachter example; clarify the sparsity code.

## [1.2.1] — 2025-10-13

- Finish removing the deprecated `SparseDiffTools` dependency.

## [1.2.0] — 2025-10-13

- Remove the deprecated `SparseDiffTools` dependency.
- Add He–Krishnamurthy and Brunnermeier–Sannikov examples and save more intermediate outputs.

## [1.1.1] — 2024-06-17

- Update `BlockBandedMatrices`; add a manual (non-`pdesolve`) example.

## [1.1.0] — 2024-05-28

- Add support for a third state dimension.

## [1.0.3] — 2023-04-02

- Maintenance release.

## [1.0.2] — 2023-03-22

- Add the (later-removed) `:linearization` method and a two-state example; clarify example papers.

## [1.0.1] — 2021-09-13

- Documentation and syntax updates for the NamedTuple interface.

## [1.0.0] — 2021-07-29

- Switch the PDE-function interface to `NamedTuple`s for states, unknowns, and outputs.
- Invert the sign convention of the time derivative and simplify upwinding.

## Earlier (0.x)

See the [git history](https://github.com/matthieugomez/EconPDEs.jl/commits/main) for releases
prior to 1.0.0.
