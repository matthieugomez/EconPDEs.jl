# Changelog

All notable changes to EconPDEs.jl are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/), and the project follows
[semantic versioning](https://semver.org/).

## [Unreleased]

- The Jacobian coloring used by the sparse finite differences is now computed in closed
  form from the stencil structure (coordinate mod 3 in each state dimension, crossed with
  the unknown index) instead of by greedy graph coloring of the assembled pattern. Same
  provably optimal number of colors, but orders of magnitude faster to build (×300 on a
  200×200 grid, ×3,000 on a 50×50×50 grid) — noticeable when `pdesolve` is called
  repeatedly, e.g. in an estimation loop. `finiteschemesolve` gains a `colorvec` keyword
  for callers who supply their own `J0` and know its coloring; by default a user-supplied
  `J0` is still greedily colored.
- Pass the residual to the solver behind a compile-once barrier (`ResidualWrapper`) when a
  sparsity pattern is available (1–3 states), so the solver stack (`finiteschemesolve`,
  `implicit_timestep`, FiniteDiff, NLsolve/NonlinearSolve internals) no longer recompiles
  for every new model function. First solve of a new model in a warm session drops from
  ~0.7s to ~0.3s; the remainder is the model's own PDE function and derivative accessors.
- Add a `PrecompileTools` workload that solves a tiny model once at package precompilation,
  caching the model-independent machinery (sparsity pattern and coloring, Jacobian cache,
  Newton solve internals). First `pdesolve` call drops from ~5s to ~0.6s; the remainder is
  code specialized on the user's model (the generated derivative accessors and the PDE
  function itself), which must compile per model.
- Better solver output: `pdesolve` prints a one-line problem summary (unknowns, grid size)
  before solving and a convergence summary (time steps, elapsed time) after. A step whose
  inner solve fails prints ✗ in the residual column (no step is taken; it is retried with
  `Δ/10`) instead of a `NaN` residual, with unusual causes (a `NaN` from the model, a
  singular Jacobian) flagged in parentheses.
- Convergence failures are reported with `@warn` even when `verbose = false`, with
  actionable hints, so solves embedded in silent loops (e.g. estimation) stay quiet on
  success but surface failures. A solve that converges exactly at the iteration limit is
  no longer mislabeled as unconverged.
- `EconPDEResult` records the solver tolerance and gains a computed `converged` property,
  also shown in its compact display.
- After a rejected step, `Δ`'s regrowth is capped at half the failed value until the
  residual falls to half its level at the failure, breaking the grow → reject → shrink
  cycle on hard problems (each rejection wastes a full round of inner Newton iterations).
- Simplify the Di Tella example: drop the hand-tuned initial `Δ` (the default now
  converges) and the printed numeric callouts.

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
