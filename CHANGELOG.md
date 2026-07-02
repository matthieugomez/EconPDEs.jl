# Changelog

All notable changes to EconPDEs.jl are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/), and the project follows
[semantic versioning](https://semver.org/).

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
