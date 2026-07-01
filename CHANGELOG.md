# Changelog

## Unreleased

- Replace block-banded Jacobian prototypes with explicit sparse local stencil
  patterns for 1D, 2D, and 3D state grids.
- Use FiniteDiff's native colored sparse Jacobian path when a sparse prototype
  is available, with `:central` retained as the central finite-difference
  variant.
- Support bounded mixed-complementarity solves even when no sparse Jacobian
  prototype is available.
- Remove the unused `method=:linearization` solver path.
- Remove direct `BlockArrays` and `BlockBandedMatrices` dependencies.
- Make `examples/runbenchmark.jl` resolve example files through
  `pkgdir(EconPDEs)` instead of the caller's current working directory.
