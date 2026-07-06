using EconPDEs

# The example scripts now import Plots for their figures, so this benchmark also needs Plots
# available in the active environment.

const EXAMPLES_DIR = joinpath(pkgdir(EconPDEs), "examples")

# PDE with 1 state variable
include(joinpath(EXAMPLES_DIR, "asset_pricing", "campbell_cochrane.jl"))
@time pdesolve(m, stategrid, guess)
# 0.002575 seconds (11.81 k allocations: 3.769 MiB)

# PDE with 2 state variables
include(joinpath(EXAMPLES_DIR, "asset_pricing", "bansal_yaron.jl"))
@time pdesolve(m, stategrid, guess)
# 0.007524 seconds (11.55 k allocations: 9.674 MiB)

# System of 4 PDEs with 1 state variable
include(joinpath(EXAMPLES_DIR, "asset_pricing", "garleanu_panageas.jl"))
@time pdesolve(m, stategrid, guess)
# 0.003823 seconds (11.01 k allocations: 7.179 MiB)

# System of 3 PDEs with 2 state variables
include(joinpath(EXAMPLES_DIR, "asset_pricing", "di_tella.jl"))
@time pdesolve(m, stategrid, guess; Δ = 1e-3)
# 0.636638 seconds (278.49 k allocations: 977.864 MiB, 14.91% gc time, 3.18% compilation time: 100% of which was recompilation)