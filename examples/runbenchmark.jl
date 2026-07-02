using EconPDEs

# The example scripts now import Plots for their figures, so this benchmark also needs Plots
# available in the active environment.

const EXAMPLES_DIR = joinpath(pkgdir(EconPDEs), "examples")

# PDE with 1 state variable
include(joinpath(EXAMPLES_DIR, "asset_pricing", "campbell_cochrane.jl"))
@time pdesolve(m, stategrid, yend)

# PDE with 2 state variables
include(joinpath(EXAMPLES_DIR, "asset_pricing", "bansal_yaron.jl"))
@time pdesolve(m, stategrid, yend)

# System of 4 PDEs with 1 state variable
include(joinpath(EXAMPLES_DIR, "asset_pricing", "garleanu_panageas.jl"))
@time pdesolve(m, stategrid, yend)

# System of 3 PDEs with 2 state variables
include(joinpath(EXAMPLES_DIR, "asset_pricing", "di_tella.jl"))
@time pdesolve(m, stategrid, yend)
