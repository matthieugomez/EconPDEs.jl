using EconPDEs

# PDE with 1 state variable
include("./AssetPricing/CampbellCochrane.jl")
@time pdesolve(m, stategrid, yend)
#   0.011112 seconds (286.45 k allocations: 12.268 MiB)


# PDE with 2 state variables
include("./AssetPricing/BansalYaron.jl")
@time pdesolve(m, stategrid, yend)
# 0.081095 seconds (485.84 k allocations: 78.910 MiB, 11.41% gc time)

# System of 4 PDEs with 1 state variable
include("./AssetPricing/GarleanuPanageas.jl")
@time pdesolve(m, stategrid, yend)
# 0.047346 seconds (167.95 k allocations: 65.049 MiB, 12.90% gc time)


# System of 3 PDEs with 2 state variables
include("./AssetPricing/DiTella.jl")
@time pdesolve(m, stategrid, yend)
# 0.912719 seconds (2.55 M allocations: 1.303 GiB, 9.01% gc time)

