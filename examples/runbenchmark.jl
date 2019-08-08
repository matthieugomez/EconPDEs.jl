using EconPDEs

# PDE with 1 state variable
include("./AssetPricing/CampbellCochrane.jl")
m = CampbellCochraneModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
pdesolve(m, stategrid, y0)
@time y, result, distance = pdesolve(m, stategrid, y0)
# 0.019167 seconds (249.93 k allocations: 15.914 MiB, 29.78% gc time)

# PDE with 2 state variables
include("./AssetPricing/BansalYaron.jl")
m = BansalYaronModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
pdesolve(m, stategrid, y0)
@time pdesolve(m, stategrid, y0)
# 0.081095 seconds (485.84 k allocations: 78.910 MiB, 11.41% gc time)

# System of 4 PDEs with 1 state variable
include("./AssetPricing/GarleanuPanageas.jl")
m = GarleanuPanageasModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
pdesolve(m, stategrid, y0)
@time pdesolve(m, stategrid, y0)
# 0.047346 seconds (167.95 k allocations: 65.049 MiB, 12.90% gc time)


# System of 3 PDEs with 2 state variables
include("./AssetPricing/DiTella.jl")
m = DiTellaModel()
stategrid = initialize_stategrid(m; xn = 80, νn = 10)
y0 = initialize_y(m, stategrid)
pdesolve(m, stategrid, y0)
@time pdesolve(m, stategrid, y0)
# 0.912719 seconds (2.55 M allocations: 1.303 GiB, 9.01% gc time)

