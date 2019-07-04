using EconPDEs

# PDE with 1 state variable
include("./AssetPricing/CampbellCochrane.jl")
m = CampbellCochraneModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
pdesolve(m, stategrid, y0)
@time y, result, distance = pdesolve(m, stategrid, y0)
# Old: 0.522157 seconds (6.75 M allocations: 1021.049 MiB, 16.41% gc time)
# With SparseDiffTools:    0.029206 seconds (361.40 k allocations: 23.384 MiB, 32.25% gc time)

# PDE with 2 state variables
include("./AssetPricing/BansalYaron.jl")
m = BansalYaronModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
pdesolve(m, stategrid, y0)
@time pdesolve(m, stategrid, y0)
# Old: 1.456456 seconds (15.41 M allocations: 3.434 GiB, 15.67% gc time)
# With SparseDiffTools:   0.112325 seconds (639.39 k allocations: 104.911 MiB, 8.37% gc time)

# System of 4 PDEs with 1 state variable
include("./AssetPricing/GarleanuPanageas.jl")
m = GarleanuPanageasModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
pdesolve(m, stategrid, y0)
@time pdesolve(m, stategrid, y0)
# Old: 1.134897 seconds (4.30 M allocations: 1.939 GiB, 10.29% gc time)
# With SparseDiffTools:   0.070427 seconds (262.73 k allocations: 88.128 MiB, 10.60% gc time)

# System of 3 PDEs with 2 state variables
include("./AssetPricing/DiTella.jl")
m = DiTellaModel()
stategrid = initialize_stategrid(m; xn = 80, νn = 10)
y0 = initialize_y(m, stategrid)
pdesolve(m, stategrid, y0)
@time pdesolve(m, stategrid, y0)
# Old:  12.711423 seconds (48.52 M allocations: 24.635 GiB, 10.62% gc time)
# With SparseDiffTools:  1.850576 seconds (4.14 M allocations: 2.418 GiB, 10.71% gc time)