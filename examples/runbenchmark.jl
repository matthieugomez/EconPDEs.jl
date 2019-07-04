using EconPDEs

# PDE with 1 state variable
include("./AssetPricing/CampbellCochrane.jl")
m = CampbellCochraneModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
pdesolve(m, stategrid, y0)
@time y, result, distance = pdesolve(m, stategrid, y0)
# Old: 0.522157 seconds (6.75 M allocations: 1021.049 MiB, 16.41% gc time)
# With SparseDiffTools:    0.024419 seconds (361.35 k allocations: 23.222 MiB, 15.98% gc time)

# PDE with 2 state variables
include("./AssetPricing/BansalYaron.jl")
m = BansalYaronModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
pdesolve(m, stategrid, y0)
@time pdesolve(m, stategrid, y0)
# Old: 1.456456 seconds (15.41 M allocations: 3.434 GiB, 15.67% gc time)
# With SparseDiffTools:     0.099714 seconds (606.58 k allocations: 101.652 MiB, 8.60% gc time)

# System of 4 PDEs with 1 state variable
include("./AssetPricing/GarleanuPanageas.jl")
m = GarleanuPanageasModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
pdesolve(m, stategrid, y0)
@time pdesolve(m, stategrid, y0)
# Old: 1.134897 seconds (4.30 M allocations: 1.939 GiB, 10.29% gc time)
# With SparseDiffTools:    0.058315 seconds (210.77 k allocations: 82.124 MiB, 10.73% gc time)

# System of 3 PDEs with 2 state variables
include("./AssetPricing/DiTella.jl")
m = DiTellaModel()
stategrid = initialize_stategrid(m; xn = 80, νn = 10)
y0 = initialize_y(m, stategrid)
pdesolve(m, stategrid, y0)
@time pdesolve(m, stategrid, y0)
# Old:  12.711423 seconds (48.52 M allocations: 24.635 GiB, 10.62% gc time)
# With SparseDiffTools:    1.744255 seconds (2.88 M allocations: 2.343 GiB, 11.40% gc time)
