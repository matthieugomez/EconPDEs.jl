using EconPDEs

# PDE with 1 state variable
include("/Users/Matthieu/Dropbox/Github/EconPDEs.jl/examples/AssetPricing/CampbellCochrane.jl")
m = CampbellCochraneModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
pdesolve(m, state, y0)
@time y, result, distance = pdesolve(m, state, y0)
# Old: 0.522157 seconds (6.75 M allocations: 1021.049 MiB, 16.41% gc time)
# With SparseDiffTools:    0.029206 seconds (361.40 k allocations: 23.384 MiB, 32.25% gc time)

# PDE with 2 state variables
include("/Users/Matthieu/Dropbox/Github/EconPDEs.jl/examples/AssetPricing/BansalYaron.jl")
m = BansalYaronModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
pdesolve(m, state, y0)
@time pdesolve(m, state, y0)
# Old: 1.456456 seconds (15.41 M allocations: 3.434 GiB, 15.67% gc time)
# With SparseDiffTools:   0.112325 seconds (639.39 k allocations: 104.911 MiB, 8.37% gc time)


# System of 4 PDEs with 1 state variable
include("/Users/Matthieu/Dropbox/Github/EconPDEs.jl/examples/AssetPricing/GarleanuPanageas.jl")
m = GarleanuPanageasModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
pdesolve(m, state, y0)
@time pdesolve(m, state, y0)
# Old: 1.134897 seconds (4.30 M allocations: 1.939 GiB, 10.29% gc time)
# With SparseDiffTools:   0.070427 seconds (262.73 k allocations: 88.128 MiB, 10.60% gc time)

# System of 3 PDEs with 2 state variables
include("/Users/Matthieu/Dropbox/Github/EconPDEs.jl/examples/AssetPricing/DiTella.jl")
m = DiTellaModel()
state = initialize_state(m; xn = 80, νn = 10)
y0 = initialize_y(m, state)
pdesolve(m, state, y0)
@time pdesolve(m, state, y0)
# Old:  12.711423 seconds (48.52 M allocations: 24.635 GiB, 10.62% gc time)
# With SparseDiffTools:  1.850576 seconds (4.14 M allocations: 2.418 GiB, 10.71% gc time)
