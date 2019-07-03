using EconPDEs

# PDE with 1 state variable
include("/Users/Matthieu/Dropbox/Github/EconPDEs.jl/examples/AssetPricing/CampbellCochrane.jl")
m = CampbellCochraneModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
pdesolve(m, state, y0)
@time y, result, distance = pdesolve(m, state, y0)
# Old: 0.522157 seconds (6.75 M allocations: 1021.049 MiB, 16.41% gc time)
# With SparseDiffTools:   0.053521 seconds (901.37 k allocations: 45.875 MiB, 29.16% gc time)

# PDE with 2 state variable
include("/Users/Matthieu/Dropbox/Github/EconPDEs.jl/examples/AssetPricing/BansalYaron.jl")
m = BansalYaronModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
pdesolve(m, state, y0)
@time pdesolve(m, state, y0)
# Old: 1.456456 seconds (15.41 M allocations: 3.434 GiB, 15.67% gc time)
# With SparseDiffTools: 0.186748 seconds (3.03 M allocations: 243.180 MiB, 14.58% gc time)


# System of 4 PDEs with 1 state variable
include("/Users/Matthieu/Dropbox/Github/EconPDEs.jl/examples/AssetPricing/GarleanuPanageas.jl")
m = GarleanuPanageasModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
pdesolve(m, state, y0)
@time pdesolve(m, state, y0)
# Old: 1.134897 seconds (4.30 M allocations: 1.939 GiB, 10.29% gc time)
# With SparseDiffTools: 0.078087 seconds (952.86 k allocations: 109.397 MiB, 18.61% gc time)

# System of 3 PDEs with 3 state variable
include("/Users/Matthieu/Dropbox/Github/EconPDEs.jl/examples/AssetPricing/DiTella.jl")
m = DiTellaModel()
state = initialize_state(m; xn = 80, νn = 10)
y0 = initialize_y(m, state)
pdesolve(m, state, y0)
@time pdesolve(m, state, y0)
# Old:  12.711423 seconds (48.52 M allocations: 24.635 GiB, 10.62% gc time)
# With SparseDiffTools:  2.341447 seconds (19.59 M allocations: 2.023 GiB, 16.00% gc time)



