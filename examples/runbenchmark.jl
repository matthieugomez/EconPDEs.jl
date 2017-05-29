using EconPDEs


include("CampbellCochrane.jl")
m = CampbellCochraneModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
pdesolve(m, grid, y0)
@time pdesolve(m, grid, y0)
#0.45s, 93MB

include("BansalYaron.jl")
m = BansalYaronModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
pdesolve(m, grid, y0)
@time pdesolve(m, grid, y0)
#1.13s, 150MB

include("GarleanuPanageas.jl")
m = GarleanuPanageasModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
pdesolve(m, grid, y0)
@time pdesolve(m, grid, y0)
#1.9s, 150MB

include("DiTella.jl")
m = DiTellaModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
pdesolve(m, grid, y0)
@time pdesolve(m, grid, y0)
# 30.820142s, 10.2 GB
pdesolve(m, grid, y0, is_algebraic = (false, false, true))

include("WangWangYang.jl")
m = WangWangYangModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
pdesolve(m, grid, y0)
@time pdesolve(m, grid, y0)
# 0.08s, 3.0MB



