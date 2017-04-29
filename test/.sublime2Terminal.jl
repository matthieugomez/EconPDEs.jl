m = GarleanuPanageasModel()
grid = StateGrid(m; n = 10)
y0 = initialize(m, grid)
result, distance = fullsolve(m, grid, y0)