using EconPDEs, Base.Test


println("Running tests:")

include("../examples/CampbellCochrane.jl")
try
	m = CampbellCochraneModel()
	state = initialize_state(m; n = 1000)
	y0 = initialize_y(m, state)
	result, distance = pde_solve(m, state, y0)
	@test distance <= 1e-5
	println("\t\033[1m\033[32mPASSED\033[0m: CampbellCochrane")
catch e
	println("\t\033[1m\033[31mFAILED\033[0m: CampbellCochrane")
	showerror(STDOUT, e, backtrace())
	rethrow(e)
end

include("../examples/BansalYaron.jl")
try
	m = BansalYaronModel()
	state = initialize_state(m; μn = 5, σn = 5)
	y0 = initialize_y(m, state)
	result, distance = pde_solve(m, state, y0)
	@test distance <= 1e-5
	println("\t\033[1m\033[32mPASSED\033[0m: BansalYaron")
catch e
	println("\t\033[1m\033[31mFAILED\033[0m: BansalYaron")
	showerror(STDOUT, e, backtrace())
	rethrow(e)
end


include("../examples/GarleanuPanageas.jl")
try
	m = GarleanuPanageasModel()
	state = initialize_state(m; n = 10)
	y0 = initialize_y(m, state)
	result, distance = pde_solve(m, state, y0)
	@test distance <= 1e-5
	println("\t\033[1m\033[32mPASSED\033[0m: GarleanuPanageas")
catch e
	println("\t\033[1m\033[31mFAILED\033[0m: GarleanuPanageas")
	showerror(STDOUT, e, backtrace())
	rethrow(e)
end


include("../examples/DiTella.jl")
try
	m = DiTellaModel()
	state = initialize_state(m ; xn = 10, νn = 3)
	y0 = initialize_y(m, state)
	result, distance = pde_solve(m, state, y0)
	@time pde_solve(m, state, y0, is_algebraic = (false, false, true))
	@test distance <= 1e-5
	println("\t\033[1m\033[32mPASSED\033[0m: DiTella")
catch e
	println("\t\033[1m\033[31mFAILED\033[0m: DiTella")
	showerror(STDOUT, e, backtrace())
	rethrow(e)
end

include("../examples/WangWangYang.jl")
try
	m = WangWangYangModel()
	state = initialize_state(m; n = 10)
	y0 = initialize_y(m, state)
	result, distance = pde_solve(m, state, y0)
	@test distance <= 1e-5
	println("\t\033[1m\033[32mPASSED\033[0m: WangWangYang")
catch e
	println("\t\033[1m\033[31mFAILED\033[0m: WangWangYang")
	showerror(STDOUT, e, backtrace())
	rethrow(e)
end


