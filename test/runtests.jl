using EconPDEs, Test


println("Running tests:")

include("../examples/AssetPricing/Habit.jl")
try
	m = HabitModel()
	state = initialize_state(m; n = 1000)
	y0 = initialize_y(m, state)
	y, a, distance = pdesolve(m, state, y0)
	@test distance <= 1e-5
	println("\t\033[1m\033[32mPASSED\033[0m: Habit")
catch e
	println("\t\033[1m\033[31mFAILED\033[0m: Habit")
	showerror(stdout, e, backtrace())
	rethrow(e)
end

include("../examples/AssetPricing/LongRunRisk.jl")
try
	m = LongRunRiskModel()
	state = initialize_state(m; μn = 5, vn = 5)
	y0 = initialize_y(m, state)
	y, a, distance = pdesolve(m, state, y0)
	@test distance <= 1e-5
	println("\t\033[1m\033[32mPASSED\033[0m: LongRunRisk")
catch e
	println("\t\033[1m\033[31mFAILED\033[0m: LongRunRisk")
	showerror(stdout, e, backtrace())
	rethrow(e)
end


include("../examples/AssetPricing/Disaster.jl")
try
	m = DisasterModel()
	state = initialize_state(m; n = 5)
	y0 = initialize_y(m, state)
	y, a, distance = pdesolve(m, state, y0)
	@test distance <= 1e-5
	println("\t\033[1m\033[32mPASSED\033[0m: Disaster")
catch e
	println("\t\033[1m\033[31mFAILED\033[0m: Disaster")
	showerror(stdout, e, backtrace())
	rethrow(e)
end


include("../examples/AssetPricing/GarleanuPanageas.jl")
try
	m = GarleanuPanageasModel()
	state = initialize_state(m; n = 10)
	y0 = initialize_y(m, state)
	y, a, distance = pdesolve(m, state, y0)
	@test distance <= 1e-5
	println("\t\033[1m\033[32mPASSED\033[0m: GarleanuPanageas")
catch e
	println("\t\033[1m\033[31mFAILED\033[0m: GarleanuPanageas")
	showerror(stdout, e, backtrace())
	rethrow(e)
end


include("../examples/AssetPricing/DiTella.jl")
try
	m = DiTellaModel()
	state = initialize_state(m ; xn = 10, νn = 3)
	y0 = initialize_y(m, state)
	y, a, distance = pdesolve(m, state, y0)
	@time pdesolve(m, state, y0, is_algebraic = OrderedDict(:pA => false, :pB => false, :p => true))
	@test distance <= 1e-5
	println("\t\033[1m\033[32mPASSED\033[0m: DiTella")
catch e
	println("\t\033[1m\033[31mFAILED\033[0m: DiTella")
	showerror(stdout, e, backtrace())
	rethrow(e)
end

include("../examples/ConsumptionProblem/WangWangYang.jl")
try
	m = WangWangYangModel()
	state = initialize_state(m; n = 10)
	y0 = initialize_y(m, state)
	y, a, distance = pdesolve(m, state, y0)
	@test distance <= 1e-5
	println("\t\033[1m\033[32mPASSED\033[0m: WangWangYang")
catch e
	println("\t\033[1m\033[31mFAILED\033[0m: WangWangYang")
	showerror(stdout, e, backtrace())
	rethrow(e)
end


include("../examples/AssetPricing/ArbitrageHoldingCosts.jl")
try
	m = ArbitrageHoldingCosts()
	state = initialize_state(m; n = 10)
	y0 = initialize_y(m, state)
	τs = range(0, stop = 1.0, length = 5)
	y, a, distance = pdesolve(m, state, y0, τs)
	@test distance <= 1e-5
	println("\t\033[1m\033[32mPASSED\033[0m: HoldingCost")
catch e
	println("\t\033[1m\033[31mFAILED\033[0m: HoldingCost")
	showerror(stdout, e, backtrace())
	rethrow(e)
end


include("../examples/ConsumptionProblem/AchdouHanLasryLionsMoll_OneAsset.jl")
try
	m = AchdouHanLasryLionsMollModel()
	state = initialize_state(m; yn = 3, an = 5)
	y0 = initialize_y(m, state)
	y, a, distance = pdesolve(m, state, y0)
	@test distance <= 1e-5
	println("\t\033[1m\033[32mPASSED\033[0m: AchdouHanLasryLionsMoll")
catch e
	println("\t\033[1m\033[31mFAILED\033[0m: AchdouHanLasryLionsMoll")
	showerror(stdout, e, backtrace())
	rethrow(e)
end

include("../examples/ConsumptionProblem/AchdouHanLasryLionsMoll_TwoAssets.jl")
try
	m = AchdouHanLasryLionsMoll_TwoAssetsModel(amax = 10.0)
	state = initialize_state(m; yn = 3, an = 10)
	y0 = initialize_y(m, state)
	y, a, distance = pdesolve(m, state, y0)
	@test distance <= 1e-5
	println("\t\033[1m\033[32mPASSED\033[0m: AchdouHanLasryLionsMoll Two Assets")
catch e
	println("\t\033[1m\033[31mFAILED\033[0m: AchdouHanLasryLionsMoll Two Assets")
	showerror(stdout, e, backtrace())
	rethrow(e)
end


include("../examples/InvestmentProblem/BoltonChenWang.jl")
try
	m = BoltonChenWangModel()
	state = initialize_state(m; n = 10)
	y0 = initialize_y(m, state)
	y, a, distance = pdesolve(m, state, y0, bc = OrderedDict(:vw => (1.5, 1.0)))
	@test distance <= 1e-5
	println("\t\033[1m\033[32mPASSED\033[0m: BoltonChenWang")
catch e
	println("\t\033[1m\033[31mFAILED\033[0m: BoltonChenWang")
	showerror(stdout, e, backtrace())
	rethrow(e)
end




