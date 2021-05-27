using Test


println("Running tests:")

for x in (:CampbellCochrane, :Wachter, :BansalYaron, :GarleanuPanageas, :DiTella, :Haddad, :ArbitrageHoldingCosts)
	try
		include("../examples/AssetPricing/$(x).jl")
		@test distance <= 1e-5
		println("\t\033[1m\033[32mPASSED\033[0m: $(x)")
	catch e
		println("\t\033[1m\033[31mFAILED\033[0m: $(x)")
		showerror(stdout, e, backtrace())
		rethrow(e)
	end
end

# Test Campbell Cochrane with fixed time grid
m = CampbellCochraneModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
τs = range(1000, stop = 0, length = 50)
y, result, distance = pdesolve(m, stategrid, y0)
y2, result, distance = pdesolve(m, stategrid, y0, τs)
@test sum(abs2, y[:p] .- y2[:p][:, end]) <= 1e-10



for x in (:WangWangYang, :AchdouHanLasryLionsMoll_OneAsset, :AchdouHanLasryLionsMoll_TwoAssets)
	try
		include("../examples/ConsumptionProblem/$(x).jl")
		@test distance <= 1e-5
		println("\t\033[1m\033[32mPASSED\033[0m: $(x)")
	catch e
		println("\t\033[1m\033[31mFAILED\033[0m: $(x)")
		showerror(stdout, e, backtrace())
		rethrow(e)
	end
end

for x in (:BoltonChenWang, )
	try
		include("../examples/InvestmentProblem/$(x).jl")
		@test distance <= 1e-5
		println("\t\033[1m\033[32mPASSED\033[0m: $(x)")
	catch e
		println("\t\033[1m\033[31mFAILED\033[0m: $(x)")
		showerror(stdout, e, backtrace())
		rethrow(e)
	end
end

for x in (:Leland, )
	try
		include("../examples/OptimalStoppingTime/$(x).jl")
		@test distance <= 1e-5
		println("\t\033[1m\033[32mPASSED\033[0m: $(x)")
	catch e
		println("\t\033[1m\033[31mFAILED\033[0m: $(x)")
		showerror(stdout, e, backtrace())
		rethrow(e)
	end
end


