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



