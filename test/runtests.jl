using Test

for x in (:CampbellCochrane, :Wachter, :BansalYaron, :GarleanuPanageas, :DiTella, :Haddad, :TuckmanVila)
	@testset begin "$x" include("../examples/AssetPricing/$(x).jl") end
end

for x in (:AchdouHanLasryLionsMoll_Diffusion, :AchdouHanLasryLionsMoll_DiffusionTwoAssets)
	@testset begin "$x" include("../examples/ConsumptionProblem/$(x).jl") end
end


for x in (:BoltonChenWang, )
	@testset begin "$x" include("../examples/InvestmentProblem/$(x).jl") end
end

for x in (:Leland, )
	@testset begin "$x" include("../examples/OptimalStoppingTime/$(x).jl") end
end


