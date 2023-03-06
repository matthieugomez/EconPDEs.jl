using Test

for x in (:CampbellCochrane, :Wachter, :BansalYaron, :GarleanuPanageas, :DiTella, :Haddad, :TuckmanVila)
	@testset "$x" include("../examples/AssetPricing/$(x).jl") 
end

for x in (:AchdouHanLasryLionsMoll_Diffusion, :AchdouHanLasryLionsMoll_DiffusionTwoAssets)
	@testset "$x" include("../examples/ConsumptionProblem/$(x).jl") 
end


for x in (:BoltonChenWang, )
	@testset "$x" include("../examples/InvestmentProblem/$(x).jl") 
end

for x in (:Leland, )
	@testset "$x" include("../examples/OptimalStoppingTime/$(x).jl") 
end


