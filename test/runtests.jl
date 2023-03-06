using Test

for x in (:CampbellCochrane, :Wachter, :BansalYaron, :GarleanuPanageas, :DiTella, :Haddad, :TuckmanVila)
	@testset "$x" begin include("../examples/AssetPricing/$(x).jl") end
end

for x in (:AchdouHanLasryLionsMoll_Diffusion, :AchdouHanLasryLionsMoll_DiffusionTwoAssets)
	@testset "$x" begin include("../examples/ConsumptionProblem/$(x).jl") end
end


for x in (:BoltonChenWang, )
	@testset "$x" begin include("../examples/InvestmentProblem/$(x).jl") end
end

for x in (:Leland, )
	@testset "$x" begin include("../examples/OptimalStoppingTime/$(x).jl") end
end


