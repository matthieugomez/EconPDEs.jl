using Test
using Logging
using EconPDEs

@testset "Monotonicity diagnostics" begin
    grid = OrderedDict(:x => range(0.0, 1.0, length = 2))
    y0 = OrderedDict(:V => collect(range(1.0, 2.0, length = 2)))

    function wrong_upwind(state, y)
        (; V, Vx_down) = y
        Vt = -(1.0 + V + Vx_down)
        return (; Vt)
    end

    @test_logs (:warn, r"Non-monotone finite-difference stencil") pdesolve(wrong_upwind, grid, y0; Δ = Inf, iterations = 1, verbose = false, check_monotonicity = true)

    grid_good = OrderedDict(:x => range(0.0, 1.0, length = 3))
    y_good = OrderedDict(:V => collect(range(1.0, 2.0, length = 3)))

    function correct_upwind(state, y)
        (; V, Vx_up) = y
        Vt = -(1.0 + V + Vx_up)
        return (; Vt)
    end

    @test_logs min_level=Logging.Warn pdesolve(correct_upwind, grid_good, y_good; Δ = Inf, iterations = 1, verbose = false, check_monotonicity = true)
end

for x in (:CampbellCochrane, :Wachter, :BansalYaron, :GarleanuPanageas, :DiTella, :Haddad, :TuckmanVila, :HeKrishnamurthy, :BrunnermeirSannikov)
	@testset "$x" begin include("../examples/AssetPricing/$(x).jl") end
end

for x in (:AchdouHanLasryLionsMoll_Diffusion, :AchdouHanLasryLionsMoll_DiffusionTwoAssets, :AchdouHanLasryLionsMoll_TwoStates, :WangWangYang)
	@testset "$x" begin include("../examples/ConsumptionProblem/$(x).jl") end
end


for x in (:BoltonChenWang, )
	@testset "$x" begin include("../examples/InvestmentProblem/$(x).jl") end
end

for x in (:NeoclassicalGrowthModel, )
	@testset "$x" begin include("../examples/GrowthModel/$(x).jl") end
end

for x in (:Leland, )
	@testset "$x" begin include("../examples/OptimalStoppingTime/$(x).jl") end
end
