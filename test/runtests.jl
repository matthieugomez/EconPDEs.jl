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
    @test_throws ArgumentError pdesolve(correct_upwind, grid_good, y_good; Δ = Inf, iterations = 1, verbose = false, method = :linearization)
end

@testset "Bounded solve without sparse prototype" begin
    function bounded_residual!(ydot, y)
        ydot[1] = y[1] - 2.0
        return ydot
    end

    y, residual_norm = finiteschemesolve(bounded_residual!, [0.5]; Δ = Inf, verbose = false, J0 = nothing, y̲ = [0.0], ȳ = [1.0])

    @test y[1] ≈ 1.0 atol = 1e-6
    @test residual_norm <= 1e-6
end

@testset "2D multi-function sparsity" begin
    grid = OrderedDict(:x => range(0.0, 1.0, length = 4), :z => range(0.0, 1.0, length = 5))
    y0 = OrderedDict(:A => ones(4, 5), :B => ones(4, 5))
    stategrid = EconPDEs.StateGrid(NamedTuple(grid))
    J = EconPDEs.sparse_jacobian(stategrid, y0)

    rows, cols, _ = EconPDEs.findnz(J)
    colors = EconPDEs.matrix_colors(J)

    @test size(J) == (40, 40)
    @test length(rows) == 520
    @test maximum(colors) <= 18
    for row in unique(rows)
        rowcols = cols[rows .== row]
        @test length(unique(colors[rowcols])) == length(rowcols)
    end
end

# The worked examples now live in examples/macro and examples/finance as Literate scripts.
# They import Plots and are executed (and thus tested) by the documentation build
# (see docs/make.jl), so they are no longer included here.
