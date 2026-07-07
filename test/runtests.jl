using Test
using Logging
using SparseArrays
using OrderedCollections: OrderedDict
using EconPDEs

#========================================================================================
Shared model: deterministic neoclassical growth (closed-form steady state to check).
========================================================================================#

const A = 0.5
const α = 0.3
const δ = 0.05
const ρ = 0.05
const γ = 2.0
const kbar = (α * A / (ρ + δ))^(1 / (1 - α))

growth_grid(n = 200) = (; k = range(0.1 * kbar, 5.0 * kbar, length = n))
growth_guess(grid) = (; v = [(A * k^α)^(1 - γ) / (1 - γ) / ρ for k in grid[:k]])

function growth_hjb(state::NamedTuple, u::NamedTuple)
    k = state.k
    cmax = 10 * A * k^α
    c_up = u.vk_up > 0 ? min(u.vk_up^(-1 / γ), cmax) : cmax
    μ_up = A * k^α - δ * k - c_up
    c_down = u.vk_down > 0 ? min(u.vk_down^(-1 / γ), cmax) : cmax
    μ_down = A * k^α - δ * k - c_down
    if μ_up > 0
        c, vk, μk = c_up, u.vk_up, μ_up
    elseif μ_down < 0
        c, vk, μk = c_down, u.vk_down, μ_down
    else
        c = A * k^α - δ * k
        vk, μk = c^(-γ), 0.0
    end
    vt = -(c^(1 - γ) / (1 - γ) + μk * vk - ρ * u.v)
    return (; vt), (; c, μk)
end

@testset "Stationary solve (neoclassical growth)" begin
    grid = growth_grid()
    result = pdesolve(growth_hjb, grid, growth_guess(grid); verbose = false)
    @test result.residual_norm <= 1e-6

    v = result.solution.v
    c = result.saved.c
    μk = result.saved.μk
    @test result.solution isa NamedTuple
    @test result.saved isa NamedTuple
    @test size(v) == size(grid.k)
    @test size(c) == size(grid.k)
    # value increasing and concave, consumption increasing
    @test all(diff(v) .> 0)
    @test all(diff(diff(v)) .< 1e-8)
    @test all(diff(c) .> 0)
    # drift crosses zero at the closed-form steady state
    icross = findfirst(<(0), μk)
    @test icross !== nothing
    @test abs(grid.k[icross] - kbar) < 0.05 * kbar
    # saved also contains the solved unknowns
    @test result.saved.v == v
    @test :saved in propertynames(result)
    @test :tolerance in propertynames(result, true)
    zero_alias = @test_deprecated result.zero
    optional_alias = @test_deprecated result.optional
    @test zero_alias === result.solution
    @test optional_alias === result.saved

    # tuple destructuring
    solution, residual_norm, saved = result
    @test solution === result.solution
    @test residual_norm === result.residual_norm
    @test saved === result.saved

    # compact show: names and sizes, not the arrays
    str = sprint(show, result)
    @test occursin("EconPDEResult", str)
    @test occursin("v (200)", str)
    @test length(str) < 500
end

@testset "OrderedDict inputs and partial bc" begin
    grid = growth_grid(100)
    result = pdesolve(growth_hjb, OrderedDict(:k => grid.k),
                      OrderedDict(:v => growth_guess(grid).v); verbose = false)
    @test result.residual_norm <= 1e-6

    # bc may cover any subset of (unknown, state) pairs; here the full set for `vk`
    result_bc = pdesolve(growth_hjb, grid, growth_guess(grid);
                         bc = (; vk = (0.0, 0.0)), verbose = false)
    @test result_bc.residual_norm <= 1e-6
end

@testset "Time-dependent solve" begin
    grid = growth_grid(100)
    guess = growth_guess(grid)
    τs = range(0, 100, length = 10)

    # two-argument function: same equation at each time
    result = pdesolve(growth_hjb, grid, guess, τs; verbose = false)
    # each solution array gains a trailing time dimension
    @test size(result.solution.v) == (length(grid.k), length(τs))
    @test result.solution.v[:, end] == collect(guess.v)   # terminal condition
    @test maximum(result.residual_norm[2:end]) <= 1e-6
    @test result.converged
    @test size(result.saved.c) == (length(grid.k), length(τs))

    # three-argument function: time-varying equation is detected and used
    ts_seen = Float64[]
    function hjb_t(state, u, t)
        isempty(ts_seen) || last(ts_seen) == t || push!(ts_seen, t)
        isempty(ts_seen) && push!(ts_seen, t)
        growth_hjb(state, u)
    end
    result_t = pdesolve(hjb_t, grid, guess, τs; verbose = false)
    @test maximum(result_t.residual_norm[2:end]) <= 1e-6
    @test length(unique(ts_seen)) > 1

    # typed time arguments are detected from actual applicability, not a broad `Number` method
    ts_seen_float = Float64[]
    function hjb_t_float(state, u, t::Float64)
        push!(ts_seen_float, t)
        growth_hjb(state, u)
    end
    result_t_float = pdesolve(hjb_t_float, grid, guess, τs; verbose = false)
    @test maximum(result_t_float.residual_norm[2:end]) <= 1e-6
    @test length(unique(ts_seen_float)) > 1

    # time-dependent solves accept the same inner-solve keywords users use for stationary solves
    result_inner = pdesolve(growth_hjb, grid, guess, τs; inner_maxiters = 20, inner_abstol = 1e-8, verbose = false)
    @test maximum(result_inner.residual_norm[2:end]) <= 1e-6

    # pseudo-transient step-size controls are stationary-only; a time grid supplies the step
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess, τs; Δ = 1.0, verbose = false)
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess, τs; scale = 2.0, verbose = false)
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess, τs; minΔ = 1e-6, verbose = false)
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess, τs; maxΔ = 10.0, verbose = false)

    no_saved_grid = (; x = range(0.0, 1.0, length = 5))
    no_saved_guess = (; v = zeros(5))
    no_saved_τs = range(0.0, 1.0, length = 3)
    no_saved_result = pdesolve(no_saved_grid, no_saved_guess, no_saved_τs; verbose = false) do state, u
        vt = -(u.v - state.x)
        (; vt)
    end
    @test no_saved_result.saved isa NamedTuple
    @test no_saved_result.saved == NamedTuple()

    # duplicate times would create a zero implicit step
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess, [0.0, 0.0, 1.0]; verbose = false)

    # an unconverged implicit step warns instead of being silently accepted
    @test_logs (:warn, r"did not converge") match_mode=:any pdesolve(growth_hjb, grid, guess, [0.0, 1.0]; maxiters = 0, verbose = false)
end

@testset "Failure reporting and verbosity" begin
    grid = growth_grid(50)
    guess = growth_guess(grid)

    # verbose = false is completely silent on success, so pdesolve can run inside
    # estimation loops without flooding the log
    result = mktemp() do path, io
        r = redirect_stdout(io) do
            pdesolve(growth_hjb, grid, guess; verbose = false)
        end
        flush(io)
        @test isempty(read(path, String))
        r
    end
    @test result.converged
    @test :converged in propertynames(result)
    @test occursin("converged:     true", sprint(show, result))

    # ...but a failed solve warns even with verbose = false
    @test_logs (:warn, r"did not converge") match_mode=:any pdesolve(growth_hjb, grid, guess; maxiters = 1, verbose = false)
    result_bad = with_logger(NullLogger()) do
        pdesolve(growth_hjb, grid, guess; maxiters = 1, verbose = false)
    end
    @test !result_bad.converged
    @test occursin("converged:     false", sprint(show, result_bad))
end

@testset "Input validation" begin
    grid = growth_grid(50)
    guess = growth_guess(grid)

    # plain Dict rejected (arbitrary iteration order)
    @test_throws ArgumentError pdesolve(growth_hjb, Dict(:k => grid.k), guess; verbose = false)
    @test_throws ArgumentError pdesolve(growth_hjb, grid, Dict(:v => guess.v); verbose = false)

    # guess shape mismatch names the unknown
    err = try pdesolve(growth_hjb, grid, (; v = zeros(3)); verbose = false); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("`v`", err.msg)

    # PDE function must return `<unknown>t`
    @test_throws ArgumentError pdesolve((s, u) -> (; wrong = 0.0), grid, guess; verbose = false)
    # ... and must return a NamedTuple
    @test_throws ArgumentError pdesolve((s, u) -> (0.0,), grid, guess; verbose = false)

    # each state dimension needs at least two grid points for finite differences
    @test_throws ArgumentError pdesolve(growth_hjb, (; k = [kbar]), (; v = [0.0]); verbose = false)

    # unknown bc entry
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess; bc = (; vx = (0.0, 0.0)), verbose = false)

    # removed solver keyword, checked upfront
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess; method = :trust_region, verbose = false)

    # is_algebraic must use the unknowns' names
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess; is_algebraic = (; w = true), verbose = false)
    # ...and accepts either a scalar Bool or a Bool array on the grid
    @test pdesolve(growth_hjb, grid, guess; is_algebraic = (; v = fill(false, size(guess.v))), verbose = false).converged
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess; is_algebraic = (; v = zeros(Int, size(guess.v))), verbose = false)

    # bounds may be keyed by unknown name at the pdesolve layer
    linear_grid = (; x = range(0.0, 1.0, length = 5))
    linear_guess = (; v = zeros(5))
    linear_pde = (s, u) -> (; vt = u.v - 2.0)
    bounded = pdesolve(linear_pde, linear_grid, linear_guess; Δ = Inf, lower_bound = (; v = fill(3.0, 5)), verbose = false)
    @test bounded.solution.v ≈ fill(3.0, 5)
    @test_throws ArgumentError pdesolve(linear_pde, linear_grid, linear_guess; lower_bound = (; w = zeros(5)), verbose = false)
    trust_region_result = pdesolve(linear_pde, linear_grid, linear_guess; Δ = Inf, alg = NonlinearSolve.TrustRegion(), verbose = false)
    @test trust_region_result.solution.v ≈ fill(2.0, 5)

    # ambiguous generated names: states (k, a) with unknowns (v, vk) collide on vka_up
    @test_throws ArgumentError pdesolve((s, u) -> (; vt = 0.0, vkt = 0.0),
                                        (; k = 0:0.1:1, a = 0:0.1:1),
                                        (; v = ones(11, 11), vk = ones(11, 11)); verbose = false)

    # removed 1.x keywords are rejected before reaching NonlinearSolve
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess; maxdist = 1e-6, verbose = false)
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess; innerdist = 1e-6, verbose = false)
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess; autodiff = :finite, verbose = false)
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess; algorithm = NonlinearSolve.TrustRegion(), verbose = false)
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess; iterations = 1, verbose = false)
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess; inner_iterations = 1, verbose = false)
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess; jac = (J, y) -> J, verbose = false)
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess; jac_prototype = spzeros(length(guess.v), length(guess.v)), verbose = false)
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess; colorvec = ones(Int, length(guess.v)), verbose = false)
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess; fdcache = nothing, verbose = false)
    @test_throws ArgumentError pdesolve(growth_hjb, grid, guess; monotonicity_check = nothing, verbose = false)
    @test_throws ArgumentError finiteschemesolve((ydot, y) -> (ydot .= y), [1.0]; method = :newton, verbose = false)
    @test_throws ArgumentError finiteschemesolve((ydot, y) -> (ydot .= y), [1.0]; autodiff = :finite, verbose = false)
    @test_throws ArgumentError finiteschemesolve((ydot, y) -> (ydot .= y), [1.0]; algorithm = NonlinearSolve.TrustRegion(), verbose = false)
    @test_throws ArgumentError finiteschemesolve((ydot, y) -> (ydot .= y), [1.0]; iterations = 1, verbose = false)
    @test_throws ArgumentError finiteschemesolve((ydot, y) -> (ydot .= y), [1.0]; inner_iterations = 1, verbose = false)
    @test_throws ArgumentError finiteschemesolve((ydot, y) -> (ydot .= y), [1.0]; maxdist = 1e-6, verbose = false)
    @test_throws ArgumentError finiteschemesolve((ydot, y) -> (ydot .= y), [1.0]; y̲ = [0.0], verbose = false)
    @test_throws ArgumentError finiteschemesolve((ydot, y) -> (ydot .= y), [1.0]; J0 = nothing, verbose = false)
    @test_throws ArgumentError finiteschemesolve((ydot, y) -> (ydot .= y), [1.0]; reformulation = :smooth, verbose = false)
    @test_throws ArgumentError finiteschemesolve((ydot, y) -> (ydot .= y), [1.0]; autoscale = false, verbose = false)
end

@testset "Mixed grid vector containers" begin
    grid = (; x = range(0.0, 1.0, length = 4),
              z = collect(range(0.0, 1.0, length = 5)))
    stategrid = EconPDEs.StateGrid(grid)
    @test size(stategrid) == (4, 5)
    @test stategrid[CartesianIndex(2, 3)] == (x = grid.x[2], z = grid.z[3])

    result = pdesolve(grid, (; v = zeros(4, 5)); Δ = Inf, verbose = false) do state, u
        vt = -(u.v - state.x - state.z)
        (; vt)
    end
    @test result.residual_norm <= 1e-8
    @test result.solution.v ≈ [x + z for x in grid.x, z in grid.z]
    @test result.saved isa NamedTuple
    @test result.saved == NamedTuple()
    @test !occursin("saved:", sprint(show, result))
    solution, residual_norm, saved = result
    @test solution === result.solution
    @test residual_norm === result.residual_norm
    @test saved === result.saved

    int_float_grid = (; x = 1:4, z = collect(range(0.0, 1.0, length = 5)))
    @test EconPDEs.StateGrid(int_float_grid)[CartesianIndex(1, 1)] == (x = 1.0, z = 0.0)
end

@testset "ND sparse coloring" begin
    grid = (; x = range(0.0, 1.0, length = 2),
              z = range(0.0, 1.0, length = 2),
              q = range(0.0, 1.0, length = 2),
              r = range(0.0, 1.0, length = 2))
    guess = (; v = zeros(2, 2, 2, 2))
    stategrid = EconPDEs.StateGrid(grid)
    J = EconPDEs.sparse_jacobian(stategrid, guess)
    colors = EconPDEs.stencil_colors(size(stategrid), length(guess))

    @test J !== nothing
    @test size(J) == (16, 16)
    @test maximum(colors) == 16

    rows, cols, _ = EconPDEs.findnz(J)
    for row in unique(rows)
        rowcols = cols[rows .== row]
        @test length(unique(colors[rowcols])) == length(rowcols)
    end

    result = pdesolve(grid, guess; Δ = Inf, verbose = false) do state, u
        vt = -(u.v - state.x - state.z - state.q - state.r)
        (; vt)
    end
    @test result.converged
    @test result.solution.v ≈ [x + z + q + r for x in grid.x, z in grid.z, q in grid.q, r in grid.r]
end

@testset "Monotonicity diagnostics" begin
    grid = OrderedDict(:x => range(0.0, 1.0, length = 2))
    y0 = OrderedDict(:V => collect(range(1.0, 2.0, length = 2)))

    function wrong_upwind(state, y)
        (; V, Vx_down) = y
        Vt = -(1.0 + V + Vx_down)
        return (; Vt)
    end

    @test_logs (:warn, r"Non-monotone finite-difference stencil") pdesolve(wrong_upwind, grid, y0; Δ = Inf, maxiters = 1, verbose = false, check_monotonicity = true)

    grid_good = OrderedDict(:x => range(0.0, 1.0, length = 3))
    y_good = OrderedDict(:V => collect(range(1.0, 2.0, length = 3)))

    function correct_upwind(state, y)
        (; V, Vx_up) = y
        Vt = -(1.0 + V + Vx_up)
        return (; Vt)
    end

    @test_logs min_level=Logging.Warn pdesolve(correct_upwind, grid_good, y_good; Δ = Inf, maxiters = 1, verbose = false, check_monotonicity = true)

    checker = EconPDEs.MonotonicityChecker(EconPDEs.StateGrid((; x = range(0.0, 1.0, length = 2))), (; V = ones(2)))
    @test EconPDEs._try_run_monotonicity_check!(checker, spzeros(2, 2), zeros(2), 1.0, Bool[]) === nothing
    @test checker.disabled
end

@testset "Sign-convention diagnostic" begin
    grid = OrderedDict(:x => range(0.0, 1.0, length = 5))
    y0 = OrderedDict(:V => collect(range(1.0, 2.0, length = 5)))

    # stationary equation 0 = 1 + Vx - V, with the sign of the returned time derivative flipped
    function flipped_sign(state, y)
        (; V, Vx_up) = y
        Vt = 1.0 + Vx_up - V
        return (; Vt)
    end

    # match_mode=:any because the truncated solve (maxiters = 1) now also warns about
    # non-convergence, even with verbose = false
    @test_logs (:warn, r"negative diagonal at every grid point") match_mode=:any pdesolve(flipped_sign, grid, y0; Δ = 0.5, maxiters = 1, verbose = false)

    # the correct convention triggers no warning under the default (finite Δ) solve
    function correct_sign(state, y)
        (; V, Vx_up) = y
        Vt = -(1.0 + Vx_up - V)
        return (; Vt)
    end

    @test_logs min_level=Logging.Warn pdesolve(correct_sign, grid, y0; verbose = false)
end

@testset "Bounded solve without sparse prototype" begin
    function bounded_residual!(ydot, y)
        ydot[1] = y[1] - 2.0
        return ydot
    end

    # the default bounded solve uses the minmax MCP reformulation
    y, residual_norm = finiteschemesolve(bounded_residual!, [0.5]; Δ = Inf, verbose = false, jac_prototype = nothing, lower_bound = [0.0], upper_bound = [1.0])

    @test y[1] ≈ 1.0 atol = 1e-6
    @test residual_norm <= 1e-6

    @test_throws ArgumentError finiteschemesolve(bounded_residual!, [0.5]; Δ = Inf, verbose = false, jac_prototype = nothing, lower_bound = [0.0], upper_bound = [1.0], reformulation = :smooth)
    @test_throws ArgumentError finiteschemesolve(bounded_residual!, [0.5]; Δ = Inf, verbose = false, jac_prototype = nothing, lower_bound = [0.0], upper_bound = [1.0], reformulation = :bogus)
end

@testset "Analytic Jacobian in finiteschemesolve" begin
    jac_calls = Ref(0)
    function residual!(ydot, y)
        ydot[1] = 2 * y[1] - 4
        ydot[2] = y[2] - 3
        return nothing
    end
    function residual_jac!(J, y)
        jac_calls[] += 1
        fill!(J, 0)
        J[1, 1] = 2
        J[2, 2] = 1
        return J
    end

    y, residual_norm = finiteschemesolve(residual!, [0.0, 0.0]; Δ = Inf, jac = residual_jac!, verbose = false)
    @test y ≈ [2.0, 3.0]
    @test residual_norm <= 1e-8
    @test jac_calls[] > 0

    y_sparse, sparse_residual_norm = finiteschemesolve(residual!, [0.0, 0.0]; Δ = Inf, jac = residual_jac!, jac_prototype = sparse([1.0 0.0; 0.0 1.0]), colorvec = [1], verbose = false)
    @test y_sparse ≈ [2.0, 3.0]
    @test sparse_residual_norm <= 1e-8

    jac_calls[] = 0
    ystep, step_residual_norm = EconPDEs.implicit_timestep(residual!, [0.0, 0.0], 1.0; is_algebraic = [false, true], jac = residual_jac!, jac_prototype = sparse([1.0 0.0; 0.0 1.0]), maxiters = 3, abstol = 1e-10, verbose = false)
    @test ystep ≈ [4 / 3, 3.0] atol = 1e-8
    @test step_residual_norm <= 1e-8
    @test jac_calls[] > 0

    y_bounded, bounded_residual_norm = finiteschemesolve(residual!, [0.0, 0.0]; Δ = Inf, jac = residual_jac!, jac_prototype = sparse([1.0 0.0; 0.0 1.0]), lower_bound = [0.0, 0.0], upper_bound = [1.0, 10.0], verbose = false)
    @test y_bounded ≈ [1.0, 3.0]
    @test bounded_residual_norm <= 1e-8
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

@testset "stencil_colors" begin
    # dimensions below 3 points exercise the reduced mixed-radix encoding
    for (s, F) in [((5,), 1), ((2,), 2), ((4, 5), 2), ((2, 5), 1), ((3, 4, 5), 2)]
        J = EconPDEs.local_stencil_jacobian(s, F)
        colors = EconPDEs.stencil_colors(s, F)
        rows, cols, _ = EconPDEs.findnz(J)
        # valid distance-2 coloring: within any row, nonzero columns have distinct colors
        for row in unique(rows)
            rowcols = cols[rows .== row]
            @test length(unique(colors[rowcols])) == length(rowcols)
        end
        # optimal: exactly the size of the largest clique of the column conflict graph
        @test maximum(colors) == F * prod(min.(3, s))
        # contiguous: no color class is empty, so no wasted function evaluation
        @test sort(unique(colors)) == 1:maximum(colors)
    end

    # a wrong-length colorvec is rejected before it reaches FiniteDiff
    @test_throws ArgumentError finiteschemesolve((ydot, y) -> (ydot .= y .- 2.0), [0.5, 0.5]; Δ = Inf, verbose = false, jac_prototype = EconPDEs.sparse([1.0 0.0; 0.0 1.0]), colorvec = [1])
end

# The worked examples live under examples/ as Literate scripts. They import Plots and are
# executed (and thus verified) by the documentation build (see docs/make.jl), so they are
# not run here.
