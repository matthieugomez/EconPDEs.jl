##############################################################################
##
## Non Linear solver using Pseudo-Transient Continuation Method
##
##############################################################################

# Implicit time step
function _nonlinear_autodiff(autodiff)
    if autodiff == :forward
        return AutoForwardDiff()
    elseif autodiff == :finite
        return AutoFiniteDiff()
    elseif autodiff == :central
        return AutoFiniteDiff(fdjtype = Val(:central))
    else
        return autodiff
    end
end

function _nonlinear_algorithm(method, autodiff; linsolve = \, concrete_jac = true)
    ad = _nonlinear_autodiff(autodiff)
    if method == :newton
        return NewtonRaphson(; autodiff = ad, linsolve = linsolve, concrete_jac = concrete_jac)
    elseif method == :trust_region
        return TrustRegion(; autodiff = ad, linsolve = linsolve, concrete_jac = concrete_jac)
    elseif method == :linearization
        return NewtonRaphson(; autodiff = ad, linsolve = linsolve, concrete_jac = concrete_jac)
    else
        throw(ArgumentError("method must be :newton, :trust_region, or :linearization"))
    end
end

_has_bounds(y̲, ȳ) = any(x -> x != -Inf, y̲) || any(x -> x != Inf, ȳ)

function _nonlinear_solve(G!, y0; jac = nothing, jac_prototype = nothing, iterations = 100, verbose = true, method = :newton, autodiff = :forward, maxdist = sqrt(eps()), linsolve = \, kwargs...)
    G_solve! = (du, u, p) -> G!(du, u)
    f = jac === nothing ?
        NonlinearFunction{true}(G_solve!; jac_prototype = jac_prototype) :
        NonlinearFunction{true}(G_solve!; jac = (J, u, p) -> jac(J, u), jac_prototype = jac_prototype)
    problem = NonlinearProblem(f, y0)
    algorithm = _nonlinear_algorithm(method, autodiff; linsolve = linsolve, concrete_jac = true)
    result = try
        solve(problem, algorithm; maxiters = iterations, abstol = maxdist, verbose = verbose, kwargs...)
    catch err
        err isa SingularException || rethrow()
        return y0, Inf
    end
    return result.u, norm(result.resid) / length(result.resid)
end

function implicit_timestep(G!, ypost, Δ; is_algebraic = fill(false, size(ypost)...), iterations = 100, verbose = true, method = :newton, autodiff = :forward, maxdist = sqrt(eps()), J0 = nothing, y̲ = fill(-Inf, length(ypost)), ȳ = fill(Inf, length(ypost)), reformulation = :smooth, autoscale = true, linsolve = \, monotonicity_check = nothing, kwargs...)
    #
    #if any(is_algebraic)
    # One option is 
    # G_helper!(ydot, y) = (G!(ydot, y) ; ydot .*= 1 .+ .!is_algebraic .* (Δ - 1) ; ydot .-= .!is_algebraic .* (ypost .- y))
    # however does not work if Δ is Inf
    G_helper!(ydot, y) = (G!(ydot, y) ; ydot .-= .!is_algebraic .* (ypost .- y) ./ Δ)
    if J0 === nothing
        if method == :linearization
            method = :newton
        end
        zero, residual_norm = _nonlinear_solve(G_helper!, ypost; iterations = iterations, verbose = verbose, method = method, autodiff = autodiff, maxdist = maxdist, linsolve = linsolve, kwargs...)
    else
        # remove forwarddiff path and use FiniteDiff everytime
        colorvec = matrix_colors(J0)
        J0_sparse = J0 isa Tridiagonal ? J0 : sparse(J0)
        bbbcache = JacobianCache(ypost, colorvec = colorvec, sparsity = J0_sparse)
        function j_helper!(J, y)
            finite_difference_jacobian!(J, G_helper!, y, bbbcache)
            _run_monotonicity_check!(monotonicity_check, J, y)
            return J
        end
        if _has_bounds(y̲, ȳ)
            # HJBVI bounds are mixed complementarity conditions, not box-constrained roots.
            result = mcpsolve(OnceDifferentiable(G_helper!, j_helper!, deepcopy(ypost), deepcopy(ypost), J0_sparse), y̲, ȳ, ypost; iterations = iterations, show_trace = verbose, ftol = maxdist, method = method, reformulation = reformulation)
            zero, residual_norm = result.zero, result.residual_norm
        elseif method == :linearization
            finite_difference_jacobian!(J0_sparse, G!, ypost, bbbcache)
            _run_monotonicity_check!(monotonicity_check, J0_sparse, ypost)
            GV = deepcopy(ypost)
            G!(GV, ypost)
            zero = (I + Δ .* J0_sparse) \ (ypost .- Δ .* (GV .- J0_sparse * ypost))
            residual_norm = 0.0
        else
            zero, residual_norm = _nonlinear_solve(G_helper!, ypost; jac = j_helper!, jac_prototype = J0_sparse, iterations = iterations, verbose = verbose, method = method, autodiff = autodiff, maxdist = maxdist, linsolve = linsolve, kwargs...)
        end
    end
    return zero, residual_norm
end

# Solve for steady state
function finiteschemesolve(G!, y0; Δ = 1.0, is_algebraic = fill(false, size(y0)...), iterations = 100, inner_iterations = 10, verbose = true, inner_verbose = false, method = :newton, autodiff = :forward, maxdist = sqrt(eps()), innerdist = sqrt(eps()), scale = 10.0, J0 = nothing, minΔ = 1e-9, y̲ = fill(-Inf, length(y0)), ȳ = fill(Inf, length(y0)), reformulation = :smooth, maxΔ = Inf, autoscale = true, linsolve = \, monotonicity_check = nothing, kwargs...)
    ypost = y0
    ydot = zero(y0)
    # check that does not return NAN or zero
    G!(ydot, ypost)
    residual_norm = norm(ydot) / length(ydot)
    isnan(residual_norm) && throw(ArgumentError("G! returns NaN with the initial value"))
    if residual_norm <= maxdist
        verbose && @warn "G! already returns zero with the initial value"
        return ypost, residual_norm
    end
    if Δ == Inf
        ypost, residual_norm = implicit_timestep(G!, y0, Δ; is_algebraic = is_algebraic, verbose = verbose, iterations = iterations,  method = method, autodiff = autodiff, maxdist = maxdist, J0 = J0, y̲ = y̲, ȳ = ȳ, linsolve = linsolve, monotonicity_check = monotonicity_check, kwargs...)
    else
        coef = 1.0
        oldresidual_norm = residual_norm
        iter = 0
        if verbose
            @printf "Iter   TimeStep   Residual\n"
            @printf "---- ---------- ----------\n"
        end
        while (iter < iterations) && (Δ >= minΔ) && (residual_norm > maxdist)
            iter += 1
            y, nlresidual_norm = implicit_timestep(G!, ypost, Δ; is_algebraic = is_algebraic, verbose = inner_verbose, iterations = inner_iterations, method = method, autodiff = autodiff, maxdist = innerdist, J0 = J0, y̲ = y̲, ȳ = ȳ, reformulation = reformulation, linsolve = linsolve, monotonicity_check = monotonicity_check, kwargs...)
            G!(ydot, y)
            if _has_bounds(y̲, ȳ)
                mask = y̲ .+ eps() .<= y .<= ȳ .- eps() # only unconstrained ydot is relevant for residual_norm calculation
                residual_norm, oldresidual_norm = norm(ydot .* (mask))/sum(mask), residual_norm
            else
                residual_norm, oldresidual_norm = norm(ydot) / length(ydot), residual_norm
            end
            residual_norm = isnan(residual_norm) ? Inf : residual_norm
            if nlresidual_norm <= innerdist
                # if the implicit time step is correctly solved
                if verbose
                    @printf "%4d %8.4e %8.4e\n" iter Δ residual_norm
                end
                coef = (residual_norm <= oldresidual_norm) ? scale * coef : 1.0
                Δ = min(Δ * coef * oldresidual_norm / residual_norm, maxΔ)
                ypost, y = y, ypost
            else
                if verbose
                    @printf "%4d %8.4e %8.4e\n" iter Δ NaN
                end
                # verbose && @show iter, Δ, NaN
                # if the implict time step is not solved
                # revert and diminish the time step
                coef = 1.0
                Δ = Δ / 10
                residual_norm = oldresidual_norm
            end
        end
    end
    verbose && (iter >= iterations) && @warn "Algorithm did not converge: Iter higher than the limit $(iterations)"
    verbose && (Δ < minΔ) && @warn "Algorithm did not converge: TimeStep lower than the limit $(minΔ)"
    return ypost, residual_norm
end
