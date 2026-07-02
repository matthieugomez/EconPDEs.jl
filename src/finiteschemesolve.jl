##############################################################################
##
## Non Linear solver using Pseudo-Transient Continuation Method
##
##############################################################################

_has_bounds(y̲, ȳ) = any(x -> x != -Inf, y̲) || any(x -> x != Inf, ȳ)

function implicit_timestep(G!, ypost, Δ; is_algebraic = fill(false, size(ypost)...), iterations = 100, verbose = true, method = :newton, autodiff = :forward, maxdist = sqrt(eps()), J0 = nothing, y̲ = fill(-Inf, length(ypost)), ȳ = fill(Inf, length(ypost)), reformulation = :smooth, autoscale = true, monotonicity_check = nothing, kwargs...)
    method in (:newton, :trust_region) || throw(ArgumentError("method must be :newton or :trust_region"))
    G_helper!(ydot, y) = (G!(ydot, y) ; ydot .-= .!is_algebraic .* (ypost .- y) ./ Δ)

    jac = nothing
    J0c = J0 === nothing ? nothing : sparse(J0)
    if J0c !== nothing
        colorvec = matrix_colors(J0c)
        fdtype = autodiff == :central ? Val(:central) : Val(:forward)
        fdcache = JacobianCache(ypost, fdtype, eltype(ypost); colorvec = colorvec, sparsity = J0c)
        function jac!(J, y)
            finite_difference_jacobian!(J, G_helper!, y, fdcache)
            _run_monotonicity_check!(monotonicity_check, J, y)
            return J
        end
        jac = jac!
    end

    if _has_bounds(y̲, ȳ)
        # HJBVI bounds are mixed complementarity conditions, not box-constrained roots.
        if jac === nothing
            nlsolve_autodiff = autodiff == :finite ? :finiteforward : autodiff
            result = mcpsolve(G_helper!, y̲, ȳ, ypost; iterations = iterations, show_trace = verbose, ftol = maxdist, method = method, reformulation = reformulation, autoscale = autoscale, autodiff = nlsolve_autodiff)
        else
            df = OnceDifferentiable(G_helper!, jac, deepcopy(ypost), deepcopy(ypost), J0c)
            result = mcpsolve(df, y̲, ȳ, ypost; iterations = iterations, show_trace = verbose, ftol = maxdist, method = method, reformulation = reformulation, autoscale = autoscale)
        end
        return result.zero, result.residual_norm
    else
        G_solve! = (du, u, p) -> G_helper!(du, u)
        if jac === nothing
            f = NonlinearFunction{true}(G_solve!; jac_prototype = J0c)
        else
            f = NonlinearFunction{true}(G_solve!; jac = (J, u, p) -> jac(J, u), jac_prototype = J0c)
        end
        problem = NonlinearProblem(f, ypost)
        if autodiff == :forward
            autodiff = AutoForwardDiff()
        elseif autodiff == :finite
            autodiff = AutoFiniteDiff()
        elseif autodiff == :central
            autodiff = AutoFiniteDiff(fdjtype = Val(:central))
        else
            autodiff = autodiff
        end
        if method == :newton
            algorithm = NewtonRaphson(; autodiff = autodiff, concrete_jac = true)
        else
            algorithm = TrustRegion(; autodiff = autodiff, concrete_jac = true)
        end
        try
            result = solve(problem, algorithm; maxiters = iterations, abstol = maxdist, verbose = verbose, kwargs...)
            return result.u, norm(result.resid) / length(result.resid)
        catch err
            # catch SingularException error because can simply mean time step too big
            err isa SingularException || rethrow()
            return ypost, Inf
        end
    end
end

"""
    finiteschemesolve(F!, y0; kwargs...)

Lower-level nonlinear solver behind `pdesolve`. Finds `y` such that `F!(ydot, y)` writes the
residual into `ydot` and returns zero, using Newton's method with pseudo-transient
continuation from the initial guess `y0`.

`pdesolve` assembles `F!` (the finite-difference residual of the PDE) and calls this function,
so most users should call `pdesolve` instead. It accepts the same solver keyword arguments —
`Δ`, `iterations`, `method`, `maxdist`, `autodiff`, `y̲`, `ȳ`, … — and returns the tuple
`(y, residual_norm)`.
"""
function finiteschemesolve(G!, y0; Δ = 1.0, is_algebraic = fill(false, size(y0)...), iterations = 100, inner_iterations = 10, verbose = true, inner_verbose = false, method = :newton, autodiff = :forward, maxdist = sqrt(eps()), innerdist = sqrt(eps()), scale = 10.0, J0 = nothing, minΔ = 1e-9, y̲ = fill(-Inf, length(y0)), ȳ = fill(Inf, length(y0)), reformulation = :smooth, maxΔ = Inf, autoscale = true, monotonicity_check = nothing, kwargs...)
    ypost = y0
    ydot = zero(y0)
    # check that does not return NAN or zero
    G!(ydot, ypost)
    residual_norm = norm(ydot) / length(ydot)
    isnan(residual_norm) && throw(ArgumentError("G! returns NaN with the initial value"))
    if residual_norm <= maxdist
        verbose && @warn "G! already returns zero with the initial value"
        return ypost, residual_norm
    elseif Δ == Inf
        # infinite time step => solve in one step
        ypost, residual_norm = implicit_timestep(G!, y0, Δ; is_algebraic = is_algebraic, verbose = verbose, iterations = iterations,  method = method, autodiff = autodiff, maxdist = maxdist, J0 = J0, y̲ = y̲, ȳ = ȳ, monotonicity_check = monotonicity_check, kwargs...)
        return ypost, residual_norm
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
            y, nlresidual_norm = implicit_timestep(G!, ypost, Δ; is_algebraic = is_algebraic, verbose = inner_verbose, iterations = inner_iterations, method = method, autodiff = autodiff, maxdist = innerdist, J0 = J0, y̲ = y̲, ȳ = ȳ, reformulation = reformulation, monotonicity_check = monotonicity_check, kwargs...)
            G!(ydot, y)
            if _has_bounds(y̲, ȳ)
                # only unconstrained ydot is relevant for residual_norm calculation
                mask = y̲ .+ eps() .<= y .<= ȳ .- eps() 
                residual_norm, oldresidual_norm = norm(ydot .* mask) / sum(mask), residual_norm
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
        if verbose
            (iter >= iterations) && @warn "Algorithm did not converge: Iter higher than the limit $(iterations)"
            (Δ < minΔ) && @warn "Algorithm did not converge: TimeStep lower than the limit $(minΔ)"
        end
        return ypost, residual_norm
    end
end
