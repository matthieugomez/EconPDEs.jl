##############################################################################
##
## Non Linear solver using Pseudo-Transient Continuation Method
##
##############################################################################

_has_bounds(lower_bound, upper_bound) = any(x -> x != -Inf, lower_bound) || any(x -> x != Inf, upper_bound)

function implicit_timestep(G!, ypost, Δ; is_algebraic = fill(false, size(ypost)...), iterations = 100, verbose = true, method = :newton, autodiff = :forward, maxdist = sqrt(eps()), J0 = nothing, lower_bound = fill(-Inf, length(ypost)), upper_bound = fill(Inf, length(ypost)), y̲ = nothing, ȳ = nothing, reformulation = :smooth, autoscale = true, monotonicity_check = nothing, kwargs...)
    method in (:newton, :trust_region) || throw(ArgumentError("method must be :newton or :trust_region"))
    # `y̲`/`ȳ` are deprecated aliases of `lower_bound`/`upper_bound`, kept so older scripts keep running
    y̲ === nothing || (Base.depwarn("the keyword `y̲` is deprecated; use `lower_bound` instead", :implicit_timestep); lower_bound = y̲)
    ȳ === nothing || (Base.depwarn("the keyword `ȳ` is deprecated; use `upper_bound` instead", :implicit_timestep); upper_bound = ȳ)
    G_helper!(ydot, y) = (G!(ydot, y) ; ydot .-= .!is_algebraic .* (ypost .- y) ./ Δ)

    jac = nothing
    J0c = J0 === nothing ? nothing : sparse(J0)
    if J0c !== nothing
        colorvec = matrix_colors(J0c)
        fdtype = autodiff == :central ? Val(:central) : Val(:forward)
        fdcache = JacobianCache(ypost, fdtype, eltype(ypost); colorvec = colorvec, sparsity = J0c)
        function jac!(J, y)
            finite_difference_jacobian!(J, G_helper!, y, fdcache)
            _try_run_monotonicity_check!(monotonicity_check, J, y, Δ, is_algebraic)
            return J
        end
        jac = jac!
    end

    if _has_bounds(lower_bound, upper_bound)
        # HJBVI bounds are mixed complementarity conditions, not box-constrained roots.
        if jac === nothing
            nlsolve_autodiff = autodiff == :finite ? :finiteforward : autodiff
            result = mcpsolve(G_helper!, lower_bound, upper_bound, ypost; iterations = iterations, show_trace = verbose, ftol = maxdist, method = method, reformulation = reformulation, autoscale = autoscale, autodiff = nlsolve_autodiff, kwargs...)
        else
            df = OnceDifferentiable(G_helper!, jac, deepcopy(ypost), deepcopy(ypost), J0c)
            result = mcpsolve(df, lower_bound, upper_bound, ypost; iterations = iterations, show_trace = verbose, ftol = maxdist, method = method, reformulation = reformulation, autoscale = autoscale, kwargs...)
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
so most users should call `pdesolve` instead. It returns the tuple `(y, residual_norm)`.

### Keyword arguments
* `Δ = 1.0`: initial pseudo-transient time step. `Δ = Inf` solves in a single Newton step.
* `scale = 10.0`: growth factor for the time step. After a successful step that reduces the
  residual, `Δ` is multiplied by `scale * (old residual / new residual)`; after a failed
  inner solve, `Δ` is divided by 10.
* `minΔ = 1e-9`, `maxΔ = Inf`: bounds on the time step. The solve stops when `Δ < minΔ`
  (typically a sign that the scheme is non-monotone or the initial guess is poor).
* `iterations = 100`: maximum number of pseudo-transient (outer) iterations.
* `maxdist = sqrt(eps())`: convergence tolerance on the residual norm.
* `inner_iterations = 10`, `innerdist = sqrt(eps())`, `inner_verbose = false`: Newton
  iteration limit, tolerance, and verbosity for each implicit time step (the inner solve).
* `method`: `:newton` (default) or `:trust_region`.
* `autodiff`: `:forward` (default), `:finite`, or `:central`. When a sparsity pattern `J0`
  is available, the Jacobian is computed by colored sparse finite differences
  (`:forward`/`:finite` forward differences, `:central` central differences); forward-mode
  AD is used only without a sparsity pattern.
* `is_algebraic`: `Bool` per entry of `y0`, marking algebraic (no time-derivative) equations.
* `J0 = nothing`: sparsity pattern of the Jacobian.
* `lower_bound`, `upper_bound`: lower/upper bounds on `y`. When any bound is finite, the
  problem is solved as a mixed complementarity problem with `NLsolve.mcpsolve`, with
  `reformulation` (`:smooth`, the default, or `:minmax`) and `autoscale = true` passed
  through to it. The Unicode keywords `y̲`/`ȳ` are deprecated aliases.
* `verbose = true`: print outer-iteration progress (`Iter`, `TimeStep`, `Residual`). A `NaN`
  residual line means the inner Newton solve failed and the time step was reduced.
"""
function finiteschemesolve(G!, y0; Δ = 1.0, is_algebraic = fill(false, size(y0)...), iterations = 100, inner_iterations = 10, verbose = true, inner_verbose = false, method = :newton, autodiff = :forward, maxdist = sqrt(eps()), innerdist = sqrt(eps()), scale = 10.0, J0 = nothing, minΔ = 1e-9, lower_bound = fill(-Inf, length(y0)), upper_bound = fill(Inf, length(y0)), y̲ = nothing, ȳ = nothing, reformulation = :smooth, maxΔ = Inf, autoscale = true, monotonicity_check = nothing, kwargs...)
    method in (:newton, :trust_region) || throw(ArgumentError("method must be :newton or :trust_region"))
    # `y̲`/`ȳ` are deprecated aliases of `lower_bound`/`upper_bound`, kept so older scripts keep running
    y̲ === nothing || (Base.depwarn("the keyword `y̲` is deprecated; use `lower_bound` instead", :finiteschemesolve); lower_bound = y̲)
    ȳ === nothing || (Base.depwarn("the keyword `ȳ` is deprecated; use `upper_bound` instead", :finiteschemesolve); upper_bound = ȳ)
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
        ypost, residual_norm = implicit_timestep(G!, y0, Δ; is_algebraic = is_algebraic, verbose = verbose, iterations = iterations,  method = method, autodiff = autodiff, maxdist = maxdist, J0 = J0, lower_bound = lower_bound, upper_bound = upper_bound, monotonicity_check = monotonicity_check, kwargs...)
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
            y, nlresidual_norm = implicit_timestep(G!, ypost, Δ; is_algebraic = is_algebraic, verbose = inner_verbose, iterations = inner_iterations, method = method, autodiff = autodiff, maxdist = innerdist, J0 = J0, lower_bound = lower_bound, upper_bound = upper_bound, reformulation = reformulation, monotonicity_check = monotonicity_check, kwargs...)
            G!(ydot, y)
            if _has_bounds(lower_bound, upper_bound)
                # only unconstrained ydot is relevant for residual_norm calculation
                mask = lower_bound .+ eps() .<= y .<= upper_bound .- eps() 
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
