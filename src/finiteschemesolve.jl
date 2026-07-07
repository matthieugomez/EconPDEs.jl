##############################################################################
##
## Non Linear solver using Pseudo-Transient Continuation Method
##
##############################################################################
# Compile-once barrier between the model-specific residual closure and the solver stack.
# `finiteschemesolve`, `implicit_timestep`, FiniteDiff, and NonlinearSolve all
# specialize on the residual's type, so without the barrier that whole chain recompiles
# for every new model function (several hundred ms per model). Hiding the closure behind
# an untyped field gives the chain a single fixed type, compiled once and cached in the
# package image by the precompile workload. The price is one dynamic dispatch per
# residual sweep over the grid, which is negligible next to the sweep itself.
struct ResidualWrapper
    f::Any
end
(r::ResidualWrapper)(ydot, y) = (r.f(ydot, y); nothing)

# The barrier applies only when a sparsity pattern is available (1–3 states): the
# Jacobian is then colored finite differences, so the solver only ever calls the
# residual on Float64 vectors. Without a pattern (4+ states) NonlinearSolve
# differentiates through the residual with ForwardDiff, and the closure stays unwrapped
# so that path specializes on it as before.

# Build everything the colored sparse finite-difference Jacobian needs: the normalized
# sparse pattern (also used as the Jacobian prototype) and the FiniteDiff cache holding the
# coloring. The pattern is fixed across pseudo-transient iterations and implicit time steps,
# so callers compute this once and pass it to every `implicit_timestep` call.
# `pdesolve` passes the closed-form coloring of its stencil pattern (`stencil_colors`);
# without one, the coloring falls back to a greedy coloring of the pattern.
function _sparse_fd_setup(J0, y, autodiff, colorvec = nothing; force_diagonal = false)
    J0 === nothing && return nothing, nothing
    J0c = sparse(J0)
    if force_diagonal
        added_diagonal = false
        for i in 1:min(size(J0c)...)
            if iszero(J0c[i, i])
                J0c[i, i] = one(eltype(J0c))
                added_diagonal = true
            end
        end
        added_diagonal && (colorvec = nothing)
    end
    if colorvec === nothing
        colorvec = matrix_colors(J0c)
    else
        length(colorvec) == size(J0c, 2) || throw(ArgumentError("`colorvec` must assign a color to each column of `J0`: got $(length(colorvec)) colors for $(size(J0c, 2)) columns"))
    end
    fdtype = autodiff == :central ? Val(:central) : Val(:forward)
    fdcache = JacobianCache(y, fdtype, eltype(y); colorvec = colorvec, sparsity = J0c)
    return J0c, fdcache
end

function implicit_timestep(G!, ypost, Δ; is_algebraic = fill(false, size(ypost)...), iterations = 100, verbose = true, method = :newton, autodiff = :forward, maxdist = sqrt(eps()), J0 = nothing, colorvec = nothing, fdcache = nothing, lower_bound = fill(-Inf, length(ypost)), upper_bound = fill(Inf, length(ypost)), y̲ = nothing, ȳ = nothing, reformulation = :minmax, autoscale = true, monotonicity_check = nothing, kwargs...)
    method in (:newton, :trust_region) || throw(ArgumentError("method must be :newton or :trust_region"))
    # `y̲`/`ȳ` are deprecated aliases of `lower_bound`/`upper_bound`, kept so older scripts keep running
    y̲ === nothing || (Base.depwarn("the keyword `y̲` is deprecated; use `lower_bound` instead", :implicit_timestep); lower_bound = y̲)
    ȳ === nothing || (Base.depwarn("the keyword `ȳ` is deprecated; use `upper_bound` instead", :implicit_timestep); upper_bound = ȳ)
    # Any non-default bound turns the implicit solve into a mixed complementarity problem.
    has_bounds = any(x -> x != -Inf, lower_bound) || any(x -> x != Inf, upper_bound)
    function G_helper!(ydot, y)
        G!(ydot, y)
        ydot .-= .!is_algebraic .* (ypost .- y) ./ Δ
        return nothing
    end

    jac = nothing
    if fdcache === nothing
        J0c, fdcache = _sparse_fd_setup(J0, ypost, autodiff, colorvec; force_diagonal = has_bounds)
    else
        J0c = fdcache.sparsity
    end
    if fdcache !== nothing
        function jac!(J, y)
            finite_difference_jacobian!(J, G_helper!, y, fdcache)
            _try_run_monotonicity_check!(monotonicity_check, J, y, Δ, is_algebraic)
            return J
        end
        jac = jac!
    end

    solve_residual! = G_helper!
    solve_jac = jac
    if has_bounds
        if reformulation == :smooth
            @warn "`reformulation = :smooth` is no longer supported for bounded solves; using `:minmax` instead." maxlog = 1
        elseif reformulation != :minmax
            throw(ArgumentError("only `reformulation = :minmax` is supported for bounded solves"))
        end
        # HJBVI bounds are mixed complementarity conditions. The minmax residual is zero
        # exactly when F_i(y) = 0 in the interior, F_i(y) >= 0 at a lower bound, or
        # F_i(y) <= 0 at an upper bound.
        function mcp_residual!(ydot, y)
            G_helper!(ydot, y)
            @inbounds for i in eachindex(ydot, y, lower_bound, upper_bound)
                ydot[i] = min(max(ydot[i], y[i] - upper_bound[i]), y[i] - lower_bound[i])
            end
            return nothing
        end

        solve_residual! = mcp_residual!
        if fdcache !== nothing
            function mcp_jac!(J, y)
                if monotonicity_check !== nothing
                    finite_difference_jacobian!(J, G_helper!, y, fdcache)
                    _try_run_monotonicity_check!(monotonicity_check, J, y, Δ, is_algebraic)
                end
                finite_difference_jacobian!(J, mcp_residual!, y, fdcache)
                return J
            end
            solve_jac = mcp_jac!
        end
    end

    G_solve! = (du, u, p) -> solve_residual!(du, u)
    if solve_jac === nothing
        f = NonlinearFunction{true}(G_solve!; jac_prototype = J0c)
    else
        f = NonlinearFunction{true}(G_solve!; jac = (J, u, p) -> solve_jac(J, u), jac_prototype = J0c)
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

# Adaptive pseudo-transient time step, as described in the `finiteschemesolve` docstring
# under `scale`: an accepted step grows Δ by `scale * (old residual / new residual)`, with
# `scale` compounding across consecutive improving steps (`coef`). A rejected step shrinks
# Δ by 10 and caps the next step at half the failed value (`Δ_cap`) — otherwise Δ would
# regrow straight back to the value that just failed and the solve would cycle between
# growth and rejection. The cap relaxes by `scale` on every accepted step: whether a larger
# Δ will succeed depends on how close the iterate is to the solution, not on the residual
# level, so after a step or two a retry usually succeeds — and if it fails again the cap
# simply re-tightens. Holding the cap any longer pins the solve at a small fixed Δ, where
# each step costs a full inner solve but barely reduces the residual.
mutable struct TimeStepController
    Δ::Float64
    scale::Float64  # base growth factor for accepted steps
    maxΔ::Float64
    coef::Float64   # compounded growth factor across improving steps
    Δ_cap::Float64  # regrowth cap set by a rejected step, relaxed by each accepted step
end
TimeStepController(Δ, scale, maxΔ) = TimeStepController(Δ, scale, maxΔ, 1.0, maxΔ)

# the inner solve at Δ succeeded, moving the residual from `oldresidual_norm` to `residual_norm`
function accept!(stepper::TimeStepController, oldresidual_norm, residual_norm)
    stepper.coef = residual_norm <= oldresidual_norm ? stepper.scale * stepper.coef : 1.0
    stepper.Δ *= stepper.coef * oldresidual_norm / residual_norm
    if stepper.Δ > min(stepper.Δ_cap, stepper.maxΔ)
        # pinned at the cap: stop compounding the acceleration factor
        stepper.Δ = min(stepper.Δ_cap, stepper.maxΔ)
        stepper.coef = 1.0
    end
    stepper.Δ_cap = min(stepper.Δ_cap * stepper.scale, stepper.maxΔ)
    return stepper.Δ
end

# the inner solve at Δ failed (the residual is unchanged)
function reject!(stepper::TimeStepController)
    stepper.coef = 1.0
    stepper.Δ_cap = stepper.Δ / 2
    stepper.Δ /= 10
    return stepper.Δ
end

"""
    finiteschemesolve(G!, y0; kwargs...)

Lower-level nonlinear solver behind `pdesolve`. Finds `y` such that `G!(ydot, y)` writes the
residual into `ydot` and returns zero, using a nonlinear solver (Newton by default) with
pseudo-transient continuation from the initial guess `y0`.

`pdesolve` assembles `G!` (the finite-difference residual of the PDE) and calls this function,
so most users should call `pdesolve` instead. It returns the tuple `(y, residual_norm)`.

### Keyword arguments
* `Δ = 1.0`: initial pseudo-transient time step. `Δ = Inf` solves the stationary residual
  in one nonlinear solve, with no continuation.
* `scale = 10.0`: growth factor for the time step. After a successful step that reduces the
  residual, `Δ` is multiplied by `scale * (old residual / new residual)`, and `scale`
  compounds across consecutive improving steps. After a failed inner solve, `Δ` is divided
  by 10 and its regrowth is capped at half the failed `Δ`, so a `Δ` that just failed is not
  immediately retried; the cap relaxes by `scale` with each subsequent successful step.
* `minΔ = 1e-9`, `maxΔ = Inf`: bounds on the time step. The solve stops when `Δ < minΔ`
  (typically a sign that the scheme is non-monotone or the initial guess is poor).
* `iterations = 100`: maximum number of pseudo-transient (outer) iterations.
* `maxdist = sqrt(eps())`: convergence tolerance on the residual norm.
* `inner_iterations = 10`, `innerdist = sqrt(eps())`, `inner_verbose = false`: iteration
  limit, tolerance, and verbosity for each implicit time step (the inner solve).
* `method`: `:newton` (default) or `:trust_region`.
* `autodiff`: `:forward` (default), `:finite`, or `:central`. When a sparsity pattern `J0`
  is available, the Jacobian is computed by colored sparse finite differences
  (`:forward`/`:finite` forward differences, `:central` central differences); forward-mode
  AD is used only without a sparsity pattern.
* `is_algebraic`: `Bool` per entry of `y0`, marking algebraic (no time-derivative) equations.
* `J0 = nothing`: sparsity pattern of the Jacobian.
* `colorvec = nothing`: column coloring of `J0` for the colored finite differences
  (no two columns of the same color may share a nonzero row). Defaults to a greedy
  coloring of `J0`; pass one when it is known in closed form, as `pdesolve` does for
  its stencil pattern.
* `lower_bound`, `upper_bound`: lower/upper bounds on `y`. When any bound is finite, the
  problem is solved with the minmax mixed-complementarity reformulation. The Unicode
  keywords `y̲`/`ȳ` are deprecated aliases.
* `verbose = true`: print one line per pseudo-transient time step (`Iter`, `TimeStep`,
  `Residual`) and a final convergence summary. A ✗ in place of the residual means the
  inner nonlinear solve failed at that `Δ`: no step is taken (so there is no residual to
  show) and it is retried with `Δ/10` — the unusual causes, a `NaN` from the model or a
  singular Jacobian, are flagged in parentheses. With `verbose = false` a successful
  solve prints nothing; convergence failures are always reported with `@warn`.
"""
function finiteschemesolve(G!, y0; Δ = 1.0, is_algebraic = fill(false, size(y0)...), iterations = 100, inner_iterations = 10, verbose = true, inner_verbose = false, method = :newton, autodiff = :forward, maxdist = sqrt(eps()), innerdist = sqrt(eps()), scale = 10.0, J0 = nothing, colorvec = nothing, minΔ = 1e-9, lower_bound = fill(-Inf, length(y0)), upper_bound = fill(Inf, length(y0)), y̲ = nothing, ȳ = nothing, reformulation = :minmax, maxΔ = Inf, autoscale = true, monotonicity_check = nothing, kwargs...)
    method in (:newton, :trust_region) || throw(ArgumentError("method must be :newton or :trust_region"))
    # `y̲`/`ȳ` are deprecated aliases of `lower_bound`/`upper_bound`, kept so older scripts keep running
    y̲ === nothing || (Base.depwarn("the keyword `y̲` is deprecated; use `lower_bound` instead", :finiteschemesolve); lower_bound = y̲)
    ȳ === nothing || (Base.depwarn("the keyword `ȳ` is deprecated; use `upper_bound` instead", :finiteschemesolve); upper_bound = ȳ)
    # Bounds are static across pseudo-transient iterations; compute this once after alias normalization.
    has_bounds = any(x -> x != -Inf, lower_bound) || any(x -> x != Inf, upper_bound)
    tstart = time()
    ypost = y0
    ydot = zero(y0)
    # the sparsity pattern is fixed, so build the coloring and Jacobian cache once here
    # rather than inside every implicit time step
    J0c, fdcache = _sparse_fd_setup(J0, y0, autodiff, colorvec; force_diagonal = has_bounds)
    # check that does not return NAN or zero
    G!(ydot, ypost)
    residual_norm = norm(ydot) / length(ydot)
    isnan(residual_norm) && throw(ArgumentError("G! returns NaN with the initial value"))
    if residual_norm <= maxdist
        verbose && @warn "G! already returns zero with the initial value"
        return ypost, residual_norm
    elseif Δ == Inf
        # infinite time step => solve in one step
        ypost, residual_norm = implicit_timestep(G!, y0, Δ; is_algebraic = is_algebraic, verbose = verbose, iterations = iterations,  method = method, autodiff = autodiff, maxdist = maxdist, J0 = J0c, fdcache = fdcache, lower_bound = lower_bound, upper_bound = upper_bound, reformulation = reformulation, autoscale = autoscale, monotonicity_check = monotonicity_check, kwargs...)
        if residual_norm <= maxdist
            elapsed = time() - tstart
            # Format fast solves in milliseconds so they do not display as "0.0s".
            verbose && @printf "Converged (%s)\n" (elapsed < 1 ? @sprintf("%.0fms", 1000 * elapsed) : @sprintf("%.1fs", elapsed))
        else
            # warn even when verbose = false: solves embedded in silenced loops
            # (e.g. estimation) must still surface failures
            @warn "did not converge: residual $(@sprintf("%.2e", residual_norm)) > tolerance $(@sprintf("%.2e", maxdist)). `Δ = Inf` disables pseudo-transient continuation — try a finite initial time step `Δ`, a better initial guess, or more `iterations`."
        end
        return ypost, residual_norm
    else
        stepper = TimeStepController(Δ, scale, maxΔ)
        iter = 0
        if verbose
            @printf "Iter TimeStep Residual\n"
            @printf "---- -------- --------\n"
        end
        while (iter < iterations) && (stepper.Δ >= minΔ) && (residual_norm > maxdist)
            iter += 1
            y, nlresidual_norm = implicit_timestep(G!, ypost, stepper.Δ; is_algebraic = is_algebraic, verbose = inner_verbose, iterations = inner_iterations, method = method, autodiff = autodiff, maxdist = innerdist, J0 = J0c, fdcache = fdcache, lower_bound = lower_bound, upper_bound = upper_bound, reformulation = reformulation, autoscale = autoscale, monotonicity_check = monotonicity_check, kwargs...)
            G!(ydot, y)
            if has_bounds
                # only unconstrained ydot is relevant for residual_norm calculation
                mask = lower_bound .+ eps() .<= y .<= upper_bound .- eps()
                residual_norm, oldresidual_norm = norm(ydot .* mask) / sum(mask), residual_norm
            else
                residual_norm, oldresidual_norm = norm(ydot) / length(ydot), residual_norm
            end
            residual_norm = isnan(residual_norm) ? Inf : residual_norm
            if nlresidual_norm <= innerdist
                # the implicit time step is correctly solved: accept the point and grow Δ
                if verbose
                    @printf "%4d %8.2e %8.2e\n" iter stepper.Δ residual_norm
                end
                accept!(stepper, oldresidual_norm, residual_norm)
                ypost, y = y, ypost
            else
                if verbose
                    # ✗ sits in the residual column: the inner solve failed, so there is
                    # no new point to evaluate; only the unusual causes get a parenthetical
                    reason = isnan(nlresidual_norm) ? " (model returned NaN)" :
                             isinf(nlresidual_norm) ? " (singular Jacobian)" : ""
                    @printf "%4d %8.2e        ✗%s\n" iter stepper.Δ reason
                end
                # the implicit time step is not solved: revert and shrink Δ
                reject!(stepper)
                residual_norm = oldresidual_norm
            end
        end
        if residual_norm <= maxdist
            elapsed = time() - tstart
            # Format fast solves in milliseconds so they do not display as "0.0s".
            verbose && @printf "Converged after %d time steps (%s)\n" iter (elapsed < 1 ? @sprintf("%.0fms", 1000 * elapsed) : @sprintf("%.1fs", elapsed))
        elseif stepper.Δ < minΔ
            # warn even when verbose = false: solves embedded in silenced loops
            # (e.g. estimation) must still surface failures
            @warn "did not converge after $iter time steps: the pseudo-transient time step fell below `minΔ` = $minΔ (residual $(@sprintf("%.2e", residual_norm)) > tolerance $(@sprintf("%.2e", maxdist))). This usually means a non-monotone finite-difference scheme (wrong upwind direction) or a poor initial guess — try `check_monotonicity = true` or a better initial guess."
        else
            @warn "did not converge after $iter time steps: residual $(@sprintf("%.2e", residual_norm)) > tolerance $(@sprintf("%.2e", maxdist)). If the residual was still falling, increase `iterations`; otherwise try a better initial guess or `check_monotonicity = true`."
        end
        return ypost, residual_norm
    end
end
