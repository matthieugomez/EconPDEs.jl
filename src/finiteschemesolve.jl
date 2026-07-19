##############################################################################
##
## Non Linear solver using Pseudo-Transient Continuation Method
##
##############################################################################
# Compile-once barrier between the model-specific residual closure and the solver stack.
# `finiteschemesolve`, `implicit_timestep`, FiniteDiff, and the nonlinear solver all
# specialize on the residual's type, so without the barrier that whole chain recompiles
# for every new model function (several hundred ms per model). Hiding the closure behind
# an untyped field gives the chain a single fixed type, compiled once and cached in the
# package image by the precompile workload. The price is one dynamic dispatch per
# residual sweep over the grid, which is negligible next to the sweep itself.
struct ResidualWrapper
    f::Any
end
(r::ResidualWrapper)(ydot, y) = (r.f(ydot, y); nothing)

function _reject_removed_solver_keywords(kwargs)
    haskey(kwargs, :maxdist) && throw(ArgumentError("`maxdist` was removed in EconPDEs 2.0; use `abstol` instead."))
    haskey(kwargs, :innerdist) && throw(ArgumentError("`innerdist` was removed in EconPDEs 2.0; use `inner_abstol` instead."))
    haskey(kwargs, :y̲) && throw(ArgumentError("`y̲` was removed in EconPDEs 2.0; use `lower_bound` instead."))
    haskey(kwargs, :ȳ) && throw(ArgumentError("`ȳ` was removed in EconPDEs 2.0; use `upper_bound` instead."))
    haskey(kwargs, :method) && throw(ArgumentError("`method` was removed in EconPDEs 2.0; use `alg = :newton` or `alg = :trust_region`, or load NonlinearSolve and pass one of its algorithm objects."))
    haskey(kwargs, :autodiff) && throw(ArgumentError("`autodiff` was removed in EconPDEs 2.0; EconPDEs configures its Jacobian internally, or you can load NonlinearSolve and configure a NonlinearSolve algorithm through `alg`."))
    haskey(kwargs, :algorithm) && throw(ArgumentError("`algorithm` was renamed in EconPDEs 2.0; use `alg` instead."))
    haskey(kwargs, :iterations) && throw(ArgumentError("`iterations` was renamed in EconPDEs 2.0; use `maxiters` instead."))
    haskey(kwargs, :inner_iterations) && throw(ArgumentError("`inner_iterations` was renamed in EconPDEs 2.0; use `inner_maxiters` instead."))
    haskey(kwargs, :J0) && throw(ArgumentError("`J0` was renamed in EconPDEs 2.0; use `jac_prototype` instead."))
    haskey(kwargs, :reformulation) && throw(ArgumentError("`reformulation` was removed in EconPDEs 2.0; bounded solves always use the minmax mixed-complementarity residual."))
    haskey(kwargs, :autoscale) && throw(ArgumentError("`autoscale` was removed in EconPDEs 2.0."))
    return nothing
end

# The barrier applies only when a sparsity pattern is available: the Jacobian is then
# colored finite differences, so the solver only ever calls the residual on Float64
# vectors. Without a pattern, the backend differentiates through the residual, and the
# closure stays unwrapped so that path specializes on it as before.

# Normalize a user-supplied Jacobian prototype to the sparse matrix handed to the nonlinear
# solver. If the implicit step or MCP reformulation can add diagonal entries,
# make those structural nonzeros explicit before the nonlinear solver allocates a Jacobian.
function _normalize_jac_prototype(jac_prototype, colorvec = nothing; force_diagonal = false)
    jac_prototype === nothing && return nothing, colorvec
    J0c = sparse(jac_prototype)
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
    return J0c, colorvec
end

# Build everything the colored sparse finite-difference Jacobian needs: the normalized
# sparse pattern (also used as the Jacobian prototype) and the FiniteDiff cache holding the
# coloring. The pattern is fixed across pseudo-transient iterations and implicit time steps,
# so callers compute this once and pass it to every `implicit_timestep` call.
# `pdesolve` passes the closed-form coloring of its stencil pattern (`stencil_colors`);
# without one, the coloring falls back to a greedy coloring of the pattern.
function _sparse_fd_setup(jac_prototype, y, colorvec = nothing; force_diagonal = false)
    J0c, colorvec = _normalize_jac_prototype(jac_prototype, colorvec; force_diagonal = force_diagonal)
    J0c === nothing && return nothing, nothing
    if colorvec === nothing
        colorvec = matrix_colors(J0c)
    else
        length(colorvec) == size(J0c, 2) || throw(ArgumentError("`colorvec` must assign a color to each column of `jac_prototype`: got $(length(colorvec)) colors for $(size(J0c, 2)) columns"))
    end
    fdcache = JacobianCache(y, Val(:forward), eltype(y); colorvec = colorvec, sparsity = J0c)
    return J0c, fdcache
end

_has_dynamic_rows(is_algebraic::Bool) = !is_algebraic
_has_dynamic_rows(is_algebraic) = any(!, is_algebraic)

# A row is relaxed when it is dynamic and its pseudo-time step is finite; only relaxed
# rows receive the 1/Δ diagonal, so only they force structural diagonal entries.
_has_relaxed_rows(Δ::Real, is_algebraic) = isfinite(Δ) && _has_dynamic_rows(is_algebraic)
function _has_relaxed_rows(Δ::AbstractVector, is_algebraic)
    for i in eachindex(Δ)
        if is_algebraic isa Bool
            dynamic_row = !is_algebraic
        else
            dynamic_row = !is_algebraic[i]
        end
        dynamic_row && isfinite(Δ[i]) && return true
    end
    return false
end

function _add_implicit_diagonal!(J, Δ, is_algebraic)
    @inbounds for i in 1:min(size(J)...)
        if is_algebraic isa Bool
            dynamic_row = !is_algebraic
        else
            dynamic_row = !is_algebraic[i]
        end
        dynamic_row || continue
        Δi = Δ isa Real ? Δ : Δ[i]
        # ±Inf entries are undamped rows; a negative Δi is a forward-relaxed row and
        # contributes a negative diagonal, matching its negative residual diagonal
        isfinite(Δi) || continue
        J[i, i] += inv(Δi)
    end
    return J
end

function _apply_minmax_mcp_jacobian!(J, residual, y, lower_bound, upper_bound)
    @inbounds for i in eachindex(residual, y, lower_bound, upper_bound)
        lower_residual = y[i] - lower_bound[i]
        upper_residual = y[i] - upper_bound[i]
        if !(upper_residual <= residual[i] <= lower_residual)
            J[i, :] .= zero(eltype(J))
            J[i, i] = one(eltype(J))
        end
    end
    return J
end

# Symbols select NLsolve algorithms in core; package extensions add methods for algorithm
# objects from optional backends.
function _newtonstep(residual!, jac, jac_prototype, y0, alg::Symbol; maxiters, abstol, verbose, kwargs...)
    alg in (:newton, :trust_region) || throw(ArgumentError("NLsolve `alg` must be `:newton` or `:trust_region`; got `$(repr(alg))`"))
    try
        if jac === nothing
            result = nlsolve(residual!, y0; method = alg, autodiff = :forward,
                iterations = maxiters, ftol = abstol, show_trace = verbose, kwargs...)
        elseif jac_prototype === nothing
            result = nlsolve(residual!, jac, y0; method = alg,
                iterations = maxiters, ftol = abstol, show_trace = verbose, kwargs...)
        else
            df = OnceDifferentiable(residual!, jac, copy(y0), copy(y0), copy(jac_prototype))
            result = nlsolve(df, y0; method = alg, iterations = maxiters,
                ftol = abstol, show_trace = verbose, kwargs...)
        end
        final_residual = similar(result.zero)
        residual!(final_residual, result.zero)
        return result.zero, norm(final_residual) / length(final_residual)
    catch err
        # A singular Jacobian can simply mean that the pseudo-time step is too large.
        err isa SingularException || rethrow()
        return y0, Inf
    end
end

function _newtonstep(residual!, jac, jac_prototype, y0, alg; kwargs...)
    throw(ArgumentError("unsupported nonlinear solver algorithm `$(typeof(alg))`; load NonlinearSolve.jl before passing a NonlinearSolve algorithm object"))
end

function implicit_timestep(G!, ypost, Δ; is_algebraic = fill(false, size(ypost)...), lower_bound = fill(-Inf, length(ypost)), upper_bound = fill(Inf, length(ypost)), abstol = sqrt(eps()), verbose = true, alg = :newton, maxiters = 100, jac = nothing, jac_prototype = nothing, colorvec = nothing, fdcache = nothing, monotonicity_check = nothing, kwargs...)
    _reject_removed_solver_keywords(kwargs)
    # Any non-default bound turns the implicit solve into a mixed complementarity problem.
    has_bounds = any(x -> x != -Inf, lower_bound) || any(x -> x != Inf, upper_bound)
    force_diagonal = has_bounds || _has_relaxed_rows(Δ, is_algebraic)
    function G_helper!(ydot, y)
        G!(ydot, y)
        ydot .-= .!is_algebraic .* (ypost .- y) ./ Δ
        return nothing
    end

    implicit_jac = nothing
    if fdcache === nothing
        if jac === nothing
            J0c, fdcache = _sparse_fd_setup(jac_prototype, ypost, colorvec; force_diagonal = force_diagonal)
        else
            J0c, _ = _normalize_jac_prototype(jac_prototype, nothing; force_diagonal = force_diagonal)
        end
    else
        J0c = fdcache.sparsity
    end
    if jac !== nothing
        function analytic_implicit_jac!(J, y)
            jac(J, y)
            _add_implicit_diagonal!(J, Δ, is_algebraic)
            _try_run_monotonicity_check!(monotonicity_check, J, y, Δ, is_algebraic)
            return J
        end
        implicit_jac = analytic_implicit_jac!
    elseif fdcache !== nothing
        function fd_implicit_jac!(J, y)
            finite_difference_jacobian!(J, G_helper!, y, fdcache)
            _try_run_monotonicity_check!(monotonicity_check, J, y, Δ, is_algebraic)
            return J
        end
        implicit_jac = fd_implicit_jac!
    end

    solve_residual! = G_helper!
    solve_jac = implicit_jac
    if has_bounds
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
        if jac !== nothing
            base_residual = similar(ypost)
            function analytic_mcp_jac!(J, y)
                implicit_jac(J, y)
                G_helper!(base_residual, y)
                _apply_minmax_mcp_jacobian!(J, base_residual, y, lower_bound, upper_bound)
                return J
            end
            solve_jac = analytic_mcp_jac!
        elseif fdcache !== nothing
            function fd_mcp_jac!(J, y)
                if monotonicity_check !== nothing
                    finite_difference_jacobian!(J, G_helper!, y, fdcache)
                    _try_run_monotonicity_check!(monotonicity_check, J, y, Δ, is_algebraic)
                end
                finite_difference_jacobian!(J, mcp_residual!, y, fdcache)
                return J
            end
            solve_jac = fd_mcp_jac!
        end
    end

    return _newtonstep(solve_residual!, solve_jac, J0c, ypost, alg;
        maxiters = maxiters, abstol = abstol, verbose = verbose, kwargs...)
end

# Adaptive pseudo-transient time step, using a switched-evolution-relaxation (SER)
# residual ratio plus EconPDEs-specific acceleration and rejected-step caps. Keeping the
# update rule in this controller lets `finiteschemesolve` focus on solving/checking
# implicit steps.
mutable struct SERTimeStepController
    Δ::Float64
    scale::Float64  # base growth factor for accepted steps
    maxΔ::Float64
    acceleration::Float64 # compounded growth factor across improving steps
    Δ_cap::Float64        # regrowth cap set by a rejected step, relaxed by each accepted step
end
SERTimeStepController(Δ, scale, maxΔ) = SERTimeStepController(Δ, scale, maxΔ, 1.0, maxΔ)

# The inner solve at Δ succeeded, moving the residual from `oldresidual_norm` to
# `residual_norm`. Grow Δ by a SER-style residual ratio. Consecutive improving steps
# accelerate by `scale`; a step that worsens the residual resets that acceleration.
function accept_step!(stepper::SERTimeStepController, oldresidual_norm, residual_norm)
    if residual_norm <= oldresidual_norm
        stepper.acceleration *= stepper.scale
    else
        stepper.acceleration = 1.0
    end

    proposed_Δ = stepper.Δ * stepper.acceleration * oldresidual_norm / residual_norm
    Δ_limit = min(stepper.Δ_cap, stepper.maxΔ)
    if proposed_Δ > Δ_limit
        # Pinned at the cap: stop compounding the acceleration factor.
        stepper.Δ = Δ_limit
        stepper.acceleration = 1.0
    else
        stepper.Δ = proposed_Δ
    end
    stepper.Δ_cap = min(stepper.Δ_cap * stepper.scale, stepper.maxΔ)
    return stepper.Δ
end

# The inner solve at Δ failed, so the residual is unchanged. Shrink Δ and cap future
# regrowth at half the failed value; accepted steps relax that cap by `scale`.
function reject_step!(stepper::SERTimeStepController)
    stepper.acceleration = 1.0
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
* `is_algebraic`: `Bool` per entry of `y0`, marking algebraic (no time-derivative) equations.
* `lower_bound`, `upper_bound`: lower/upper bounds on `y`. When any bound is finite, the
  problem is solved with the minmax mixed-complementarity reformulation.
* `abstol = sqrt(eps())`: convergence tolerance on the residual norm.
* `verbose = true`: print one line per pseudo-transient time step (`Iter`, `TimeStep`,
  `Residual`) and a final convergence summary. A ✗ in place of the residual means the
  inner nonlinear solve failed at that `Δ`: no step is taken (so there is no residual to
  show) and it is retried with `Δ/10` — the unusual causes, a `NaN` from the model or a
  singular Jacobian, are flagged in parentheses. With `verbose = false` a successful
  solve prints nothing; convergence failures are still reported with `@warn`.
* `warn = true`: report convergence failures with `@warn` even when `verbose = false`,
  so solves embedded in silenced loops (e.g. estimation) still surface failures. Pass
  `warn = false` when the caller inspects convergence itself and a failure is an
  expected outcome (e.g. probing whether a guess is good enough for a one-shot solve).
* `alg = :newton`: NLsolve method used for each nonlinear solve; `:trust_region` is also
  supported. To use a NonlinearSolve.jl algorithm instead, load NonlinearSolve and pass
  an algorithm object, for example `alg = NonlinearSolve.TrustRegion()`.
* `maxiters = 100`: maximum number of pseudo-transient (outer) iterations.
* `Δ = 1.0`: initial pseudo-transient time step. `Δ = Inf` solves the stationary residual
  in one nonlinear solve, with no continuation. May also be a vector with one entry per
  element of `y0`. The sign of an entry sets the direction of that equation's false
  transient: positive marches it backward in time (the value-function convention, which
  requires a positive residual-Jacobian diagonal), negative forward (e.g. a market-clearing
  equation written in its natural tâtonnement direction `rt = r_implied - r`, whose
  diagonal is negative); `±Inf` entries disable damping for that equation. The entries fix
  a relative profile: the continuation scales all of them by one common adaptive factor.
* `scale = 10.0`: growth factor for the time step. After a successful step that reduces the
  residual, `Δ` is multiplied by `scale * (old residual / new residual)`, and `scale`
  compounds across consecutive improving steps. After a failed inner solve, `Δ` is divided
  by 10 and its regrowth is capped at half the failed `Δ`, so a `Δ` that just failed is not
  immediately retried; the cap relaxes by `scale` with each subsequent successful step.
* `minΔ = 1e-9`, `maxΔ = Inf`: bounds on the time step. The solve stops when `Δ < minΔ`
  (typically a sign that the scheme is non-monotone or the initial guess is poor). For a
  vector `Δ`, the bounds — like the `TimeStep` column of the verbose output — refer to
  the smallest finite `|Δ|` entry.
* `max_residual_growth = 10.0`: guard on accepted steps. A step whose inner solve
  converges but whose stationary residual exceeds `max_residual_growth` times the
  previous residual (a blowup, or a `NaN`, which counts as infinite growth) is reverted
  and retried with a smaller `Δ` instead of accepted. The default is deliberately loose:
  a false transient may legitimately increase the residual along the way, and rejecting
  every worsening step can deadlock the continuation at `minΔ`.
* `inner_maxiters = 10`, `inner_abstol = sqrt(eps())`, `inner_verbose = false`: iteration
  limit, tolerance, and verbosity for each implicit time step (the inner solve).
* `jac = nothing`: optional in-place Jacobian of the stationary residual `G!`, with
  signature `jac(J, y)`. EconPDEs adds the pseudo-time diagonal implied by `Δ` and
  `is_algebraic` before handing the Jacobian to the nonlinear solver.
* `jac_prototype = nothing`: prototype/sparsity pattern for the residual Jacobian. This is
  metadata for the nonlinear solver, not an algorithm choice. When `jac` is omitted,
  a sparse prototype makes EconPDEs compute the Jacobian by colored finite differences;
  when `jac` is supplied, the prototype controls the allocated Jacobian type.
* `colorvec = nothing`: column coloring of `jac_prototype` for the colored finite differences
  (no two columns of the same color may share a nonzero row). Defaults to a greedy
  coloring of `jac_prototype`; pass one when it is known in closed form, as `pdesolve`
  does for its stencil pattern.
"""
function finiteschemesolve(G!, y0; is_algebraic = fill(false, size(y0)...), lower_bound = fill(-Inf, length(y0)), upper_bound = fill(Inf, length(y0)), abstol = sqrt(eps()), verbose = true, warn = true, alg = :newton, maxiters = 100, Δ = 1.0, scale = 10.0, minΔ = 1e-9, maxΔ = Inf, max_residual_growth = 10.0, inner_maxiters = 10, inner_abstol = sqrt(eps()), inner_verbose = false, jac = nothing, jac_prototype = nothing, colorvec = nothing, monotonicity_check = nothing, kwargs...)
    _reject_removed_solver_keywords(kwargs)
    # Δ may be one pseudo-time step for all equations or one per equation. The sign of an
    # entry sets the direction of that equation's false transient: positive marches it
    # backward in time (the value-function convention), negative forward (e.g. tâtonnement
    # price adjustment written as `rt = r_implied - r`).
    if Δ isa AbstractVector
        length(Δ) == length(y0) || throw(ArgumentError("`Δ` must be a real number or have one entry per element of `y0`: got length $(length(Δ)) for $(length(y0)) unknowns"))
    elseif !(Δ isa Real)
        throw(ArgumentError("`Δ` must be a real number or a vector with one entry per element of `y0`; the per-unknown NamedTuple form is only accepted by `pdesolve`, which expands it over the grid"))
    end
    any(x -> isnan(x) || iszero(x), Δ) && throw(ArgumentError("`Δ` entries must be nonzero and not NaN (use `Inf` for no damping, a negative entry for forward relaxation)"))
    # Bounds are static across pseudo-transient iterations.
    has_bounds = any(x -> x != -Inf, lower_bound) || any(x -> x != Inf, upper_bound)
    tstart = time()
    ypost = y0
    ydot = zero(y0)
    # the sparsity pattern is fixed, so build the coloring and Jacobian cache once here
    # rather than inside every implicit time step
    force_diagonal = has_bounds || _has_relaxed_rows(Δ, is_algebraic)
    if jac === nothing
        J0c, fdcache = _sparse_fd_setup(jac_prototype, y0, colorvec; force_diagonal = force_diagonal)
    else
        J0c, _ = _normalize_jac_prototype(jac_prototype, nothing; force_diagonal = force_diagonal)
        fdcache = nothing
    end
    # check that does not return NAN or zero
    G!(ydot, ypost)
    residual_norm = norm(ydot) / length(ydot)
    isnan(residual_norm) && throw(ArgumentError("G! returns NaN with the initial value"))
    if residual_norm <= abstol
        verbose && @warn "G! already returns zero with the initial value"
        return ypost, residual_norm
    end
    if all(isinf, Δ)
        # infinite time step => solve in one step (the sign of an infinite entry is moot)
        ypost, residual_norm = implicit_timestep(G!, y0, Inf; is_algebraic = is_algebraic, lower_bound = lower_bound, upper_bound = upper_bound, abstol = abstol, verbose = verbose, alg = alg, maxiters = maxiters, jac = jac, jac_prototype = J0c, fdcache = fdcache, monotonicity_check = monotonicity_check, kwargs...)
        if residual_norm <= abstol
            elapsed = time() - tstart
            # Format fast solves in milliseconds so they do not display as "0.0s".
            verbose && @printf "Converged (%s)\n" (elapsed < 1 ? @sprintf("%.0fms", 1000 * elapsed) : @sprintf("%.1fs", elapsed))
        else
            # warn even when verbose = false (unless warn = false): solves embedded in
            # silenced loops (e.g. estimation) must still surface failures
            warn && @warn "did not converge: residual $(@sprintf("%.2e", residual_norm)) > tolerance $(@sprintf("%.2e", abstol)). `Δ = Inf` disables pseudo-transient continuation — try a finite initial time step `Δ`, a better initial guess, or more `maxiters`."
        end
        return ypost, residual_norm
    else
        # The controller adapts one positive multiplier over a fixed signed profile Δrel,
        # gauged so that stepper.Δ (the bounded, printed step) is the smallest finite |Δ|.
        if Δ isa AbstractVector
            Δbase = minimum(abs(x) for x in Δ if isfinite(x))
            Δrel = Δ ./ Δbase
        else
            Δbase = abs(Δ)
            Δrel = sign(Δ)
        end
        stepper = SERTimeStepController(Δbase, scale, maxΔ)
        iter = 0
        if verbose
            @printf "Iter TimeStep Residual\n"
            @printf "---- -------- --------\n"
        end
        while (iter < maxiters) && (stepper.Δ >= minΔ) && (residual_norm > abstol)
            iter += 1
            y, nlresidual_norm = implicit_timestep(G!, ypost, stepper.Δ .* Δrel; is_algebraic = is_algebraic, lower_bound = lower_bound, upper_bound = upper_bound, abstol = inner_abstol, verbose = inner_verbose, alg = alg, maxiters = inner_maxiters, jac = jac, jac_prototype = J0c, fdcache = fdcache, monotonicity_check = monotonicity_check, kwargs...)
            G!(ydot, y)
            if has_bounds
                # only unconstrained ydot is relevant for residual_norm calculation
                mask = lower_bound .+ eps() .<= y .<= upper_bound .- eps()
                residual_norm, oldresidual_norm = norm(ydot .* mask) / sum(mask), residual_norm
            else
                residual_norm, oldresidual_norm = norm(ydot) / length(ydot), residual_norm
            end
            residual_norm = isnan(residual_norm) ? Inf : residual_norm
            if nlresidual_norm <= inner_abstol && residual_norm <= max_residual_growth * oldresidual_norm
                # the implicit time step is correctly solved: accept the point and grow Δ
                if verbose
                    @printf "%4d %8.2e %8.2e\n" iter stepper.Δ residual_norm
                end
                accept_step!(stepper, oldresidual_norm, residual_norm)
                ypost, y = y, ypost
            elseif nlresidual_norm <= inner_abstol
                # the inner solve converged but the step grossly worsened the stationary
                # residual (a blowup or NaN, not the mild non-monotonicity a false
                # transient is allowed): revert the point and shrink Δ. Mild increases
                # must stay accepted — near a locally expansive transient (e.g. strong
                # price feedback), rejecting every worsening step deadlocks Δ at minΔ,
                # whereas riding out the drift lets the controller regrow Δ and escape.
                if verbose
                    growth = residual_norm / oldresidual_norm
                    growth_text = isfinite(growth) ? @sprintf("%.1f×", growth) : "∞"
                    @printf "%4d %8.2e        ✗ (residual grew %s > max_residual_growth)\n" iter stepper.Δ growth_text
                end
                reject_step!(stepper)
                residual_norm = oldresidual_norm
            else
                if verbose
                    # ✗ sits in the residual column: the inner solve failed, so there is
                    # no new point to evaluate; only the unusual causes get a parenthetical
                    reason = isnan(nlresidual_norm) ? " (model returned NaN)" :
                             isinf(nlresidual_norm) ? " (singular Jacobian)" : ""
                    @printf "%4d %8.2e        ✗%s\n" iter stepper.Δ reason
                end
                # the implicit time step is not solved: revert and shrink Δ
                reject_step!(stepper)
                residual_norm = oldresidual_norm
            end
        end
        if residual_norm <= abstol
            elapsed = time() - tstart
            # Format fast solves in milliseconds so they do not display as "0.0s".
            verbose && @printf "Converged after %d time steps (%s)\n" iter (elapsed < 1 ? @sprintf("%.0fms", 1000 * elapsed) : @sprintf("%.1fs", elapsed))
        elseif stepper.Δ < minΔ
            # warn even when verbose = false (unless warn = false): solves embedded in
            # silenced loops (e.g. estimation) must still surface failures
            warn && @warn "did not converge after $iter time steps: the pseudo-transient time step fell below `minΔ` = $minΔ (residual $(@sprintf("%.2e", residual_norm)) > tolerance $(@sprintf("%.2e", abstol))). This usually means a non-monotone finite-difference scheme (wrong upwind direction) or a poor initial guess — try `check_monotonicity = true` or a better initial guess."
        else
            warn && @warn "did not converge after $iter time steps: residual $(@sprintf("%.2e", residual_norm)) > tolerance $(@sprintf("%.2e", abstol)). If the residual was still falling, increase `maxiters`; otherwise try a better initial guess or `check_monotonicity = true`."
        end
        return ypost, residual_norm
    end
end
