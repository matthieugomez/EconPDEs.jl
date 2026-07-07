function _reject_pdesolve_lower_level_keywords(kwargs)
    haskey(kwargs, :jac) && throw(ArgumentError("`jac` is only supported by `finiteschemesolve`; `pdesolve` constructs the discretized residual internally."))
    haskey(kwargs, :jac_prototype) && throw(ArgumentError("`jac_prototype` is only supported by `finiteschemesolve`; `pdesolve` builds the stencil sparsity pattern internally."))
    haskey(kwargs, :colorvec) && throw(ArgumentError("`colorvec` is only supported by `finiteschemesolve`; `pdesolve` computes the stencil coloring internally."))
    haskey(kwargs, :fdcache) && throw(ArgumentError("`fdcache` is internal and cannot be passed to `pdesolve`."))
    haskey(kwargs, :monotonicity_check) && throw(ArgumentError("`monotonicity_check` is internal; use `check_monotonicity` instead."))
    return nothing
end

function _reject_time_dependent_stationary_keywords(kwargs)
    for keyword in (:Δ, :scale, :minΔ, :maxΔ)
        haskey(kwargs, keyword) && throw(ArgumentError("`$keyword` only applies to stationary pseudo-transient continuation; time-dependent `pdesolve(..., τs)` uses the spacing of `τs` as the time step."))
    end
    return nothing
end

"""
    pdesolve(pde, grid, guess [, τs]; kwargs...)

Solve a system of nonlinear ODEs/PDEs — typically Hamilton–Jacobi–Bellman equations —
by finite differences. Handles an arbitrary number of coupled unknown functions on a
state grid with any positive number of state variables.

### Positional arguments
* `pde`: the local equation, a function `(state, u) -> out` (or `(state, u, t) -> out` for a
  time-dependent equation), where
    - `state` is a `NamedTuple` with the current grid point (one entry per state variable),
    - `u` is a `NamedTuple` with each unknown function and its finite-difference derivatives
      at that point (e.g. `v`, `vk_up`, `vk_down`, `vkk`, …),
    - `out` is a `NamedTuple` with one time derivative per unknown (e.g. `(; vt)`).
  `pde` may also return a second `NamedTuple` of objects to save on the grid; these are
  returned in `result.saved`.
* `grid`: a `NamedTuple` (or `OrderedDict`) mapping each state variable name to an
  `AbstractVector` (its grid), e.g. `(; k = range(...))`.
* `guess`: a `NamedTuple` (or `OrderedDict`) mapping each unknown name to an array of initial
  values with the same shape as the grid. For a time-dependent problem, it is the terminal
  value at `τs[end]`.
* `τs` (optional): an increasing time grid. The equation is then solved backward from
  `τs[end]`, and each solution array gains a trailing time dimension:
  `result.solution.v[.., i]` is the solution at time `τs[i]`.

### Keyword arguments
* `bc`: boundary derivatives, a `NamedTuple` (or `OrderedDict`) mapping
  `Symbol(unknown, state)` to a `(lower, upper)` tuple (scalars, or arrays matching the
  boundary slice). Entries may be given for any subset of (unknown, state) pairs; omitted
  pairs default to reflecting boundaries (zero first derivative at the boundary). Both
  entries are `∂v/∂state` oriented in the increasing-state direction.
* `is_algebraic`: a `NamedTuple` (or `OrderedDict`) of `Bool`s marking equations that are
  algebraic (no time derivative) rather than PDEs. Defaults to all `false`.
* `lower_bound`, `upper_bound`: lower/upper bounds on the unknowns for HJB variational
  inequalities (optimal stopping), solved as a mixed complementarity problem. Pass either
  arrays in the same order as the unknowns or a `NamedTuple`/`OrderedDict` with the same
  names as `guess`. Default to `-Inf`/`Inf` (unbounded).
* `alg = NonlinearSolve.NewtonRaphson()`: NonlinearSolve algorithm used for each
  nonlinear solve. EconPDEs exports the `NonlinearSolve` module, so pass any compatible
  algorithm object, for example `alg = NonlinearSolve.TrustRegion()`.
* `abstol`: convergence tolerance on the residual. Defaults to `sqrt(eps())`.
* `maxiters`: maximum number of pseudo-transient iterations. Defaults to 100.
* `Δ`: initial pseudo-transient time step for stationary problems. Defaults to `1.0`;
  pass `Δ = Inf` to solve the stationary residual in one nonlinear solve, with no
  continuation. Not accepted when a time grid `τs` is supplied; then the spacing of `τs`
  is the time step.
* `verbose`: print a one-line problem summary, convergence progress, and a final
  convergence summary. Defaults to `true`. With `verbose = false` a successful solve
  prints nothing — convergence failures are always reported with `@warn` — so `pdesolve`
  can run inside a loop (e.g. an estimation) without flooding the log.
* Further stationary solver knobs (`scale`, `minΔ`, `maxΔ`) are forwarded to
  [`finiteschemesolve`](@ref). Further inner-solve knobs (`inner_maxiters`,
  `inner_abstol`, `inner_verbose`) are also accepted for time-dependent problems.
* `check_monotonicity`: if true, warn when the assembled residual Jacobian has same-variable
  spatial off-diagonal entries with the wrong monotonicity sign. Defaults to false. A flipped
  sign convention — returning `RHS - ρv` instead of `-(RHS - ρv)`, which makes the Jacobian
  diagonal negative at every grid point — is detected and warned about regardless of this
  option (whenever `Δ` is finite).
* `monotonicity_tol`: tolerance for the monotonicity check. Defaults to 1e-6.
* `monotonicity_max_warnings`: maximum number of monotonicity warnings to print. Defaults to 5.

Returns an `EconPDEResult` with fields `solution` (a `NamedTuple` with the solved unknowns),
`residual_norm`, and `saved` (a `NamedTuple` with the saved objects, together with the solved
unknowns, or an empty `NamedTuple` if the PDE saves nothing). `result.converged` reports
whether the residual met the tolerance (across all times, for a time-dependent problem).
"""
function pdesolve(pde, @nospecialize(grid), @nospecialize(guess), τs::Union{Nothing, AbstractVector} = nothing; is_algebraic = nothing, bc = nothing, lower_bound = nothing, upper_bound = nothing, abstol = sqrt(eps()), verbose = true, alg = NonlinearSolve.NewtonRaphson(), maxiters = 100, check_monotonicity = false, monotonicity_tol = 1e-6, monotonicity_max_warnings = 5, kwargs...)
    _reject_removed_solver_keywords(kwargs)
    _reject_pdesolve_lower_level_keywords(kwargs)
    # `grid`, `guess`, `is_algebraic`, and `bc` may be passed either as an OrderedDict or as a
    # NamedTuple. Normalize to a NamedTuple so everything below runs on one uniform representation.
    grid = _asnamedtuple(grid)
    guess = _asnamedtuple(guess)
    _check_generated_names(collect(keys(grid)), collect(keys(guess)))
    stategrid = StateGrid(grid)
    S = size(stategrid)
    for (name, v) in pairs(guess)
        size(v) == S || throw(ArgumentError("the initial guess (e.g. terminal value) for `$name` has size $(size(v)) but the state grid has size $S"))
    end
    is_algebraic = _fill_is_algebraic(is_algebraic, guess, S)
    # Concatenate the per-unknown arrays into one array with a trailing unknown dimension.
    guess_array = cat(values(guess)...; dims = ndims(first(values(guess))) + 1)
    is_algebraic_array = cat(values(is_algebraic)...; dims = ndims(first(values(is_algebraic))) + 1)
    bc_array = _bc_array(bc, guess, grid)
    lower_bound_array = _bound_array(lower_bound, guess, -Inf, :lower_bound)
    upper_bound_array = _bound_array(upper_bound, guess, Inf, :upper_bound)
    # create sparsity pattern and its coloring (closed-form, from the stencil structure)
    jac_prototype = sparse_jacobian(stategrid, guess)
    colorvec = jac_prototype === nothing ? nothing : stencil_colors(S, length(guess))
    # the sign-convention check runs by default; the stencil (monotonicity) warnings are opt-in
    monotonicity_check = MonotonicityChecker(stategrid, guess; tol = monotonicity_tol, max_warnings = monotonicity_max_warnings, check_stencils = check_monotonicity)
    if τs === nothing
        return _solve_stationary(pde, stategrid, guess, guess_array, is_algebraic_array, bc_array, lower_bound_array, upper_bound_array, jac_prototype, colorvec, monotonicity_check; abstol = abstol, verbose = verbose, alg = alg, maxiters = maxiters, kwargs...)
    else
        _reject_time_dependent_stationary_keywords(kwargs)
        _check_time_grid(τs)
        return _solve_backward(pde, stategrid, guess, τs, guess_array, is_algebraic_array, bc_array, lower_bound_array, upper_bound_array, jac_prototype, colorvec, monotonicity_check; abstol = abstol, verbose = verbose, alg = alg, maxiters = maxiters, kwargs...)
    end
end

# Stationary problem: one pseudo-transient continuation solve of the stationary residual.
function _solve_stationary(pde, stategrid::StateGrid, @nospecialize(guess), guess_array, is_algebraic_array, bc_array, lower_bound_array, upper_bound_array, jac_prototype, colorvec, monotonicity_check; abstol, verbose, alg, maxiters, kwargs...)
    ynames = tuple(keys(guess)...)
    solutionnames = Val(ynames)
    saved = _init_saved(pde, stategrid, solutionnames, ynames, guess_array, bc_array)
    verbose && println("Solving for ", _problem_description(guess, size(stategrid)))
    residual! = (ydot, yvec) -> pde!(pde, stategrid, solutionnames, ydot, yvec, bc_array, size(guess_array))
    G! = jac_prototype === nothing ? residual! : ResidualWrapper(residual!)
    y_array, residual_norm = finiteschemesolve(G!, vec(guess_array); is_algebraic = vec(is_algebraic_array), jac_prototype = jac_prototype, colorvec = colorvec, lower_bound = vec(lower_bound_array), upper_bound = vec(upper_bound_array), abstol = abstol, verbose = verbose, alg = alg, maxiters = maxiters, monotonicity_check = monotonicity_check, kwargs...)
    y_array = reshape(y_array, size(guess_array)...)
    # split the concatenated solution array back into one named array per unknown
    solution = NamedTuple{ynames}(ntuple(k -> collect(selectdim(y_array, ndims(y_array), k)), length(ynames)))
    if saved === nothing
        saved = NamedTuple()
    else
        _fill_saved!(saved, pde, stategrid, solutionnames, y_array, bc_array)
        saved = merge(solution, saved)
    end
    return EconPDEResult(solution, residual_norm, saved, abstol)
end

# Time-dependent problem: starting from the terminal condition at τs[end], take one
# implicit time step per interval of the time grid, recording the solution at each time.
function _solve_backward(pde, stategrid::StateGrid, @nospecialize(guess), τs, guess_array, is_algebraic_array, bc_array, lower_bound_array, upper_bound_array, jac_prototype, colorvec, monotonicity_check; abstol, verbose, alg, maxiters = 100, inner_maxiters = maxiters, inner_abstol = abstol, inner_verbose = false, kwargs...)
    ynames = tuple(keys(guess)...)
    solutionnames = Val(ynames)
    S = size(stategrid)
    N = length(S)
    # one array per unknown, with a trailing time dimension: solution.v[.., i] is the value at τs[i]
    solution = NamedTuple{ynames}(ntuple(_ -> Array{Float64}(undef, S..., length(τs)), length(ynames)))
    # the PDE function may or may not take a time argument; check once, not at every time step
    has_time = _pde_has_time_argument(pde, stategrid, solutionnames, guess_array, bc_array, τs[end])
    pde_terminal = has_time ? ((state, u) -> pde(state, u, τs[end])) : pde
    # scratch arrays on the state grid, refilled at each time and copied into the time slices
    saved_scratch = _init_saved(pde_terminal, stategrid, solutionnames, ynames, guess_array, bc_array)
    if saved_scratch === nothing
        saved = NamedTuple()
    else
        saved = merge(solution, map(_ -> Array{Float64}(undef, S..., length(τs)), saved_scratch))
    end
    y_array = guess_array
    residual_norms = zeros(length(τs))
    # the sparsity pattern is fixed, so build the coloring and Jacobian cache once
    # rather than inside every backward time step
    J0c, fdcache = _sparse_fd_setup(jac_prototype, vec(guess_array), colorvec)
    tstart = time()
    if verbose
        @printf "Solving for %s, backward from τ = %g to %g (%d steps)\n" _problem_description(guess, size(stategrid)) τs[end] τs[1] (length(τs) - 1)
        @printf "    Time Residual\n"
        @printf "-------- --------\n"
    end
    for iτ in length(τs):(-1):1
        for (k, name) in enumerate(ynames)
            copyto!(selectdim(solution[name], N + 1, iτ), selectdim(y_array, N + 1, k))
        end
        pde_onestep = has_time ? ((state, u) -> pde(state, u, τs[iτ])) : pde
        if saved_scratch !== nothing
            _fill_saved!(saved_scratch, pde_onestep, stategrid, solutionnames, y_array, bc_array)
            for (name, v) in pairs(saved_scratch)
                copyto!(selectdim(saved[name], N + 1, iτ), v)
            end
        end
        if iτ > 1
            residual! = (ydot, yvec) -> pde!(pde_onestep, stategrid, solutionnames, ydot, yvec, bc_array, size(guess_array))
            G! = J0c === nothing ? residual! : ResidualWrapper(residual!)
            y_array, residual_norms[iτ] = implicit_timestep(G!, vec(y_array), τs[iτ] - τs[iτ-1]; is_algebraic = vec(is_algebraic_array), verbose = inner_verbose, maxiters = inner_maxiters, alg = alg, jac_prototype = J0c, fdcache = fdcache, lower_bound = vec(lower_bound_array), upper_bound = vec(upper_bound_array), abstol = inner_abstol, monotonicity_check = monotonicity_check, kwargs...)
            if verbose
                @printf "%8g %8.2e\n" τs[iτ-1] residual_norms[iτ]
            end
            # an unconverged step would otherwise be silently accepted and propagated
            # to all earlier times; `!(x <= tol)` also catches a NaN residual
            if !(residual_norms[iτ] <= abstol)
                @warn "the implicit time step at τ = $(τs[iτ-1]) did not converge (residual norm: $(@sprintf("%.2e", residual_norms[iτ]))); the solution at this and earlier times may be inaccurate — try a finer time grid or more `inner_maxiters`"
            end
            y_array = reshape(y_array, size(guess_array)...)
        end
    end
    if verbose
        nfailed = count(iτ -> !(residual_norms[iτ] <= abstol), 2:length(τs))
        elapsed = time() - tstart
        # Format fast solves in milliseconds so they do not display as "0.0s".
        elapsed_text = elapsed < 1 ? @sprintf("%.0fms", 1000 * elapsed) : @sprintf("%.1fs", elapsed)
        if nfailed == 0
            @printf "Completed %d time steps (%s)\n" (length(τs) - 1) elapsed_text
        else
            @printf "Completed %d time steps (%s): %d steps did not converge\n" (length(τs) - 1) elapsed_text nfailed
        end
    end
    return EconPDEResult(solution, residual_norms, saved, abstol)
end

# e.g. "2 unknowns (pA, pB) on a 30×40 grid" — shared by the preambles of both solve paths
function _problem_description(@nospecialize(guess), S)
    n = length(guess)
    string(n, n == 1 ? " unknown (" : " unknowns (", join(keys(guess), ", "), ") on a ",
           length(S) == 1 ? "$(S[1])-point" : join(S, "×"), " grid")
end

function _asnamedtuple(x)
    x isa NamedTuple && return x
    # A plain Dict is rejected: its iteration order is arbitrary, and the order of names
    # determines how solution arrays are laid out.
    x isa Dict && throw(ArgumentError("pass a NamedTuple, e.g. `(; k = ...)`, or an OrderedDict — a Dict has no fixed iteration order, and the order of the names determines how arrays are laid out"))
    return NamedTuple(x)
end

function _check_time_grid(τs)
    isempty(τs) && throw(ArgumentError("the time grid `τs` must contain at least one time"))
    previous = first(τs)
    (previous isa Real && isfinite(previous)) || throw(ArgumentError("the time grid `τs` must contain finite real numbers"))
    for t in Iterators.drop(τs, 1)
        (t isa Real && isfinite(t)) || throw(ArgumentError("the time grid `τs` must contain finite real numbers"))
        t > previous || throw(ArgumentError("the time grid `τs` must be strictly increasing"))
        previous = t
    end
    return nothing
end

function _pde_has_time_argument(pde, stategrid::StateGrid, solutionnames, y_array::AbstractArray, bc_array::AbstractArray, t)
    i0 = first(eachindex(stategrid))
    state = stategrid[i0]
    derivatives = differentiate(solutionnames, stategrid, y_array, i0, bc_array)
    applicable(pde, state, derivatives, t) && return true
    applicable(pde, state, derivatives) && return false
    throw(ArgumentError("the PDE function must accept either `(state, u)` or `(state, u, t)` at the supplied grid and time types"))
end

# Derivative fields are named by concatenation (unknown * state * suffix), so some
# combinations of names generate the same field twice — e.g. states `(k, a)` with
# unknowns `(v, vk)` both generate `vka_up`. Without this check, the collision surfaces
# as a lowering error inside the @generated `differentiate` with no hint that the
# user's naming choice is the cause.
function _check_generated_names(statenames::Vector{Symbol}, ynames::Vector{Symbol})
    seen = Dict{Symbol, String}()
    function register!(name::Symbol, desc::String)
        if haskey(seen, name)
            throw(ArgumentError("ambiguous names: the field `$name` would mean both $(seen[name]) and $desc. Derivative fields are named by concatenating the unknown and state names (e.g. `v` and `k` give `vk_up`); rename one of the state variables or unknown functions so the generated names are unambiguous."))
        end
        seen[name] = desc
    end
    for s in ynames
        register!(s, "the value of unknown `$s`")
        for x in statenames
            register!(Symbol(s, x, :_up), "the forward derivative of `$s` in `$x`")
            register!(Symbol(s, x, :_down), "the backward derivative of `$s` in `$x`")
            register!(Symbol(s, x, x), "the second derivative of `$s` in `$x`")
        end
        for (i1, x1) in enumerate(statenames)
            for x2 in statenames[(i1 + 1):end]
                register!(Symbol(s, x1, x2), "the cross derivative of `$s` in `$x1`, `$x2`")
                register!(Symbol(s, x1, x2, :_up), "the directional (main-diagonal) cross derivative of `$s` in `$x1`, `$x2`")
                register!(Symbol(s, x1, x2, :_down), "the directional (anti-diagonal) cross derivative of `$s` in `$x1`, `$x2`")
            end
        end
    end
    return nothing
end

function _fill_is_algebraic(is_algebraic, guess, S)
    # Defaults to all-false, one Bool per grid point, with the same names/order as `guess`.
    is_algebraic === nothing && return map(_ -> fill(false, S), guess)
    ia = _asnamedtuple(is_algebraic)
    issetequal(keys(ia), keys(guess)) || throw(ArgumentError("`is_algebraic` must have the same names as the initial guess: got $(collect(keys(ia))), expected $(collect(keys(guess)))"))
    # reorder to match the guess, so the flags line up with the right unknowns
    NamedTuple{keys(guess)}(map(n -> _algebraic_array(ia[n], S, n), keys(guess)))
end

function _algebraic_array(value, S, name)
    value isa Bool && return fill(value, S)
    if value isa AbstractArray
        size(value) == S || throw(ArgumentError("`is_algebraic.$name` has size $(size(value)) but the state grid has size $S"))
        all(x -> x isa Bool, value) || throw(ArgumentError("`is_algebraic.$name` must be a Bool or an array of Bools"))
        return collect(value)
    end
    throw(ArgumentError("`is_algebraic.$name` must be a Bool or an array of Bools"))
end

function _bc_array(bc, guess, grid)
    # Default reflecting boundaries: zero derivative for every (unknown, state) pair.
    bc === nothing && return zeros(size(first(values(guess)))..., length(guess))
    bc = _asnamedtuple(bc)
    keys_grid = collect(keys(grid))
    valid = [Symbol(yname, s) for yname in keys(guess) for s in keys_grid]
    for key in keys(bc)
        key in valid || throw(ArgumentError("unknown `bc` entry `$key`: valid entries are $(valid), i.e. Symbol(unknown, state)"))
    end
    out = zeros(size(first(values(guess)))..., length(guess))
    for (k, yname) in enumerate(keys(guess))
        bck = selectdim(out, ndims(out), k)
        for (d, s) in enumerate(keys_grid)
            key = Symbol(yname, s)
            # entries may be given for any subset of (unknown, state) pairs;
            # omitted pairs keep the reflecting default (zero boundary derivative)
            haskey(bc, key) || continue
            lo, hi = bc[key]
            selectdim(bck, d, 1) .= lo
            selectdim(bck, d, size(bck, d)) .= hi
        end
    end
    return out
end

function _bound_array(bound, guess, default, keyword::Symbol)
    S = size(first(values(guess)))
    F = length(guess)
    if bound === nothing
        return fill(default, S..., F)
    elseif bound isa Number
        return fill(Float64(bound), S..., F)
    elseif bound isa NamedTuple || bound isa AbstractDict
        b = _asnamedtuple(bound)
        issetequal(keys(b), keys(guess)) || throw(ArgumentError("`$keyword` must have the same names as the initial guess: got $(collect(keys(b))), expected $(collect(keys(guess)))"))
        arrays = NamedTuple{keys(guess)}(map(n -> _bound_component(b[n], S, keyword, n), keys(guess)))
        return cat(values(arrays)...; dims = length(S) + 1)
    elseif bound isa AbstractArray
        if size(bound) == (S..., F)
            return Float64.(bound)
        elseif F == 1 && size(bound) == S
            return reshape(Float64.(bound), S..., 1)
        elseif length(bound) == prod(S) * F
            return reshape(Float64.(vec(bound)), S..., F)
        else
            throw(ArgumentError("`$keyword` has size $(size(bound)) but must have size $S for one unknown, size $(tuple(S..., F)) for all unknowns, or length $(prod(S) * F) when flattened"))
        end
    end
    throw(ArgumentError("`$keyword` must be a scalar, an array, or a NamedTuple/OrderedDict with the same names as the initial guess"))
end

function _bound_component(value, S, keyword::Symbol, name)
    value isa Number && return fill(Float64(value), S)
    if value isa AbstractArray
        size(value) == S || throw(ArgumentError("`$keyword.$name` has size $(size(value)) but the state grid has size $S"))
        return Float64.(value)
    end
    throw(ArgumentError("`$keyword.$name` must be a scalar or an array with the same size as the state grid"))
end

# Call the PDE function once at the first grid point to learn whether it returns a second
# NamedTuple of objects to save on the grid; if so, allocate one array per saved object.
function _init_saved(pde, stategrid::StateGrid, solutionnames, ynames, y_array::AbstractArray, bc_array::AbstractArray)
    i0 = first(eachindex(stategrid))
    derivatives = differentiate(solutionnames, stategrid, y_array, i0, bc_array)
    result = pde(stategrid[i0], derivatives)
    residual, optional = _split_pde_output(result)
    _check_residual_names(residual, ynames)
    optional === nothing && return nothing
    return map(_ -> Array{Float64}(undef, size(stategrid)), optional)
end

function _split_pde_output(result)
    # Single-return models give residuals; two-return models give residuals plus saved objects.
    result isa NamedTuple && return result, nothing
    result isa Tuple{<:NamedTuple, <:NamedTuple} && return result
    throw(ArgumentError("the PDE function must return a NamedTuple of time derivatives, e.g. `(; vt)`, optionally followed by a second NamedTuple of objects to save on the grid — got a value of type $(typeof(result))"))
end

function _check_residual_names(residual::NamedTuple, ynames)
    expected = map(n -> Symbol(n, :t), ynames)
    missing_names = setdiff(expected, keys(residual))
    isempty(missing_names) || throw(ArgumentError("the PDE function must return one time derivative per unknown, named `Symbol(unknown, :t)`: expected $(collect(expected)), got $(collect(keys(residual)))"))
    extra = setdiff(keys(residual), expected)
    isempty(extra) || @warn "the PDE function returns fields that do not correspond to any unknown; they are ignored (return a second NamedTuple to save objects on the grid)" ignored = collect(extra) expected = collect(expected)
    return nothing
end

# evaluate the PDE function on the whole grid and record its saved objects
function _fill_saved!(@nospecialize(saved), pde, stategrid::StateGrid, solutionnames, y_array::AbstractArray, bc_array::AbstractArray)
    for i in eachindex(stategrid)
        solution = differentiate(solutionnames, stategrid, y_array, i, bc_array)
        outi = pde(stategrid[i], solution)[2]
        for (k, v) in zip(values(saved), values(outi))
            k[i] = v
        end
    end
end

# create pde! that accepts and returns AbstractVector rather than AbstractArrays
function pde!(pde, stategrid::StateGrid, solutionnames, ydot::AbstractVector, y::AbstractVector, bc_array::AbstractArray, ysize::NTuple)
    y_array = reshape(y, ysize...)
    ydot_array = reshape(ydot, ysize...)
    vec(pde!(pde, stategrid, solutionnames, ydot_array, y_array, bc_array))
end

function pde!(pde, stategrid::StateGrid, solutionnames, ydot_array::AbstractArray, y_array::AbstractArray, bc_array::AbstractArray)
    for i in eachindex(stategrid)
        solution = differentiate(solutionnames, stategrid, y_array, i, bc_array)
        outi = pde(stategrid[i], solution)
        if isa(outi[1], Number)
            _write_residual!(ydot_array, solutionnames, outi, i)
        else
            _write_residual!(ydot_array, solutionnames, outi[1], i)
        end
    end
    return ydot_array
end

# write the returned time derivatives (`vt`, …) into slot i of each unknown's slice of ydot
@generated function _write_residual!(ydot_array::AbstractArray, ::Val{solnames}, outi::NamedTuple, i::CartesianIndex) where {solnames}
    N = length(solnames)
    quote
         $(Expr(:meta, :inline))
         $(Expr(:block, [Expr(:call, :setindex!, :ydot_array, Expr(:call, :getproperty, :outi, Meta.quot(Symbol(solnames[k], :t))), :i, k) for k in 1:N]...))
    end
end
