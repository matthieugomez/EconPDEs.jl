"""
    pdesolve(pde, grid, guess [, τs]; kwargs...)

Solve a system of nonlinear ODEs/PDEs — typically Hamilton–Jacobi–Bellman equations —
by finite differences. Handles an arbitrary number of coupled unknown functions on a
state grid of one, two, or three state variables.

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
  `τs[end]`, and `result.zero[i]` is the solution at time `τs[i]`.

### Keyword arguments
* `bc`: boundary derivatives, a `NamedTuple` (or `OrderedDict`) mapping
  `Symbol(unknown, state)` to a `(lower, upper)` tuple (scalars, or arrays matching the
  boundary slice). Entries may be given for any subset of (unknown, state) pairs; omitted
  pairs default to reflecting boundaries (zero first derivative at the boundary). Both
  entries are `∂v/∂state` oriented in the increasing-state direction.
* `is_algebraic`: a `NamedTuple` (or `OrderedDict`) of `Bool`s marking equations that are
  algebraic (no time derivative) rather than PDEs. Defaults to all `false`.
* `lower_bound`, `upper_bound`: lower/upper bounds on the unknowns for HJB variational
  inequalities (optimal stopping), solved as a mixed complementarity problem. Default to
  `-Inf`/`Inf` (unbounded). The Unicode keywords `y̲`/`ȳ` are deprecated aliases.
* `method`: `:newton` (default) or `:trust_region`.
* `maxdist`: convergence tolerance on the residual. Defaults to `sqrt(eps())`.
* `iterations`: maximum number of pseudo-transient iterations. Defaults to 100.
* `Δ`: initial pseudo-transient time step. Defaults to `1.0`; pass `Δ = Inf` to solve the
  stationary residual in one nonlinear solve, with no continuation.
* `autodiff`: `:forward` (default), `:finite`, or `:central`. With one to three state
  variables the Jacobian is always computed by colored sparse finite differences
  (`:forward` and `:finite` both use forward differences; `:central` uses central
  differences); forward-mode AD is used only when no sparsity pattern is available
  (four or more states).
* `verbose`: print a one-line problem summary, convergence progress, and a final
  convergence summary. Defaults to `true`. With `verbose = false` a successful solve
  prints nothing — convergence failures are always reported with `@warn` — so `pdesolve`
  can run inside a loop (e.g. an estimation) without flooding the log.
* Further solver knobs (`scale`, `minΔ`, `maxΔ`, `inner_iterations`, `innerdist`,
  `inner_verbose`, `reformulation`, `autoscale`) are forwarded to
  [`finiteschemesolve`](@ref); see its docstring, especially when a solve stalls.
* `check_monotonicity`: if true, warn when the assembled residual Jacobian has same-variable
  spatial off-diagonal entries with the wrong monotonicity sign. Defaults to false. A flipped
  sign convention — returning `RHS - ρv` instead of `-(RHS - ρv)`, which makes the Jacobian
  diagonal negative at every grid point — is detected and warned about regardless of this
  option (whenever `Δ` is finite).
* `monotonicity_tol`: tolerance for the monotonicity check. Defaults to 1e-6.
* `monotonicity_max_warnings`: maximum number of monotonicity warnings to print. Defaults to 5.

Returns an `EconPDEResult` with fields `zero` (the solved unknowns), `residual_norm`,
and `saved` (the saved objects, together with the solved unknowns). `result.converged`
reports whether the residual met the tolerance (across all times, for a time-dependent
problem). `result.optional` is kept as a backward-compatible alias for `result.saved`.
"""
function pdesolve(pde, @nospecialize(grid), @nospecialize(guess), τs::Union{Nothing, AbstractVector} = nothing; is_algebraic = nothing, bc = nothing, maxdist = sqrt(eps()), autodiff = :forward, verbose = true, check_monotonicity = false, monotonicity_tol = 1e-6, monotonicity_max_warnings = 5, kwargs...)
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
    # concatenate the per-unknown arrays into one array with a trailing unknown dimension
    guess_array = catlast(values(guess))
    is_algebraic_array = catlast(values(is_algebraic))
    bc_array = _bc_array(bc, guess, grid)
    # create sparsity pattern and its coloring (closed-form, from the stencil structure)
    J0 = sparse_jacobian(stategrid, guess)
    colorvec = J0 === nothing ? nothing : stencil_colors(S, length(guess))
    # the sign-convention check runs by default; the stencil (monotonicity) warnings are opt-in
    monotonicity_check = MonotonicityChecker(stategrid, guess; tol = monotonicity_tol, max_warnings = monotonicity_max_warnings, check_stencils = check_monotonicity)
    if τs === nothing
        return _solve_stationary(pde, stategrid, guess, guess_array, is_algebraic_array, bc_array, J0, colorvec, monotonicity_check; maxdist = maxdist, autodiff = autodiff, verbose = verbose, kwargs...)
    else
        issorted(τs) || throw(ArgumentError("The set of times must be increasing."))
        return _solve_backward(pde, stategrid, guess, τs, guess_array, is_algebraic_array, bc_array, J0, colorvec, monotonicity_check; maxdist = maxdist, autodiff = autodiff, verbose = verbose, kwargs...)
    end
end

# Stationary problem: one pseudo-transient continuation solve of the stationary residual.
function _solve_stationary(pde, stategrid::StateGrid, @nospecialize(guess), guess_array, is_algebraic_array, bc_array, J0, colorvec, monotonicity_check; maxdist, autodiff, verbose, kwargs...)
    ynames = tuple(keys(guess)...)
    solutionnames = Val(ynames)
    y = OrderedDict(first(p) => collect(last(p)) for p in pairs(guess))
    saved = _init_saved(pde, stategrid, solutionnames, ynames, guess_array, bc_array)
    verbose && println("Solving for ", _problem_description(guess, size(stategrid)))
    G! = _residual_barrier((ydot, yvec) -> pde!(pde, stategrid, solutionnames, ydot, yvec, bc_array, size(guess_array)), J0)
    y_array, residual_norm = finiteschemesolve(G!, vec(guess_array); is_algebraic = vec(is_algebraic_array), J0 = J0, colorvec = colorvec, maxdist = maxdist, autodiff = autodiff, verbose = verbose, monotonicity_check = monotonicity_check, kwargs...)
    y_array = reshape(y_array, size(guess_array)...)
    _copy_solution!(y, y_array)
    if saved !== nothing
        _fill_saved!(saved, pde, stategrid, solutionnames, y_array, bc_array)
        saved = merge(y, saved)
    end
    return EconPDEResult(y, residual_norm, saved, maxdist)
end

# Time-dependent problem: starting from the terminal condition at τs[end], take one
# implicit time step per interval of the time grid, recording the solution at each time.
function _solve_backward(pde, stategrid::StateGrid, @nospecialize(guess), τs, guess_array, is_algebraic_array, bc_array, J0, colorvec, monotonicity_check; maxdist, autodiff, verbose, kwargs...)
    ynames = tuple(keys(guess)...)
    solutionnames = Val(ynames)
    ys = [OrderedDict(first(p) => collect(last(p)) for p in pairs(guess)) for τ in τs]
    # the PDE function may or may not take a time argument; check once, not at every time step
    has_time = hasmethod(pde, Tuple{NamedTuple, NamedTuple, Number})
    pde_at = τ -> (has_time ? ((state, u) -> pde(state, u, τ)) : pde)
    saved = _init_saved(pde_at(τs[end]), stategrid, solutionnames, ynames, guess_array, bc_array)
    saveds = (saved === nothing) ? nothing : [deepcopy(saved) for τ in τs]
    y_array = guess_array
    residual_norms = zeros(length(τs))
    # the sparsity pattern is fixed, so build the coloring and Jacobian cache once
    # rather than inside every backward time step
    J0c, fdcache = _sparse_fd_setup(J0, vec(guess_array), autodiff, colorvec)
    tstart = time()
    if verbose
        @printf "Solving for %s, backward from τ = %g to %g (%d steps)\n" _problem_description(guess, size(stategrid)) τs[end] τs[1] (length(τs) - 1)
        @printf "    Time Residual\n"
        @printf "-------- --------\n"
    end
    for iτ in length(τs):(-1):1
        _copy_solution!(ys[iτ], y_array)
        pde_onestep = pde_at(τs[iτ])
        if saved !== nothing
            _fill_saved!(saveds[iτ], pde_onestep, stategrid, solutionnames, y_array, bc_array)
            saveds[iτ] = merge(ys[iτ], saveds[iτ])
        end
        if iτ > 1
            G! = _residual_barrier((ydot, yvec) -> pde!(pde_onestep, stategrid, solutionnames, ydot, yvec, bc_array, size(guess_array)), J0c)
            y_array, residual_norms[iτ] = implicit_timestep(G!, vec(y_array), τs[iτ] - τs[iτ-1]; is_algebraic = vec(is_algebraic_array), verbose = false, J0 = J0c, fdcache = fdcache, maxdist = maxdist, autodiff = autodiff, monotonicity_check = monotonicity_check, kwargs...)
            if verbose
                @printf "%8g %8.2e\n" τs[iτ-1] residual_norms[iτ]
            end
            # an unconverged step would otherwise be silently accepted and propagated
            # to all earlier times; `!(x <= tol)` also catches a NaN residual
            if !(residual_norms[iτ] <= maxdist)
                @warn "the implicit time step at τ = $(τs[iτ-1]) did not converge (residual norm: $(@sprintf("%.2e", residual_norms[iτ]))); the solution at this and earlier times may be inaccurate — try a finer time grid or more `iterations`"
            end
            y_array = reshape(y_array, size(guess_array)...)
        end
    end
    if verbose
        nfailed = count(iτ -> !(residual_norms[iτ] <= maxdist), 2:length(τs))
        if nfailed == 0
            @printf "Completed %d time steps (%s)\n" (length(τs) - 1) _elapsed(time() - tstart)
        else
            @printf "Completed %d time steps (%s): %d steps did not converge\n" (length(τs) - 1) _elapsed(time() - tstart) nfailed
        end
    end
    return EconPDEResult(ys, residual_norms, saveds, maxdist)
end

# e.g. "2 unknowns (pA, pB) on a 30×40 grid" — shared by the preambles of both solve paths
function _problem_description(@nospecialize(guess), S)
    n = length(guess)
    string(n, n == 1 ? " unknown (" : " unknowns (", join(keys(guess), ", "), ") on a ",
           length(S) == 1 ? "$(S[1])-point" : join(S, "×"), " grid")
end

catlast(iter) = cat(iter...; dims = ndims(first(iter)) + 1)

# Accept either an OrderedDict (or any symbol-keyed collection) or a NamedTuple.
# A plain Dict is rejected: its iteration order is arbitrary, and the order of names
# determines how the solution arrays are laid out, so accepting one risks a silent
# transposition of the guess (or of the result) on square grids.
_asnamedtuple(x::NamedTuple) = x
_asnamedtuple(x::Dict) = throw(ArgumentError("pass a NamedTuple, e.g. `(; k = ...)`, or an OrderedDict — a Dict has no fixed iteration order, and the order of the names determines how arrays are laid out"))
_asnamedtuple(x) = NamedTuple(x)

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

# Expand the `is_algebraic` flags into a NamedTuple of arrays (one Bool per grid point),
# with the same names and order as `guess`. Defaults to all-false when not provided.
_fill_is_algebraic(::Nothing, guess, S) = map(_ -> fill(false, S), guess)
function _fill_is_algebraic(is_algebraic, guess, S)
    ia = _asnamedtuple(is_algebraic)
    issetequal(keys(ia), keys(guess)) || throw(ArgumentError("`is_algebraic` must have the same names as the initial guess: got $(collect(keys(ia))), expected $(collect(keys(guess)))"))
    # reorder to match the guess, so the flags line up with the right unknowns
    NamedTuple{keys(guess)}(map(n -> fill(ia[n], S), keys(guess)))
end

_bc_array(::Nothing, guess, grid) = zeros(size(first(values(guess)))..., length(guess))
function _bc_array(bc, guess, grid)
    bc = _asnamedtuple(bc)
    keys_grid = collect(keys(grid))
    valid = [Symbol(yname, s) for yname in keys(guess) for s in keys_grid]
    for key in keys(bc)
        key in valid || throw(ArgumentError("unknown `bc` entry `$key`: valid entries are $(valid), i.e. Symbol(unknown, state)"))
    end
    out = _bc_array(nothing, guess, grid)
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

# Call the PDE function once at the first grid point to learn whether it returns a second
# NamedTuple of objects to save on the grid; if so, allocate one array per saved object.
function _init_saved(pde, stategrid::StateGrid, solutionnames, ynames, y_array::AbstractArray, bc_array::AbstractArray)
    i0 = first(eachindex(stategrid))
    derivatives = differentiate(solutionnames, stategrid, y_array, i0, bc_array)
    result = pde(stategrid[i0], derivatives)
    residual, optional = _split_pde_output(result)
    _check_residual_names(residual, ynames)
    optional === nothing && return nothing
    return OrderedDict(key => Array{Float64}(undef, size(stategrid)) for key in keys(optional))
end

# Discriminate the two allowed return shapes by type: a single-return model gives a
# NamedTuple of time derivatives, a two-return model gives (residual, optional).
# (pde! performs the same discrimination per grid point via `isa(outi[1], Number)`.)
_split_pde_output(result::NamedTuple) = (result, nothing)
_split_pde_output(result::Tuple{NamedTuple, NamedTuple}) = result
_split_pde_output(result) = throw(ArgumentError("the PDE function must return a NamedTuple of time derivatives, e.g. `(; vt)`, optionally followed by a second NamedTuple of objects to save on the grid — got a value of type $(typeof(result))"))

function _check_residual_names(residual::NamedTuple, ynames)
    expected = map(n -> Symbol(n, :t), ynames)
    missing_names = setdiff(expected, keys(residual))
    isempty(missing_names) || throw(ArgumentError("the PDE function must return one time derivative per unknown, named `Symbol(unknown, :t)`: expected $(collect(expected)), got $(collect(keys(residual)))"))
    extra = setdiff(keys(residual), expected)
    isempty(extra) || @warn "the PDE function returns fields that do not correspond to any unknown; they are ignored (return a second NamedTuple to save objects on the grid)" ignored = collect(extra) expected = collect(expected)
    return nothing
end

# copy the concatenated solution array back into the per-unknown arrays of the OrderedDict
function _copy_solution!(@nospecialize(y), y_array::AbstractArray)
    for (i, v) in enumerate(values(y))
        v[:] = selectdim(y_array, ndims(y_array), i)
    end
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
