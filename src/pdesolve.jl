"""
    pdesolve(f, grid, guess [, τs]; kwargs...)

Solve a system of nonlinear ODEs/PDEs — typically Hamilton–Jacobi–Bellman equations —
by finite differences. Handles an arbitrary number of coupled unknown functions on a
state grid of one, two, or three state variables.

### Positional arguments
* `f`: the local equation, a function `(state, u) -> out` (or `(state, u, t) -> out` for a
  time-dependent equation), where
    - `state` is a `NamedTuple` with the current grid point (one entry per state variable),
    - `u` is a `NamedTuple` with each unknown function and its finite-difference derivatives
      at that point (e.g. `v`, `vk_up`, `vk_down`, `vkk`, …),
    - `out` is a `NamedTuple` with one time derivative per unknown (e.g. `(; vt)`).
  `f` may also return a second `NamedTuple` of objects to save on the grid; these are
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
function pdesolve(apm, @nospecialize(grid), @nospecialize(guess), τs::Union{Nothing, AbstractVector} = nothing; is_algebraic = nothing, bc = nothing, verbose = true, check_monotonicity = false, monotonicity_tol = 1e-6, monotonicity_max_warnings = 5, kwargs...)
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
    if τs isa AbstractVector
        issorted(τs) || throw(ArgumentError("The set of times must be increasing."))
        ys = [OrderedDict(first(p) => collect(last(p)) for p in pairs(guess)) for t in τs]
    else
        y = OrderedDict(first(p) => collect(last(p)) for p in pairs(guess))
    end
    # convert to Matrix
    guess_M = catlast(values(guess))
    is_algebraic_M = catlast(values(is_algebraic))
    bc_M = _Array_bc(bc, guess, grid)

    # create sparsity
    J0 = sparse_jacobian(stategrid, guess)
    # the sign-convention check runs by default; the stencil (monotonicity) warnings are opt-in
    monotonicity_check = MonotonicityChecker(stategrid, guess; tol = monotonicity_tol, max_warnings = monotonicity_max_warnings, check_stencils = check_monotonicity)
    ynames = tuple(keys(guess)...)
    solutionnames = Val(ynames)

    # iterate on time
    if τs isa AbstractVector
        apm_onestep = hasmethod(apm, Tuple{NamedTuple, NamedTuple, Number}) ? (state, grid) -> apm(state, grid, τs[end]) : apm
        a = get_a(apm_onestep, stategrid, solutionnames, ynames, guess_M, bc_M)
        as = (a === nothing) ? nothing : [deepcopy(a) for τ in τs]
        y_M = guess_M
        residual_norms = zeros(length(τs))
        maxdist = get(kwargs, :maxdist, sqrt(eps()))
        # the sparsity pattern is fixed, so build the coloring and Jacobian cache once
        # rather than inside every backward time step
        J0c, fdcache = _sparse_fd_setup(J0, vec(guess_M), get(kwargs, :autodiff, :forward))
        tstart = time()
        if verbose
            @printf "Solving for %s, backward from τ = %g to %g (%d steps)\n" _problem_description(guess, S) τs[end] τs[1] (length(τs) - 1)
            @printf "    Time Residual\n"
            @printf "-------- --------\n"
        end
        for iτ in length(τs):(-1):1
            _setindex!(ys[iτ], y_M)
            apm_onestep = hasmethod(apm, Tuple{NamedTuple, NamedTuple, Number}) ? (state, grid) -> apm(state, grid, τs[iτ]) : apm
            if a !== nothing
                _setindex!(as[iτ], apm_onestep, stategrid, solutionnames, y_M, bc_M)
                as[iτ] = merge(ys[iτ], as[iτ])
            end
            if iτ > 1
                y_M, residual_norms[iτ] = implicit_timestep((ydot, y) -> hjb!(apm_onestep, stategrid, solutionnames, ydot, y, bc_M, size(guess_M)), vec(y_M), τs[iτ] - τs[iτ-1]; is_algebraic = vec(is_algebraic_M), verbose = false, J0 = J0c, fdcache = fdcache, monotonicity_check = monotonicity_check, kwargs...)
                if verbose
                    @printf "%8g %8.2e\n" τs[iτ-1] residual_norms[iτ]
                end
                # an unconverged step would otherwise be silently accepted and propagated
                # to all earlier times; `!(x <= tol)` also catches a NaN residual
                if !(residual_norms[iτ] <= maxdist)
                    @warn "the implicit time step at τ = $(τs[iτ-1]) did not converge (residual norm: $(@sprintf("%.2e", residual_norms[iτ]))); the solution at this and earlier times may be inaccurate — try a finer time grid or more `iterations`"
                end
                y_M = reshape(y_M, size(guess_M)...)
            end
        end
        if verbose
            nfailed = count(iτ -> !(residual_norms[iτ] <= maxdist), 2:length(τs))
            if nfailed == 0
                @printf "Completed %d time steps (%s): max residual %.2e ≤ tolerance %.2e\n" (length(τs) - 1) _elapsed(time() - tstart) maximum(residual_norms) maxdist
            else
                @printf "Completed %d time steps (%s): %d steps did not converge\n" (length(τs) - 1) _elapsed(time() - tstart) nfailed
            end
        end
        return EconPDEResult(ys, residual_norms, as, maxdist)
    else
        a = get_a(apm, stategrid, solutionnames, ynames, guess_M, bc_M)
        maxdist = get(kwargs, :maxdist, sqrt(eps()))
        verbose && println("Solving for ", _problem_description(guess, S))
        y_M, residual_norm = finiteschemesolve((ydot, y) -> hjb!(apm, stategrid, solutionnames, ydot, y, bc_M, size(guess_M)), vec(guess_M); is_algebraic = vec(is_algebraic_M),  J0 = J0, verbose = verbose, monotonicity_check = monotonicity_check, kwargs... )
        y_M = reshape(y_M, size(guess_M)...)
        _setindex!(y, y_M)
        if a !== nothing
            _setindex!(a, apm, stategrid, solutionnames, y_M, bc_M)
            a = merge(y, a)
        end
        return EconPDEResult(y, residual_norm, a, maxdist)
    end
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

_Array_bc(::Nothing, guess, grid) = zeros(size(first(values(guess)))..., length(guess))
function _Array_bc(bc, guess, grid)
    bc = _asnamedtuple(bc)
    keys_grid = collect(keys(grid))
    valid = [Symbol(yname, s) for yname in keys(guess) for s in keys_grid]
    for key in keys(bc)
        key in valid || throw(ArgumentError("unknown `bc` entry `$key`: valid entries are $(valid), i.e. Symbol(unknown, state)"))
    end
    bc_M = _Array_bc(nothing, guess, grid)
    for (k, yname) in enumerate(keys(guess))
        bck = selectdim(bc_M, ndims(bc_M), k)
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
    return bc_M
end


function get_a(apm, stategrid::StateGrid, solutionnames, ynames, y_M::AbstractArray, bc_M::AbstractArray)
    i0 = first(eachindex(stategrid))
    derivatives = differentiate(solutionnames, stategrid, y_M, i0, bc_M)
    result = apm(stategrid[i0], derivatives)
    residual, optional = _split_pde_output(result)
    _check_residual_names(residual, ynames)
    optional === nothing && return nothing
    return OrderedDict(a_key => Array{Float64}(undef, size(stategrid)) for a_key in keys(optional))
end

# Discriminate the two allowed return shapes by type: a single-return model gives a
# NamedTuple of time derivatives, a two-return model gives (residual, optional).
# (hjb! performs the same discrimination per grid point via `isa(outi[1], Number)`.)
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


function sparse_jacobian(stategrid::StateGrid, @nospecialize(guess))
    s = size(stategrid)
    if 1 <= ndims(stategrid) <= 3
        return local_stencil_jacobian(s, length(guess))
    else
        return nothing
    end
end

function local_stencil_jacobian(s::NTuple{N, Int}, F::Int) where {N}
    l = prod(s)
    nnz_max = l * F^2 * 3^N
    rows = Vector{Int}(undef, nnz_max)
    cols = Vector{Int}(undef, nnz_max)
    k = 0
    state_indices = CartesianIndices(s)
    for fout in 1:F, ci in state_indices
        row = LinearIndices(s)[ci] + (fout - 1) * l
        stencil_ranges = ntuple(d -> max(ci[d] - 1, 1):min(ci[d] + 1, s[d]), N)
        for fin in 1:F, neighbor in CartesianIndices(stencil_ranges)
            k += 1
            rows[k] = row
            cols[k] = LinearIndices(s)[neighbor] + (fin - 1) * l
        end
    end
    return sparse(rows[1:k], cols[1:k], ones(k), l * F, l * F)
end


function _setindex!(@nospecialize(y), y_M::AbstractArray)
    for (i, v) in enumerate(values(y))
        v[:] = selectdim(y_M, ndims(y_M), i)
    end
end

function _setindex!(@nospecialize(a), apm, stategrid::StateGrid, solutionnames, y_M::AbstractArray, bc_M::AbstractArray)
    for i in eachindex(stategrid)
         solution = differentiate(solutionnames, stategrid, y_M, i, bc_M)
         outi = apm(stategrid[i], solution)[2]
         for (k, v) in zip(values(a), values(outi))
            k[i] = v
         end
     end
 end


 # create hjb! that accepts and returns AbstractVector rather than AbstractArrays
 function hjb!(apm, stategrid::StateGrid, solutionnames, ydot::AbstractVector, y::AbstractVector, bc_M::AbstractArray, ysize::NTuple)
     y_M = reshape(y, ysize...)
     ydot_M = reshape(ydot, ysize...)
     vec(hjb!(apm, stategrid, solutionnames, ydot_M, y_M, bc_M))
 end

 function hjb!(apm, stategrid::StateGrid, solutionnames, ydot_M::AbstractArray, y_M::AbstractArray, bc_M::AbstractArray)
     for i in eachindex(stategrid)
         solution = differentiate(solutionnames, stategrid, y_M, i, bc_M)
         outi = apm(stategrid[i], solution)
         if isa(outi[1], Number)
            _setindex!(ydot_M, solutionnames, outi, i)
        else
            _setindex!(ydot_M, solutionnames, outi[1], i)
        end
     end
     return ydot_M
 end

 @generated function _setindex!(ydot_M::AbstractArray, ::Val{solnames}, outi::NamedTuple, i::CartesianIndex) where {solnames}
     N = length(solnames)
     quote
          $(Expr(:meta, :inline))
          $(Expr(:block, [Expr(:call, :setindex!, :ydot_M, Expr(:call, :getproperty, :outi, Meta.quot(Symbol(solnames[k], :t))), :i, k) for k in 1:N]...))
     end
 end
