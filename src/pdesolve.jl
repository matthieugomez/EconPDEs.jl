"""
    pdesolve(f, grid, yend [, τs]; kwargs...)

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
  `f` may also return a second `NamedTuple` of extra objects to save on the grid; these are
  returned in `result.optional`.
* `grid`: a `NamedTuple` (or `OrderedDict`) mapping each state variable name to an
  `AbstractVector` (its grid), e.g. `(; k = range(...))`.
* `yend`: a `NamedTuple` (or `OrderedDict`) mapping each unknown name to an array of initial
  values with the same shape as the grid. For a time-dependent problem, it is the terminal
  value at `τs[end]`.
* `τs` (optional): an increasing time grid. The equation is then solved backward from
  `τs[end]`, and `result.zero[i]` is the solution at time `τs[i]`.

### Keyword arguments
* `bc`: boundary derivatives, a `NamedTuple` (or `OrderedDict`) mapping each
  `Symbol(unknown, state)` to a `(lower, upper)` tuple. Defaults to reflecting boundaries
  (zero outward first derivative).
* `is_algebraic`: a `NamedTuple` (or `OrderedDict`) of `Bool`s marking equations that are
  algebraic (no time derivative) rather than PDEs. Defaults to all `false`.
* `y̲`, `ȳ`: lower/upper bounds for HJB variational inequalities (optimal stopping), solved as a
  mixed complementarity problem. Default to `-Inf`/`Inf` (unbounded).
* `method`: `:newton` (default) or `:trust_region`.
* `maxdist`: convergence tolerance on the residual. Defaults to `sqrt(eps())`.
* `iterations`: maximum number of pseudo-transient iterations. Defaults to 100.
* `Δ`: initial pseudo-transient time step. Defaults to `1.0`; pass `Δ = Inf` to solve in a
  single Newton step (no continuation).
* `autodiff`: `:forward` (default), `:finite`, or `:central`.
* `verbose`: print convergence progress. Defaults to `true`.
* `check_monotonicity`: if true, warn when the assembled residual Jacobian has same-variable
  spatial off-diagonal entries with the wrong monotonicity sign. Defaults to false.
* `monotonicity_tol`: tolerance for the monotonicity check. Defaults to 1e-6.
* `monotonicity_max_warnings`: maximum number of monotonicity warnings to print. Defaults to 5.

Returns an `EconPDEResult` with fields `zero` (the solved unknowns), `residual_norm`,
and `optional` (the saved objects, together with the solved unknowns).
"""
function pdesolve(apm, @nospecialize(grid), @nospecialize(yend), τs::Union{Nothing, AbstractVector} = nothing; is_algebraic = nothing, bc = nothing, verbose = true, check_monotonicity = false, monotonicity_tol = 1e-6, monotonicity_max_warnings = 5, kwargs...)
    # `grid`, `yend`, `is_algebraic`, and `bc` may be passed either as an OrderedDict or as a
    # NamedTuple. Normalize to a NamedTuple so everything below runs on one uniform representation.
    grid = _asnamedtuple(grid)
    yend = _asnamedtuple(yend)
    stategrid = StateGrid(grid)
    S = size(stategrid)
    all(size(v) == S for v in values(yend)) || throw(ArgumentError("The length of initial guess (e.g. terminal value) does not equal the length of the state space"))
    is_algebraic = _fill_is_algebraic(is_algebraic, yend, S)
    if τs isa AbstractVector
        issorted(τs) || throw(ArgumentError("The set of times must be increasing."))
        ys = [OrderedDict(first(p) => collect(last(p)) for p in pairs(yend)) for t in τs]
    else
        y = OrderedDict(first(p) => collect(last(p)) for p in pairs(yend))
    end
    # convert to Matrix
    yend_M = catlast(values(yend))
    is_algebraic_M = catlast(values(is_algebraic))
    bc_M = _Array_bc(bc, yend, grid)

    # create sparsity
    J0 = sparse_jacobian(stategrid, yend)
    monotonicity_check = check_monotonicity ? MonotonicityChecker(stategrid, yend; tol = monotonicity_tol, max_warnings = monotonicity_max_warnings) : nothing
    Tsolution = Type{tuple(keys(yend)...)}

    # iterate on time
    if τs isa AbstractVector
        apm_onestep = hasmethod(apm, Tuple{NamedTuple, NamedTuple, Number}) ? (state, grid) -> apm(state, grid, τs[end]) : apm
        a = get_a(apm_onestep, stategrid, Tsolution, yend_M, bc_M)
        as = (a === nothing) ? nothing : [deepcopy(a) for τ in τs]
        y_M = yend_M
        residual_norms = zeros(length(τs))
        if verbose
            @printf "    Time Residual\n"
            @printf "-------- --------\n"
        end
        for iτ in length(τs):(-1):1
            _setindex!(ys[iτ], y_M)
            apm_onestep = hasmethod(apm, Tuple{NamedTuple, NamedTuple, Number}) ? (state, grid) -> apm(state, grid, τs[iτ]) : apm
            if a !== nothing
                _setindex!(as[iτ], apm_onestep, stategrid, Tsolution, y_M, bc_M)
                as[iτ] = merge(ys[iτ], as[iτ])
            end
            if iτ > 1
                y_M, residual_norms[iτ] = implicit_timestep((ydot, y) -> hjb!(apm_onestep, stategrid, Tsolution, ydot, y, bc_M, size(yend_M)), vec(y_M), τs[iτ] - τs[iτ-1]; is_algebraic = vec(is_algebraic_M), verbose = false, J0 = J0, monotonicity_check = monotonicity_check, kwargs...)
                if verbose
                    @printf "%8g   %8.4e\n" τs[iτ-1] residual_norms[iτ]
                end
                y_M = reshape(y_M, size(yend_M)...)
            end
        end
        return EconPDEResult(ys, residual_norms, as)
    else
        a = get_a(apm, stategrid, Tsolution, yend_M, bc_M)
        y_M, residual_norm = finiteschemesolve((ydot, y) -> hjb!(apm, stategrid, Tsolution, ydot, y, bc_M, size(yend_M)), vec(yend_M); is_algebraic = vec(is_algebraic_M),  J0 = J0, verbose = verbose, monotonicity_check = monotonicity_check, kwargs... )
        y_M = reshape(y_M, size(yend_M)...)
        _setindex!(y, y_M)
        if a !== nothing
            _setindex!(a, apm, stategrid, Tsolution, y_M, bc_M)
            a = merge(y, a)
        end
        return EconPDEResult(y, residual_norm, a)
    end
end




catlast(iter) = cat(iter...; dims = ndims(first(iter)) + 1)

# Accept either an OrderedDict (or any symbol-keyed collection) or a NamedTuple.
_asnamedtuple(x::NamedTuple) = x
_asnamedtuple(x) = NamedTuple(x)

# Expand the `is_algebraic` flags into a NamedTuple of arrays (one Bool per grid point),
# with the same names and order as `yend`. Defaults to all-false when not provided.
_fill_is_algebraic(::Nothing, yend, S) = map(_ -> fill(false, S), yend)
function _fill_is_algebraic(is_algebraic, yend, S)
    ia = _asnamedtuple(is_algebraic)
    keys(ia) == keys(yend) || throw(ArgumentError("the terminal guess yend and the is_algebraic keyword argument must have the same names"))
    map(v -> fill(v, S), ia)
end

_Array_bc(::Nothing, yend, grid) = zeros(size(first(values(yend)))..., length(yend))
function _Array_bc(bc, yend, grid)
    keys_grid = collect(keys(grid))
    bc_M = _Array_bc(nothing, yend, grid)
    for (k, yname) in enumerate(keys(yend))
        if length(keys_grid) == 1
            bc_M[1, k], bc_M[end, k] = bc[Symbol(yname, keys_grid[1])]
        elseif length(keys_grid) == 2
            bc_M[1, :,  k],  bc_M[end, :, k] = bc[Symbol(yname, keys_grid[1])]
            bc_M[:, 1,  k], bc_M[:, end,  k] = bc[Symbol(yname, keys_grid[2])]
        elseif length(keys_grid) == 3
            bc_M[1, :, :, k], bc_M[end, :, :, k] = bc[Symbol(yname, keys_grid[1])]
            bc_M[:, 1, :, k], bc_M[:, end, :, k] = bc[Symbol(yname, keys_grid[2])]
            bc_M[:, :, 1, k], bc_M[:, :, end, k] = bc[Symbol(yname, keys_grid[3])]
        end
    end
    return bc_M
end


function get_a(apm, stategrid::StateGrid, Tsolution, y_M::AbstractArray, bc_M::AbstractArray)
    i0 = first(eachindex(stategrid))
    derivatives = differentiate(Tsolution, stategrid, y_M, i0, bc_M)
    result = apm(stategrid[i0], derivatives)
    # Discriminate by type, not by length: a single-return model gives a NamedTuple
    # of time derivatives whose first entry is a Number, while a two-return model gives
    # (residual, optional) whose first entry is the residual NamedTuple. Using length
    # here would misclassify any model with several unknowns (see hjb! for the same test).
    if isa(result[1], Number)
        return nothing
    else
        return OrderedDict(a_key => Array{Float64}(undef, size(stategrid)) for a_key in keys(result[2]))
    end
end


function sparse_jacobian(stategrid::StateGrid, @nospecialize(yend))
    s = size(stategrid)
    if 1 <= ndims(stategrid) <= 3
        return local_stencil_jacobian(s, length(yend))
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

function _setindex!(@nospecialize(a), apm, stategrid::StateGrid, Tsolution, y_M::AbstractArray, bc_M::AbstractArray)
    for i in eachindex(stategrid)
         solution = differentiate(Tsolution, stategrid, y_M, i, bc_M)
         outi = apm(stategrid[i], solution)[2]
         for (k, v) in zip(values(a), values(outi))
            k[i] = v
         end
     end
 end


 # create hjb! that accepts and returns AbstractVector rather than AbstractArrays
 function hjb!(apm, stategrid::StateGrid, Tsolution, ydot::AbstractVector, y::AbstractVector, bc_M::AbstractArray, ysize::NTuple)
     y_M = reshape(y, ysize...)
     ydot_M = reshape(ydot, ysize...)
     vec(hjb!(apm, stategrid, Tsolution, ydot_M, y_M, bc_M))
 end

 function hjb!(apm, stategrid::StateGrid, Tsolution, ydot_M::AbstractArray, y_M::AbstractArray, bc_M::AbstractArray)
     for i in eachindex(stategrid)
         solution = differentiate(Tsolution, stategrid, y_M, i, bc_M)
         outi = apm(stategrid[i], solution)
         if isa(outi[1], Number)
            _setindex!(ydot_M, Tsolution, outi, i)
        else
            _setindex!(ydot_M, Tsolution, outi[1], i)
        end
     end
     return ydot_M
 end

 @generated function _setindex!(ydot_M::AbstractArray, ::Type{Tsolution}, outi::NamedTuple, i::CartesianIndex) where {Tsolution}
     N = length(Tsolution.parameters[1])
     quote
          $(Expr(:meta, :inline))
          $(Expr(:block, [Expr(:call, :setindex!, :ydot_M, Expr(:call, :getproperty, :outi, Meta.quot(Symbol(Tsolution.parameters[1][k], :t))), :i, k) for k in 1:N]...))
     end
 end
