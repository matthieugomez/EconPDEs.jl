mutable struct MonotonicityChecker
    stategrid
    ynames::Vector{Symbol}
    ysize::Tuple
    tol::Float64
    max_warnings::Int
    warned::Set{Tuple{Int, Int}}
    suppressed::Bool
end

function MonotonicityChecker(stategrid::StateGrid, @nospecialize(yend); tol = 1e-6, max_warnings = 5)
    return MonotonicityChecker(
        stategrid,
        collect(keys(yend)),
        tuple(size(stategrid)..., length(yend)),
        tol,
        max_warnings,
        Set{Tuple{Int, Int}}(),
        false,
    )
end

_run_monotonicity_check!(::Nothing, J, y) = nothing
_run_monotonicity_check!(checker::Function, J, y) = checker(J, y)
_run_monotonicity_check!(checker::MonotonicityChecker, J, y) = _check_monotonicity!(checker, J)

_sparse_for_check(J::SparseMatrixCSC) = J
_sparse_for_check(J) = sparse(J)

function _check_monotonicity!(checker::MonotonicityChecker, J)
    checker.max_warnings <= 0 && return nothing
    (checker.suppressed && length(checker.warned) >= checker.max_warnings) && return nothing

    cart = CartesianIndices(checker.ysize)
    size(J, 1) == length(cart) || return nothing
    size(J, 2) == length(cart) || return nothing

    rows, cols, vals = findnz(_sparse_for_check(J))
    for k in eachindex(vals)
        val = vals[k]
        val > checker.tol || continue

        row = rows[k]
        col = cols[k]
        row == col && continue
        (row, col) in checker.warned && continue

        rowI = cart[row]
        colI = cart[col]
        _is_spatial_neighbor(checker, rowI, colI) || continue

        if length(checker.warned) >= checker.max_warnings
            _warn_suppressed!(checker)
            return nothing
        end

        push!(checker.warned, (row, col))
        _warn_nonmonotone(checker, rowI, colI, val)
    end
    return nothing
end

function _is_spatial_neighbor(checker::MonotonicityChecker, rowI::CartesianIndex, colI::CartesianIndex)
    nstate = ndims(checker.stategrid)
    rowI[nstate + 1] == colI[nstate + 1] || return false

    changed = false
    for d in 1:nstate
        delta = colI[d] - rowI[d]
        abs(delta) <= 1 || return false
        changed |= delta != 0
    end
    return changed
end

function _warn_nonmonotone(checker::MonotonicityChecker, rowI::CartesianIndex, colI::CartesianIndex, residual_weight)
    varname = checker.ynames[rowI[ndims(checker.stategrid) + 1]]
    @warn "Non-monotone finite-difference stencil detected: positive residual off-diagonal implies a negative generator weight" variable=varname state=_format_state(checker, rowI) neighbor=_format_neighbor(checker, rowI, colI) residual_jacobian=residual_weight generator_weight=-residual_weight stencil_hint=_format_stencil_hint(checker, rowI, colI, varname)
    return nothing
end

function _warn_suppressed!(checker::MonotonicityChecker)
    checker.suppressed && return nothing
    checker.suppressed = true
    @warn "Further monotonicity warnings suppressed" limit=checker.max_warnings
    return nothing
end

function _format_state(checker::MonotonicityChecker, I::CartesianIndex)
    names = keys(checker.stategrid.x)
    parts = Vector{String}(undef, ndims(checker.stategrid))
    for d in 1:ndims(checker.stategrid)
        parts[d] = string(names[d], "=", _format_number(checker.stategrid.x[d][I[d]]))
    end
    return join(parts, ", ")
end

function _format_neighbor(checker::MonotonicityChecker, rowI::CartesianIndex, colI::CartesianIndex)
    names = keys(checker.stategrid.x)
    parts = String[]
    for d in 1:ndims(checker.stategrid)
        delta = colI[d] - rowI[d]
        if delta > 0
            push!(parts, string(names[d], "+"))
        elseif delta < 0
            push!(parts, string(names[d], "-"))
        end
    end
    return join(parts, ", ")
end

function _format_stencil_hint(checker::MonotonicityChecker, rowI::CartesianIndex, colI::CartesianIndex, varname::Symbol)
    names = keys(checker.stategrid.x)
    changed = Int[]
    deltas = Int[]
    for d in 1:ndims(checker.stategrid)
        delta = colI[d] - rowI[d]
        if delta != 0
            push!(changed, d)
            push!(deltas, delta)
        end
    end

    if length(changed) == 1
        d = changed[1]
        suffix = deltas[1] > 0 ? "_up" : "_down"
        return Symbol(varname, names[d], suffix)
    elseif length(changed) == 2
        d1, d2 = changed
        suffix = sign(deltas[1]) == sign(deltas[2]) ? "_up" : "_down"
        return Symbol(varname, names[d1], names[d2], suffix)
    else
        return Symbol(varname, "_neighbor")
    end
end

_format_number(x::Real) = @sprintf("%.6g", x)
_format_number(x) = sprint(show, x)
