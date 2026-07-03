mutable struct MonotonicityChecker
    stategrid
    ynames::Vector{Symbol}
    ysize::Tuple
    tol::Float64
    max_warnings::Int
    warned::Set{Tuple{Int, Int}}
    suppressed::Bool
    check_stencils::Bool
    sign_checked::Bool
    disabled::Bool
end

function MonotonicityChecker(stategrid::StateGrid, @nospecialize(guess); tol = 1e-6, max_warnings = 5, check_stencils = true)
    return MonotonicityChecker(
        stategrid,
        collect(keys(guess)),
        tuple(size(stategrid)..., length(guess)),
        tol,
        max_warnings,
        Set{Tuple{Int, Int}}(),
        false,
        check_stencils,
        false,
        false,
    )
end

_run_monotonicity_check!(::Nothing, J, y, Δ, is_algebraic) = nothing
function _run_monotonicity_check!(checker::MonotonicityChecker, J, y, Δ, is_algebraic)
    _check_sign_convention!(checker, J, Δ, is_algebraic)
    checker.check_stencils && _check_monotonicity!(checker, J)
    return nothing
end

_try_run_monotonicity_check!(::Nothing, J, y, Δ, is_algebraic) = nothing
function _try_run_monotonicity_check!(checker::MonotonicityChecker, J, y, Δ, is_algebraic)
    checker.disabled && return nothing
    try
        _run_monotonicity_check!(checker, J, y, Δ, is_algebraic)
    catch err
        checker.disabled = true
        @debug "Monotonicity diagnostic failed; continuing without diagnostics" exception = (err, catch_backtrace())
    end
    return nothing
end

# Under the `vt = -(RHS - ρv)` convention, the diagonal of the residual Jacobian is positive
# (ρ plus the outflow rates of the discretized generator). If the user returns `RHS - ρv`
# instead, every diagonal entry flips negative and the false transient is unstable — the most
# common reason a solve diverges from iteration 1. Checked once, on the first assembled
# Jacobian. Skipped when Δ = Inf: a single Newton solve is sign-agnostic, so either sign is
# then legitimate.
function _check_sign_convention!(checker::MonotonicityChecker, J, Δ, is_algebraic)
    checker.sign_checked && return nothing
    isfinite(Δ) || return nothing
    checker.sign_checked = true
    n = prod(checker.ysize)
    (size(J, 1) == n && size(J, 2) == n) || return nothing
    npos = nneg = 0
    for i in 1:n
        is_algebraic[i] && continue
        # the assembled Jacobian is that of the implicit time step, which adds 1/Δ to
        # non-algebraic diagonal entries; remove it to recover the residual's own diagonal
        d = J[i, i] - 1 / Δ
        abs(d) > checker.tol || continue
        d < 0 ? (nneg += 1) : (npos += 1)
    end
    if nneg > 0 && npos == 0
        @warn "The residual Jacobian has a negative diagonal at every grid point: the PDE function most likely returns `RHS - ρv` instead of `-(RHS - ρv)`, which makes the false transient unstable (the solve typically diverges from iteration 1). Flip the sign of the returned time derivative."
    end
    return nothing
end

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
