#========================================================================================

Finite-difference derivatives of the unknowns at one grid point.

`differentiate` returns a NamedTuple with, for each unknown, its value and all its
finite-difference derivatives (`v`, `vk_up`, `vk_down`, `vkk`, cross derivatives, …)
at grid point `icar`. It is @generated so the field names — which depend on the
unknown and state names — are computed at compile time and the NamedTuple is built
without any allocation.

========================================================================================#

# Helpers for building index expressions at compile time.
# Called inside the @generated function body to construct Expr nodes.

# Build y[i1, ..., iN, k] with optional substitutions at specific dimensions.
function _y_ref(N::Int, k::Int, subs::Pair{Int}...)
    idx = Any[Symbol("i", d) for d in 1:N]
    for (d, v) in subs
        idx[d] = v
    end
    push!(idx, k)
    Expr(:ref, :y, idx...)
end

# Build bc[i1, ..., (1|end), ..., iN, k] with dimension d fixed to boundary.
function _bc_ref(N::Int, k::Int, d::Int, boundary)
    idx = Any[Symbol("i", dim) for dim in 1:N]
    idx[d] = boundary  # 1 or :end
    push!(idx, k)
    Expr(:ref, :bc, idx...)
end

@generated function differentiate(::Val{solnames}, grid::StateGrid{T1, Ndim, <: NamedTuple{Names}}, y::AbstractArray{T}, icar, bc) where {solnames, T1, Ndim, Names, T}
    statenames = Names

    # Preamble: extract indices and compute grid spacings
    preamble = Expr[]
    for d in 1:Ndim
        id = Symbol("i", d)
        gd = Symbol("grid", d)
        Δm = Symbol("Δx", d, "m")
        Δp = Symbol("Δx", d, "p")
        Δ  = Symbol("Δx", d)
        push!(preamble, :($id = icar[$d]))
        push!(preamble, :($gd = grid.x[$d]))
        push!(preamble, :(@inbounds $Δm = $gd[max($id, 2)] - $gd[max($id - 1, 1)]))
        push!(preamble, :(@inbounds $Δp = $gd[min($id + 1, size(y, $d))] - $gd[min($id, size(y, $d) - 1)]))
        push!(preamble, :($Δ = ($Δm + $Δp) / 2))
    end

    # Build derivative expressions for each solution variable
    expr = Expr[]
    for k in 1:length(solnames)
        solname = solnames[k]

        # Value
        push!(expr, Expr(:(=), solname, _y_ref(Ndim, k)))

        # First derivatives (upwind and downwind) per dimension
        for d in 1:Ndim
            sn = statenames[d]
            id = Symbol("i", d)
            Δp = Symbol("Δx", d, "p")
            Δm = Symbol("Δx", d, "m")
            y_base = _y_ref(Ndim, k)
            y_fwd  = _y_ref(Ndim, k, d => :($id + 1))
            y_bwd  = _y_ref(Ndim, k, d => :($id - 1))
            bc_hi  = _bc_ref(Ndim, k, d, :end)
            bc_lo  = _bc_ref(Ndim, k, d, 1)
            push!(expr, Expr(:(=), Symbol(solname, sn, :_up),
                :(($id < size(y, $d)) ? ($y_fwd - $y_base) / $Δp : convert($T, $bc_hi))))
            push!(expr, Expr(:(=), Symbol(solname, sn, :_down),
                :(($id > 1) ? ($y_base - $y_bwd) / $Δm : convert($T, $bc_lo))))
        end

        # Second derivatives per dimension (with boundary conditions)
        for d in 1:Ndim
            sn = statenames[d]
            id = Symbol("i", d)
            Δp = Symbol("Δx", d, "p")
            Δm = Symbol("Δx", d, "m")
            Δ  = Symbol("Δx", d)
            y_fwd  = _y_ref(Ndim, k, d => :($id + 1))
            y_base = _y_ref(Ndim, k)
            y_bwd  = _y_ref(Ndim, k, d => :($id - 1))
            # left boundary (id == 1)
            y_at_2  = _y_ref(Ndim, k, d => 2)
            y_at_1  = _y_ref(Ndim, k, d => 1)
            bc_lo   = _bc_ref(Ndim, k, d, 1)
            # right boundary (id == end)
            y_at_end   = _y_ref(Ndim, k, d => :(size(y, $d)))
            y_at_endm1 = _y_ref(Ndim, k, d => :(size(y, $d) - 1))
            bc_hi      = _bc_ref(Ndim, k, d, :end)
            interior = :($y_fwd / ($Δp * $Δ) + $y_bwd / ($Δm * $Δ) - 2 * $y_base / ($Δp * $Δm))
            left_bc  = :($y_at_2 / ($Δp * $Δ) + ($y_at_1 - $bc_lo * $Δm) / ($Δm * $Δ) - 2 * $y_at_1 / ($Δp * $Δm))
            right_bc = :(($y_at_end + $bc_hi * $Δp) / ($Δp * $Δ) + $y_at_endm1 / ($Δm * $Δ) - 2 * $y_at_end / ($Δp * $Δm))
            push!(expr, Expr(:(=), Symbol(solname, sn, sn),
                :((1 < $id < size(y, $d)) ? $interior : (($id == 1) ? $left_bc : $right_bc))))
        end

        # Cross derivatives ∂²y/∂x_{d1}∂x_{d2} for each pair of dimensions.
        #
        # We return THREE discretizations, mirroring how first derivatives expose an
        # `_up`/`_down` pair that the model selects on the sign of the drift:
        #
        #   * `y<s1><s2>`            symmetric central difference. 2nd-order accurate but
        #                           NOT monotone — it places negative weights on the
        #                           off-diagonal of the discrete generator as soon as the
        #                           cross term is present, breaking the discrete maximum
        #                           principle for strongly correlated states (oscillations,
        #                           Newton convergence trouble).
        #   * `y<s1><s2>_up`        directional stencil on the MAIN diagonal
        #                           (i+1,j+1)/(i-1,j-1). Use when the cross-diffusion
        #                           coefficient a₁₂ = cov(dx_{d1}, dx_{d2})/dt ≥ 0.
        #   * `y<s1><s2>_down`      directional stencil on the ANTI diagonal
        #                           (i+1,j-1)/(i-1,j+1). Use when a₁₂ < 0.
        #
        # In the model, pick per grid point just like first-order upwinding — using the
        # sign of whatever multiplies the cross derivative in the drift, e.g.
        #       yxz = (σx * σz >= 0) ? yxz_up : yxz_down
        # Each directional stencil is the average of a one-sided forward–forward and
        # backward–backward corner difference, so it is only 1st-order accurate, but its
        # axial terms combine with the (central) own-second-derivatives `y<s><s>` so that
        # every off-diagonal generator weight stays ≥ 0 when the diffusion matrix is
        # diagonally dominant (aᵢᵢ ≥ |a₁₂|) — a monotone, convergent scheme. All three
        # stencils stay within the ±1 neighborhood, so the Jacobian sparsity is unchanged.
        for d1 in 1:Ndim, d2 in (d1+1):Ndim
            sn1 = statenames[d1]
            sn2 = statenames[d2]
            id1 = Symbol("i", d1)
            id2 = Symbol("i", d2)
            Δ1  = Symbol("Δx", d1)
            Δ2  = Symbol("Δx", d2)
            Δ1p = Symbol("Δx", d1, "p")
            Δ1m = Symbol("Δx", d1, "m")
            Δ2p = Symbol("Δx", d2, "p")
            Δ2m = Symbol("Δx", d2, "m")
            y_00 = _y_ref(Ndim, k)
            y_p0 = _y_ref(Ndim, k, d1 => :(min($id1 + 1, size(y, $d1))))
            y_m0 = _y_ref(Ndim, k, d1 => :(max($id1 - 1, 1)))
            y_0p = _y_ref(Ndim, k, d2 => :(min($id2 + 1, size(y, $d2))))
            y_0m = _y_ref(Ndim, k, d2 => :(max($id2 - 1, 1)))
            y_pp = _y_ref(Ndim, k, d1 => :(min($id1 + 1, size(y, $d1))), d2 => :(min($id2 + 1, size(y, $d2))))
            y_pm = _y_ref(Ndim, k, d1 => :(min($id1 + 1, size(y, $d1))), d2 => :(max($id2 - 1, 1)))
            y_mp = _y_ref(Ndim, k, d1 => :(max($id1 - 1, 1)),            d2 => :(min($id2 + 1, size(y, $d2))))
            y_mm = _y_ref(Ndim, k, d1 => :(max($id1 - 1, 1)),            d2 => :(max($id2 - 1, 1)))
            # central: symmetric, 2nd-order, not monotone
            push!(expr, Expr(:(=), Symbol(solname, sn1, sn2),
                :(($y_pp - $y_pm - $y_mp + $y_mm) / (4 * $Δ1 * $Δ2))))
            # upwind for a₁₂ ≥ 0: average of forward–forward and backward–backward corners
            push!(expr, Expr(:(=), Symbol(solname, sn1, sn2, :_up),
                :(($y_pp - $y_p0 - $y_0p + $y_00) / (2 * $Δ1p * $Δ2p)
                + ($y_00 - $y_m0 - $y_0m + $y_mm) / (2 * $Δ1m * $Δ2m))))
            # upwind for a₁₂ < 0: average of forward–backward and backward–forward corners
            push!(expr, Expr(:(=), Symbol(solname, sn1, sn2, :_down),
                :(-($y_pm - $y_p0 - $y_0m + $y_00) / (2 * $Δ1p * $Δ2m)
                  -($y_00 - $y_m0 - $y_0p + $y_mp) / (2 * $Δ1m * $Δ2p))))
        end
    end

    quote
        $(Expr(:meta, :inline))
        $(preamble...)
        @inbounds $(Expr(:tuple, expr...))
    end
end
