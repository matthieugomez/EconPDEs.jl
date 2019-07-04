
#========================================================================================
Generates first order and second order derivatibes (using upwinding)
========================================================================================#

# Case with 1 state variable
@generated function derive(::Type{Tsolution}, grid::StateGrid{1, Tstate}, y::AbstractArray, icar, bc::Nothing, drift = (0.0,)) where {Tsolution, Tstate}
    N = length(Tsolution.parameters[1])
    statename = Tstate[1]
    expr = Expr[]
    for k in 1:N
        solname = Tsolution.parameters[1][k]
        push!(expr, Expr(:(=), solname, :(y[i, $k])))
        push!(expr, Expr(:(=), Symbol(solname, statename), :(μx >= 0.0 ? (y[min(i + 1, size(y, 1)), $k] - y[i, $k]) / Δxp : (y[i, $k] - y[max(i - 1, 1), $k]) / Δxm)))
        push!(expr, Expr(:(=), Symbol(solname, statename, statename), :(y[min(i + 1, size(y, 1)), $k] / (Δxp * Δx)+ y[max(i - 1, 1), $k] / (Δxm * Δx) - 2 * y[i, $k] / (Δxp * Δxm))))
    end
    out = Expr(:tuple, expr...)
    quote
        $(Expr(:meta, :inline))
        i = icar[1]
        μx = drift[1]
        grida = grid.x[1]
        Δxm = grida[max(i, 2)] - grida[max(i-1, 1)]
        Δxp = grida[min(i+1, size(y, 1))] - grida[min(i, size(y, 1) - 1)]
        Δx = (Δxm + Δxp) / 2
        @inbounds $out
    end
end


# Case with 2 state variables
@generated function derive(::Type{Tsolution}, grid::StateGrid{2, Tstate}, y::AbstractArray, icar, bc::Nothing, drift = (0.0, 0.0)) where {Tsolution, Tstate}
    N = length(Tsolution.parameters[1])
    statename1 = Tstate[1]
    statename2 = Tstate[2]
    expr = Expr[]
    for k in 1:N
        solname = Tsolution.parameters[1][k]
        push!(expr, Expr(:(=), solname, :(y[i1, i2, $k])))
        push!(expr, Expr(:(=), Symbol(solname, statename1), :(μx1 >= 0.0 ? (y[min(i1 + 1, size(y, 1)), i2, $k] - y[i1, i2, $k]) / Δx1p : (y[i1, i2, $k] - y[max(i1 - 1, 1), i2, $k]) / Δx1m)))
        push!(expr, Expr(:(=), Symbol(solname, statename2), :(μx2 >= 0.0 ? (y[i1, min(i2 + 1, size(y, 2)), $k] - y[i1, i2, $k]) /Δx2p : (y[i1, i2, $k] - y[i1, max(i2 - 1, 1), $k]) / Δx2m)))
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename1), :(y[min(i1 + 1, size(y, 1)), i2, $k] / (Δx1p * Δx1) + y[max(i1 - 1, 1), i2, $k] / (Δx1m * Δx1) - 2 * y[i1, i2, $k] / (Δx1p * Δx1m))))
        push!(expr, Expr(:(=), Symbol(solname, statename2, statename2), :(y[i1, min(i2 + 1, size(y, 2)), $k] / (Δx2p * Δx2) + y[i1, max(i2 - 1, 1), $k] / (Δx2m * Δx2) - 2 * y[i1, i2, $k] / (Δx2p * Δx2m))))
        # cross deriative. Maybe the sign should depend on σx1_x2. If positively correlated, I take the average of cross forward derivative  and cross backward derivative. If negatively correlated, I take the average of cross forward and backward.
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename2), :((y[min(i1 + 1, size(y, 1)), min(i2 + 1, size(y, 2)), $k] - y[min(i1 + 1, size(y, 1)), max(i2 - 1, 1), $k] - y[max(i1 - 1, 1), min(i2 + 1, size(y, 2)), $k] + y[max(i1 - 1, 1), max(i2 - 1, 1), $k]) / (4 * Δx1 * Δx2))))
    end
    out = Expr(:tuple, expr...)
    quote
        $(Expr(:meta, :inline))
        i1, i2 = icar[1], icar[2]
        μx1, μx2 = drift[1], drift[2]
        grid1, grid2 = grid.x[1], grid.x[2]
        Δx1m = grid1[max(i1, 2)] - grid1[max(i1-1, 1)]
        Δx1p = grid1[min(i1+1, size(y, 1))] - grid1[min(i1, size(y, 1) - 1)]
        Δx1 = (Δx1m + Δx1p) / 2
        Δx2m = grid2[max(i2, 2)] - grid2[max(i2-1, 1)]
        Δx2p = grid2[min(i2+1, size(y, 2))] - grid2[min(i2, size(y, 2) - 1)]
        Δx2 = (Δx2m + Δx2p) / 2
        $out
    end
end


#========================================================================================

Derive with boundary conditions
Ideally one would liek to replace previous lines with bc = 0 but it is slower. Not sure why. Maybe type instability.

========================================================================================#


# Case with 1 state variable
@generated function derive(::Type{Tsolution}, grid::StateGrid{1, Tstate}, y::AbstractArray, icar, bc, drift = (0.0,)) where {Tsolution, Tstate}
    N = length(Tsolution.parameters[1])
    statename = Tstate[1]
    expr = Expr[]
    for k in 1:N
        solname = Tsolution.parameters[1][k]
        push!(expr, Expr(:(=), solname, :(y[i, $k])))
        push!(expr, Expr(:(=), Symbol(solname, statename), :((μx >= 0.0) ? ((i < size(y, 1)) ? (y[i+1, $k] - y[i, $k]) / Δxp : bc[end, $k]) : ((i > 1) ? (y[i, $k] - y[i-1, $k]) / Δxm : bc[1, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename, statename), :((1 < i < size(y, 1)) ? (y[i + 1, $k] / (Δxp * Δx) + y[i - 1, $k] / (Δxm * Δx) - 2 * y[i, $k] / (Δxp * Δxm)) : ((i == 1) ? (y[2, $k] / (Δxp * Δx) + (y[1, $k] - bc[1, $k] * Δxm) / (Δxm * Δx) - 2 * y[1, $k] / (Δxp * Δxm)) : ((y[end, $k] + bc[end, $k] * Δxp) / (Δxp * Δx) + y[end - 1, $k] / (Δxm * Δx) - 2 * y[end, $k] / (Δxp * Δxm))))))
    end
    out = Expr(:tuple, expr...)
    quote
        $(Expr(:meta, :inline))
        i = icar[1]
        μx = drift[1]
        grida = grid.x[1]
        Δxm = grida[max(i, 2)] - grida[max(i-1, 1)]
        Δxp = grida[min(i+1, size(y, 1))] - grida[min(i, size(y, 1) - 1)]
        Δx = (Δxm + Δxp) / 2
        $out
    end
end

# Case with 2 state variables
@generated function derive(::Type{Tsolution}, grid::StateGrid{2, Tstate}, y::AbstractArray, icar, bc, drift = (0.0, 0.0)) where {Tsolution, Tstate}
    N = length(Tsolution.parameters[1])
    statename1 = Tstate[1]
    statename2 = Tstate[2]
    expr = Expr[]
    for k in 1:N
        solname = Tsolution.parameters[1][k]
        push!(expr, Expr(:(=), solname, :(y[i1, i2, $k])))
        push!(expr, Expr(:(=), Symbol(solname, statename1), :((μx1 >= 0.0) ? ((i1 < size(y, 1)) ? (y[i1+1, i2, $k] - y[i1, i2, $k]) / Δx1p : bc[end, i2, $k]) : ((i1 > 1) ? (y[i1, i2, $k] - y[i1-1, i2, $k]) / Δx1m : bc[1, i2, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename2), :((μx2 >= 0.0) ? ((i2 < size(y, 2)) ? (y[i1, i2+1, $k] - y[i1, i2, $k]) / Δx2p : bc[i1, end, $k]) : ((i2 > 1) ? (y[i1, i2, $k] - y[i1, i2-1, $k]) / Δx2m : bc[i1, 1, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename1), :((1 < i1 < size(y, 1)) ? (y[i1 + 1, i2, $k] / (Δx1p * Δx1) + y[i1 - 1, i2, $k] / (Δx1m * Δx1) - 2 * y[i1, i2, $k] / (Δx1p * Δx1m)) : ((i1 == 1) ? (y[2, i2, $k] / (Δx1p * Δx1) + (y[1, i2, $k] - bc[1, i2, $k] * Δx1m) / (Δx1m * Δx1) - 2 * y[1, i2, $k] / (Δx1p * Δx1m)) : ((y[end, i2, $k] + bc[end, i2, $k] * Δx1p) / (Δx1p * Δx1) + y[end - 1, i2, $k] / (Δx1m * Δx1) - 2 * y[end, i2, $k] / (Δx1p * Δx1m))))))
        push!(expr, Expr(:(=), Symbol(solname, statename2, statename2), :((1 < i2 < size(y, 2)) ? (y[i1, i2 + 1, $k] / (Δx2p * Δx2) + y[i1, i2 - 1, $k] / (Δx2m * Δx2) - 2 * y[i1, i2, $k] / (Δx2p * Δx2m)) : ((i2 == 1) ? (y[i1, 2, $k] / (Δx2p * Δx2) + (y[i1, 1, $k] - bc[i1, 1, $k] * Δx2m) / (Δx2m * Δx2) - 2 * y[i1, 1, $k] / (Δx2p * Δx2m)) : ((y[i1, end, $k] + bc[i1, end, $k] * Δx2p) / (Δx2p * Δx2) + y[i1, end - 1, $k] / (Δx2m * Δx2) - 2 * y[i1, end, $k] / (Δx2p * Δx2m))))))
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename2), :((y[min(i1 + 1, size(y, 1)), min(i2 + 1, size(y, 2)), $k] - y[min(i1 + 1, size(y, 1)), max(i2 - 1, 1), $k] - y[max(i1 - 1, 1), min(i2 + 1, size(y, 2)), $k] + y[max(i1 - 1, 1), max(i2 - 1, 1), $k]) / (4 * Δx1 * Δx2))))
    end
    out = Expr(:tuple, expr...)
    quote
        $(Expr(:meta, :inline))
        i1, i2 = icar[1], icar[2]
        μx1, μx2 = drift[1], drift[2]
        grid1, grid2 = grid.x[1], grid.x[2]
        Δx1m = grid1[max(i1, 2)] - grid1[max(i1-1, 1)]
        Δx1p = grid1[min(i1+1, size(y, 1))] - grid1[min(i1, size(y, 1) - 1)]
        Δx1 = (Δx1m + Δx1p) / 2
        Δx2m = grid2[max(i2, 2)] - grid2[max(i2-1, 1)]
        Δx2p = grid2[min(i2+1, size(y, 2))] - grid2[min(i2, size(y, 2) - 1)]
        Δx2 = (Δx2m + Δx2p) / 2
        $out
    end
end
