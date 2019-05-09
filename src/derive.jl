
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
        push!(expr, Expr(:(=), Symbol(solname, statename), :(μx >= 0.0 ? (y[min(i + 1, size(y, 1)), $k] - y[i, $k]) * invΔxp[i] : (y[i, $k] - y[max(i - 1, 1), $k]) * invΔxm[i])))
        push!(expr, Expr(:(=), Symbol(solname, statename, statename), :(y[min(i + 1, size(y, 1)), $k] * invΔxp[i] * invΔx[i] + y[max(i - 1, 1), $k] * invΔxm[i] * invΔx[i] - 2 * y[i, $k] * invΔxp[i] * invΔxm[i])))
    end
    out = Expr(:tuple, expr...)
    quote
        $(Expr(:meta, :inline))
        i = icar[1]
        μx = drift[1]
        invΔx = grid.invΔx[1]
        invΔxm = grid.invΔxm[1]
        invΔxp = grid.invΔxp[1]
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
        push!(expr, Expr(:(=), Symbol(solname, statename1), :(μx1 >= 0.0 ? (y[min(i1 + 1, size(y, 1)), i2, $k] - y[i1, i2, $k]) * invΔx1p[i1] : (y[i1, i2, $k] - y[max(i1 - 1, 1), i2, $k]) * invΔx1m[i1])))
        push!(expr, Expr(:(=), Symbol(solname, statename2), :(μx2 >= 0.0 ? (y[i1, min(i2 + 1, size(y, 2)), $k] - y[i1, i2, $k]) * invΔx2p[i2] : (y[i1, i2, $k] - y[i1, max(i2 - 1, 1), $k]) * invΔx2m[i2])))
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename1), :(y[min(i1 + 1, size(y, 1)), i2, $k] * invΔx1p[i1] * invΔx1[i1] + y[max(i1 - 1, 1), i2, $k] * invΔx1m[i1] * invΔx1[i1] - 2 * y[i1, i2, $k] * invΔx1p[i1] * invΔx1m[i1])))
        push!(expr, Expr(:(=), Symbol(solname, statename2, statename2), :(y[i1, min(i2 + 1, size(y, 2)), $k] * invΔx2p[i2] * invΔx2[i2] + y[i1, max(i2 - 1, 1), $k] * invΔx2m[i2] * invΔx2[i2] - 2 * y[i1, i2, $k] * invΔx2p[i2] * invΔx2m[i2])))
        # cross deriative. Maybe the sign should depend on σx1_x2. If positively correlated, I take the average of cross forward derivative  and cross backward derivative. If negatively correlated, I take the average of cross forward and backward.
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename2), :((y[min(i1 + 1, size(y, 1)), min(i2 + 1, size(y, 2)), $k] - y[min(i1 + 1, size(y, 1)), max(i2 - 1, 1), $k] - y[max(i1 - 1, 1), min(i2 + 1, size(y, 2)), $k] + y[max(i1 - 1, 1), max(i2 - 1, 1), $k]) * invΔx1[i1] * invΔx2[i2] / 4)))
    end
    out = Expr(:tuple, expr...)
    quote
        $(Expr(:meta, :inline))
        i1, i2 = icar[1], icar[2]
        μx1, μx2 = drift[1], drift[2]
        invΔx1m, invΔx2m = grid.invΔxm[1], grid.invΔxm[2]
        invΔx1p, invΔx2p = grid.invΔxp[1], grid.invΔxp[2]
        invΔx1, invΔx2 = grid.invΔx[1], grid.invΔx[2]
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
        push!(expr, Expr(:(=), Symbol(solname, statename), :((μx >= 0.0) ? ((i < size(y, 1)) ? (y[i+1, $k] - y[i, $k]) * invΔxp[i] : bc[end, $k]) : ((i > 1) ? (y[i, $k] - y[i-1, $k]) * invΔxm[i] : bc[1, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename, statename), :((1 < i < size(y, 1)) ? (y[i + 1, $k] * invΔxp[i] * invΔx[i] + y[i - 1, $k] * invΔxm[i] * invΔx[i] - 2 * y[i, $k] * invΔxp[i] * invΔxm[i]) : ((i == 1) ? (y[2, $k] * invΔxp[1] * invΔx[1] + (y[1, $k] - bc[1, $k] / invΔxm[1]) * invΔxm[1] * invΔx[1] - 2 * y[1, $k] * invΔxp[1] * invΔxm[1]) : ((y[end, $k] + bc[end, $k] / invΔxp[end]) * invΔxp[end] * invΔx[end] + y[end - 1, $k] * invΔxm[end] * invΔx[end] - 2 * y[end, $k] * invΔxp[end] * invΔxm[end])))))
    end
    out = Expr(:tuple, expr...)
    quote
        $(Expr(:meta, :inline))
        i = icar[1]
        μx = drift[1]
        invΔx = grid.invΔx[1]
        invΔxm = grid.invΔxm[1]
        invΔxp = grid.invΔxp[1]
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
        push!(expr, Expr(:(=), Symbol(solname, statename1), :((μx1 >= 0.0) ? ((i1 < size(y, 1)) ? (y[i1+1, i2, $k] - y[i1, i2, $k]) * invΔx1p[i1] : bc[end, i2, $k]) : ((i1 > 1) ? (y[i1, i2, $k] - y[i1-1, i2, $k]) * invΔx1m[i1] : bc[1, i2, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename2), :((μx2 >= 0.0) ? ((i2 < size(y, 2)) ? (y[i1, i2+1, $k] - y[i1, i2, $k]) * invΔx2p[i2] : bc[i1, end, $k]) : ((i2 > 1) ? (y[i1, i2, $k] - y[i1, i2-1, $k]) * invΔx2m[i2] : bc[i1, 1, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename1), :((1 < i1 < size(y, 1)) ? (y[i1 + 1, i2, $k] * invΔx1p[i1] * invΔx1[i1] + y[i1 - 1, i2, $k] * invΔx1m[i1] * invΔx1[i1] - 2 * y[i1, i2, $k] * invΔx1p[i1] * invΔx1m[i1]) : ((i1 == 1) ? (y[2, i2, $k] * invΔx1p[1] * invΔx1[1] + (y[1, i2, $k] - bc[1, i2, $k] / invΔx1m[1]) * invΔx1m[1] * invΔx1[1] - 2 * y[1, i2, $k] * invΔx1p[1] * invΔx1m[1]) : ((y[end, i2, $k] + bc[end, i2, $k] / invΔx1p[end]) * invΔx1p[end] * invΔx1[end] + y[end - 1, i2, $k] * invΔx1m[end] * invΔx1[end] - 2 * y[end, i2, $k] * invΔx1p[end] * invΔx1m[end])))))
        push!(expr, Expr(:(=), Symbol(solname, statename2, statename2), :((1 < i2 < size(y, 2)) ? (y[i1, i2 + 1, $k] * invΔx2p[i2] * invΔx2[i2] + y[i1, i2 - 1, $k] * invΔx2m[i2] * invΔx2[i2] - 2 * y[i1, i2, $k] * invΔx2p[i2] * invΔx2m[i2]) : ((i2 == 1) ? (y[i1, 2, $k] * invΔx2p[1] * invΔx2[1] + (y[i1, 1, $k] - bc[i1, 1, $k] / invΔx2m[1]) * invΔx2m[1] * invΔx2[1] - 2 * y[i1, 1, $k] * invΔx2p[1] * invΔx2m[1]) : ((y[i1, end, $k] + bc[i1, end, $k] / invΔx2p[end]) * invΔx2p[end] * invΔx2[end] + y[i1, end - 1, $k] * invΔx2m[end] * invΔx2[end] - 2 * y[i1, end, $k] * invΔx2p[end] * invΔx2m[end])))))
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename2), :((y[min(i1 + 1, size(y, 1)), min(i2 + 1, size(y, 2)), $k] - y[min(i1 + 1, size(y, 1)), max(i2 - 1, 1), $k] - y[max(i1 - 1, 1), min(i2 + 1, size(y, 2)), $k] + y[max(i1 - 1, 1), max(i2 - 1, 1), $k]) * invΔx1[i1] * invΔx2[i2] / 4)))
    end
    out = Expr(:tuple, expr...)
    quote
        $(Expr(:meta, :inline))
        i1, i2 = icar[1], icar[2]
        μx1, μx2 = drift[1], drift[2]
        invΔx1m, invΔx2m = grid.invΔxm[1], grid.invΔxm[2]
        invΔx1p, invΔx2p = grid.invΔxp[1], grid.invΔxp[2]
        invΔx1, invΔx2 = grid.invΔx[1], grid.invΔx[2]
        $out
    end
end
