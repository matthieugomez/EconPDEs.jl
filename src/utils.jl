"""
StrateGrid(x::NamedTuple)

Let x = (k1 = v1, k2 = v2, ...., kN = vN) be a NamedTuple of AbstractVectors for state space, 
StateGrid(x) returns an AbstractArray M such that 
M[i1, i2..., iN] = (k1 = v1[i1], k2 = v2[i2], ...., kN = vN[iN]
"""
struct StateGrid{T, N, C <: NamedTuple} <: AbstractArray{T, N}
    x::C
end
function StateGrid(x::NamedTuple{Names, <: NTuple{N, <: AbstractVector{T}}}) where {Names, N, T}
    StateGrid{T, N, typeof(x)}(x)
end
function Base.eltype(stategrid::StateGrid{T, N, <: NamedTuple{Names, V}}) where {T, N, Names, V}
    NamedTuple{Names, NTuple{N, T}}
end
Base.ndims(stategrid::StateGrid{T, N}) where {T, N} = N
Base.size(stategrid::StateGrid{T, N}) where {T, N} = ntuple(i -> length(stategrid.x[i]), N)
Base.eachindex(stategrid::StateGrid) = CartesianIndices(size(stategrid))
function Base.getindex(stategrid::StateGrid, args::CartesianIndex)
    eltype(stategrid)(ntuple(i -> stategrid.x[i][args[i]], ndims(stategrid)))
end
Base.getindex(stategrid::StateGrid, x::Symbol) = stategrid.x[x]

#========================================================================================

Derive

========================================================================================#
# 1 state variable
@generated function differentiate(::Type{Tsolution}, grid::StateGrid{T1, 1, <: NamedTuple{N}}, y::AbstractArray{T}, icar, bc) where {Tsolution, T1, N, T}
    statename = N[1]
    expr = Expr[]
    for k in 1:length(Tsolution.parameters[1])
        solname = Tsolution.parameters[1][k]
        push!(expr, Expr(:(=), solname, :(y[i, $k])))
        push!(expr, Expr(:(=), Symbol(solname, statename, :_up), :((i < size(y, 1)) ? (y[i+1, $k] - y[i, $k]) / Δxp : convert($T, bc[end, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename, :_down), :((i > 1) ? (y[i, $k] - y[i-1, $k]) / Δxm : convert($T, bc[1, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename, statename), :((1 < i < size(y, 1)) ? (y[i + 1, $k] / (Δxp * Δx) + y[i - 1, $k] / (Δxm * Δx) - 2 * y[i, $k] / (Δxp * Δxm)) : ((i == 1) ? (y[2, $k] / (Δxp * Δx) + (y[1, $k] - bc[1, $k] * Δxm) / (Δxm * Δx) - 2 * y[1, $k] / (Δxp * Δxm)) : ((y[end, $k] + bc[end, $k] * Δxp) / (Δxp * Δx) + y[end - 1, $k] / (Δxm * Δx) - 2 * y[end, $k] / (Δxp * Δxm))))))
    end
    quote
        $(Expr(:meta, :inline))
        i = icar[1]
        grida = grid.x[1]
        @inbounds Δxm = grida[max(i, 2)] - grida[max(i-1, 1)]
        @inbounds Δxp = grida[min(i+1, size(y, 1))] - grida[min(i, size(y, 1) - 1)]
        Δx = (Δxm + Δxp) / 2
        @inbounds $(Expr(:tuple, expr...))
    end
end

# 2 state variables
@generated function differentiate(::Type{Tsolution}, grid::StateGrid{T1, 2, <: NamedTuple{N}}, y::AbstractArray{T}, icar, bc) where {Tsolution, T1, N, T}
    statename1 = N[1]
    statename2 = N[2]
    expr = Expr[]
    for k in 1:length(Tsolution.parameters[1])
        solname = Tsolution.parameters[1][k]
        push!(expr, Expr(:(=), solname, :(y[i1, i2, $k])))
        push!(expr, Expr(:(=), Symbol(solname, statename1, :_up), :((i1 < size(y, 1)) ? (y[i1+1, i2, $k] - y[i1, i2, $k]) / Δx1p : convert($T, bc[end, i2, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename1, :_down), :((i1 > 1) ? (y[i1, i2, $k] - y[i1-1, i2, $k]) / Δx1m : convert($T, bc[1, i2, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename2, :_up), :((i2 < size(y, 2)) ? (y[i1, i2+1, $k] - y[i1, i2, $k]) / Δx2p : convert($T, bc[i1, end, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename2, :_down),  :((i2 > 1) ? (y[i1, i2, $k] - y[i1, i2-1, $k]) / Δx2m : convert($T, bc[i1, 1, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename1), :((1 < i1 < size(y, 1)) ? (y[i1 + 1, i2, $k] / (Δx1p * Δx1) + y[i1 - 1, i2, $k] / (Δx1m * Δx1) - 2 * y[i1, i2, $k] / (Δx1p * Δx1m)) : ((i1 == 1) ? (y[2, i2, $k] / (Δx1p * Δx1) + (y[1, i2, $k] - bc[1, i2, $k] * Δx1m) / (Δx1m * Δx1) - 2 * y[1, i2, $k] / (Δx1p * Δx1m)) : ((y[end, i2, $k] + bc[end, i2, $k] * Δx1p) / (Δx1p * Δx1) + y[end - 1, i2, $k] / (Δx1m * Δx1) - 2 * y[end, i2, $k] / (Δx1p * Δx1m))))))
        push!(expr, Expr(:(=), Symbol(solname, statename2, statename2), :((1 < i2 < size(y, 2)) ? (y[i1, i2 + 1, $k] / (Δx2p * Δx2) + y[i1, i2 - 1, $k] / (Δx2m * Δx2) - 2 * y[i1, i2, $k] / (Δx2p * Δx2m)) : ((i2 == 1) ? (y[i1, 2, $k] / (Δx2p * Δx2) + (y[i1, 1, $k] - bc[i1, 1, $k] * Δx2m) / (Δx2m * Δx2) - 2 * y[i1, 1, $k] / (Δx2p * Δx2m)) : ((y[i1, end, $k] + bc[i1, end, $k] * Δx2p) / (Δx2p * Δx2) + y[i1, end - 1, $k] / (Δx2m * Δx2) - 2 * y[i1, end, $k] / (Δx2p * Δx2m))))))
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename2), :((y[min(i1 + 1, size(y, 1)), min(i2 + 1, size(y, 2)), $k] - y[min(i1 + 1, size(y, 1)), max(i2 - 1, 1), $k] - y[max(i1 - 1, 1), min(i2 + 1, size(y, 2)), $k] + y[max(i1 - 1, 1), max(i2 - 1, 1), $k]) / (4 * Δx1 * Δx2))))
    end
    quote
        $(Expr(:meta, :inline))
        i1, i2 = icar[1], icar[2]
        grid1, grid2 = grid.x[1], grid.x[2]
        @inbounds Δx1m = grid1[max(i1, 2)] - grid1[max(i1-1, 1)]
        @inbounds Δx1p = grid1[min(i1+1, size(y, 1))] - grid1[min(i1, size(y, 1) - 1)]
        Δx1 = (Δx1m + Δx1p) / 2
        @inbounds Δx2m = grid2[max(i2, 2)] - grid2[max(i2-1, 1)]
        @inbounds Δx2p = grid2[min(i2+1, size(y, 2))] - grid2[min(i2, size(y, 2) - 1)]
        Δx2 = (Δx2m + Δx2p) / 2
        @inbounds $(Expr(:tuple, expr...))
    end
end

# 3 state variables
@generated function differentiate(::Type{Tsolution}, grid::StateGrid{T1, 3, <: NamedTuple{N}}, y::AbstractArray{T}, icar, bc) where {Tsolution, T1, N, T}
    statename1 = N[1]
    statename2 = N[2]
    statename3 = N[3]
    expr = Expr[]
    for k in 1:length(Tsolution.parameters[1])
        solname = Tsolution.parameters[1][k]
        push!(expr, Expr(:(=), solname, :(y[i1, i2, i3, $k])))
        push!(expr, Expr(:(=), Symbol(solname, statename1, :_up), :((i1 < size(y, 1)) ? (y[i1+1, i2, i3, $k] - y[i1, i2, i3, $k]) / Δx1p : convert($T, bc[end, i2, i3, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename1, :_down), :((i1 > 1) ? (y[i1, i2, i3, $k] - y[i1-1, i2, i3, $k]) / Δx1m : convert($T, bc[1, i2, i3, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename2, :_up), :((i2 < size(y, 2)) ? (y[i1, i2+1, i3, $k] - y[i1, i2, i3, $k]) / Δx2p : convert($T, bc[i1, end, i3, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename2, :_down),  :((i2 > 1) ? (y[i1, i2, i3, $k] - y[i1, i2-1, i3, $k]) / Δx2m : convert($T, bc[i1, 1, i3, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename3, :_up), :((i3 < size(y, 3)) ? (y[i1, i2, i3+1, $k] - y[i1, i2, i3, $k]) / Δx3p : convert($T, bc[i1, i2, end, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename3, :_down),  :((i3 > 1) ? (y[i1, i2, i3, $k] - y[i1, i2, i3-1, $k]) / Δx3m : convert($T, bc[i1, i2, 1, $k]))))
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename1), :((1 < i1 < size(y, 1)) ? (y[i1 + 1, i2, i3, $k] / (Δx1p * Δx1) + y[i1 - 1, i2, i3, $k] / (Δx1m * Δx1) - 2 * y[i1, i2, i3, $k] / (Δx1p * Δx1m)) : ((i1 == 1) ? (y[2, i2, i3, $k] / (Δx1p * Δx1) + (y[1, i2, i3, $k] - bc[1, i2, $k] * Δx1m) / (Δx1m * Δx1) - 2 * y[1, i2, i3, $k] / (Δx1p * Δx1m)) : ((y[end, i2, i3, $k] + bc[end, i2, i3, $k] * Δx1p) / (Δx1p * Δx1) + y[end - 1, i2, i3, $k] / (Δx1m * Δx1) - 2 * y[end, i2, i3, $k] / (Δx1p * Δx1m))))))
        push!(expr, Expr(:(=), Symbol(solname, statename2, statename2), :((1 < i2 < size(y, 2)) ? (y[i1, i2 + 1, i3, $k] / (Δx2p * Δx2) + y[i1, i2 - 1, i3, $k] / (Δx2m * Δx2) - 2 * y[i1, i2, i3, $k] / (Δx2p * Δx2m)) : ((i2 == 1) ? (y[i1, 2, i3, $k] / (Δx2p * Δx2) + (y[i1, 1, i3, $k] - bc[i1, 1, i3, $k] * Δx2m) / (Δx2m * Δx2) - 2 * y[i1, 1, i3, $k] / (Δx2p * Δx2m)) : ((y[i1, end, i3, $k] + bc[i1, end, i3, $k] * Δx2p) / (Δx2p * Δx2) + y[i1, end - 1, i3, $k] / (Δx2m * Δx2) - 2 * y[i1, end, i3, $k] / (Δx2p * Δx2m))))))
        push!(expr, Expr(:(=), Symbol(solname, statename3, statename3), :((1 < i3 < size(y, 3)) ? (y[i1, i2, i3 + 1, $k] / (Δx3p * Δx3) + y[i1, i2, i3 - 1, $k] / (Δx3m * Δx3) - 2 * y[i1, i2, i3, $k] / (Δx3p * Δx3m)) : ((i3 == 1) ? (y[i1, i2, 2, $k] / (Δx3p * Δx3) + (y[i1, i2, 1, $k] - bc[i1, i2, 1, $k] * Δx3m) / (Δx3m * Δx3) - 2 * y[i1, i2, 1, $k] / (Δx3p * Δx3m)) : ((y[i1, i2, end, $k] + bc[i1, i2, end, $k] * Δx3p) / (Δx3p * Δx3) + y[i1, i2, end - 1, $k] / (Δx3m * Δx3) - 2 * y[i1, i2, end, $k] / (Δx3p * Δx3m))))))
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename2), :((y[min(i1 + 1, size(y, 1)), min(i2 + 1, size(y, 2)), $k] - y[min(i1 + 1, size(y, 1)), max(i2 - 1, 1), $k] - y[max(i1 - 1, 1), min(i2 + 1, size(y, 2)), $k] + y[max(i1 - 1, 1), max(i2 - 1, 1), $k]) / (4 * Δx1 * Δx2))))
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename3), :((y[min(i1 + 1, size(y, 1)), i2, min(i3 + 1, size(y, 3)), $k] - y[min(i1 + 1, size(y, 1)), i2, max(i3 - 1, 1), $k] - y[max(i1 - 1, 1), i2, min(i3 + 1, size(y, 3)), $k] + y[max(i1 - 1, 1), max(i3 - 1, 1), $k]) / (4 * Δx1 * Δx3))))
        push!(expr, Expr(:(=), Symbol(solname, statename2, statename3), :((y[i1, min(i2 + 1, size(y, 2)), min(i3 + 1, size(y, 3)), $k] - y[i1, min(i2 + 1, size(y, 2)), max(i3 - 1, 1), $k] - y[i1, max(i2 - 1, 1), min(i3 + 1, size(y, 3)), $k] + y[i1, max(i2 - 1, 1), max(i3 - 1, 1), $k]) / (4 * Δx2 * Δx3))))
    end
    quote
        $(Expr(:meta, :inline))
        i1, i2, i3 = icar[1], icar[2], icar[3]
        grid1, grid2, grid3 = grid.x[1], grid.x[2], grid.x[3]
        @inbounds Δx1m = grid1[max(i1, 2)] - grid1[max(i1-1, 1)]
        @inbounds Δx1p = grid1[min(i1+1, size(y, 1))] - grid1[min(i1, size(y, 1) - 1)]
        Δx1 = (Δx1m + Δx1p) / 2
        @inbounds Δx2m = grid2[max(i2, 2)] - grid2[max(i2-1, 1)]
        @inbounds Δx2p = grid2[min(i2+1, size(y, 2))] - grid2[min(i2, size(y, 2) - 1)]
        Δx2 = (Δx2m + Δx2p) / 2
        @inbounds Δx3m = grid3[max(i3, 2)] - grid3[max(i3-1, 1)]
        @inbounds Δx3p = grid3[min(i3+1, size(y, 3))] - grid3[min(i3, size(y, 3) - 1)]
        Δx3 = (Δx3m + Δx3p) / 2
        @inbounds $(Expr(:tuple, expr...))
    end
end



# 1 state variable
@generated function differentiate2(::Type{Tsolution}, grid::StateGrid{T1, 1, <: NamedTuple{N}}, y::AbstractArray{T}, bc) where {Tsolution, T1, N, T}
    statename = N[1]
    expr = Expr[]
    for k in 1:length(Tsolution.parameters[1])
        solname = Tsolution.parameters[1][k]
        push!(expr, Expr(:(=), solname, :(v[$k])))
        push!(expr, Expr(:(=), Symbol(solname, statename, :_up), :(va_up[$k])))
        push!(expr, Expr(:(=), Symbol(solname, statename, :_down), :(va_down[$k])))
        push!(expr, Expr(:(=), Symbol(solname, statename, statename), :(vaa[$k])))
    end
    quote
        $(Expr(:meta, :inline))
        K = length(Tsolution.parameters[1])
        v = [zeros(length(grid.x)) for k in 1:K]
        va_up = [zeros(length(grid.x)) for k in 1:K]
        va_up = [zeros(length(grid.x)) for k in 1:K]
        vaa = [zeros(length(grid.x)) for k in 1:K]
        for i in 1:length(grid.x)
            grida = grid.x[1]
            Δxm = grida[max(i, 2)] - grida[max(i-1, 1)]
            Δxp = grida[min(i+1, size(y, 1))] - grida[min(i, size(y, 1) - 1)]
            Δx = (Δxm + Δxp) / 2
            for k in 1:K
                v[k][i] = y[i, k]
                va_up[k][i] = (i < size(y, 1)) ? (y[i+1, k] - y[i, k]) / Δxp : convert($T, bc[end, k])
                va_down[k][i] = ((i > 1) ? (y[i, k] - y[i-1, k]) / Δxm : convert($T, bc[1, k]))
                vaa[k][i] = ((1 < i < size(y, 1)) ? (y[i + 1, k] / (Δxp * Δx) + y[i - 1, k] / (Δxm * Δx) - 2 * y[i, k] / (Δxp * Δxm)) : ((i == 1) ? (y[2, k] / (Δxp * Δx) + (y[1, k] - bc[1, k] * Δxm) / (Δxm * Δx) - 2 * y[1, k] / (Δxp * Δxm)) : ((y[end, k] + bc[end, k] * Δxp) / (Δxp * Δx) + y[end - 1, k] / (Δxm * Δx) - 2 * y[end, k] / (Δxp * Δxm))))
            end
        end
        @inbounds $(Expr(:tuple, expr...))
    end
end