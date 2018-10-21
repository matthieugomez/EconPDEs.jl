#========================================================================================

PDE solver
1. generates first order and second order derivatibes (using upwinding)
2. call finiteschemesolve on the resulting finite difference scheme

========================================================================================#

#========================================================================================

Type State Grid

========================================================================================#
struct StateGrid{N, V}
    x::NTuple{N, Vector{Float64}}
    invΔx::NTuple{N, Vector{Float64}}
    invΔxm::NTuple{N, Vector{Float64}}
    invΔxp::NTuple{N, Vector{Float64}}
end

function make_Δ(x)
    n = length(x)
    Δxm = zero(x)
    Δxm[1] = x[2] - x[1]
    for i in 2:n
        Δxm[i] = x[i] - x[i-1]
    end
    Δxp = zero(x)
    for i in 1:(n-1)
        Δxp[i] = x[i+1] - x[i]
    end
    Δxp[end] = x[n] - x[n-1]
    Δx = (Δxm .+ Δxp) / 2
    return x, 1 ./ Δx, 1 ./ Δxm, 1 ./ Δxp
end

function StateGrid(x)
    StateGrid{length(x), tuple(keys(x)...)}(
        map(i -> map(v -> make_Δ(v)[i], tuple(values(x)...)), 1:4)...)
end
Base.size(grid::StateGrid) = map(length, grid.x)
Base.eachindex(grid::StateGrid) = CartesianIndices(size(grid))
@generated function Base.getindex(grid::StateGrid{N, V}, args::CartesianIndex) where {N, V}
    quote
        $(Expr(:meta, :inline))
        $(Expr(:tuple, [Expr(:(=), V[i], :(grid.x[$i][args[$i]])) for i in 1:N]...))
    end
end
function Base.getindex(grid::StateGrid{N, V}, x::Symbol) where {N, V}
    grid.x[find(collect(V) .== x)[1]]
end

#========================================================================================

Derive

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
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename2), :((y[min(i1 + 1, size(y, 1)), min(i2 + 1, size(y, 2)), $k] - y[min(i1 + 1, size(y, 1)), max(i2 -1, 1), $k] - y[max(i1 -1, 1), min(i2 + 1, size(y, 2)), $k] + y[max(i1 -1, 1), max(i2 -1, 1), $k]) * invΔx1[i1] * invΔx2[i2])))
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

# Case with 3 state variables
#@generated function derive(::Type{Tsolution}, grid::StateGrid{3}, y::AbstractArray, icar, bc, drift = (0.0, 0.0, 0.0)) where {Tsolution}
#    N = length(Tsolution.parameters[1])
#    statename1 = Tstate[1]
#    statename2 = Tstate[2]
#    statename3 = Tstate[3]
#    expr = Expr[]
#    for k in 1:N
#        solname = Tsolution.parameters[1][k]
#        push!(expr, Expr(:(=), solname, :(y[i1, i2, i3, $k])))
#        push!(expr, Expr(:(=), Symbol(solname, statename1), :((y[i1h, i2, i3, $k] - y[i1l, i2, i3, $k]) * invΔx1[i1])))
#        push!(expr, Expr(:(=), Symbol(solname, statename2), :((y[i1, i2h, i3, $k] - y[i1, i2l, i3, $k]) * invΔx2[i2])))
#        push!(expr, Expr(:(=), Symbol(solname, statename3), :((y[i1, i2, i3h, $k] - y[i1, i2, i3l, $k]) * invΔx3[i3])))
#        push!(expr, Expr(:(=), Symbol(solname, statename1, statename1), :((y[min(i1 + 1, size(y, 1)), i2, i3, $k] + y[max(i1 - 1, 1), i2, i3, $k] - 2 * y[i1, i2, i3, $k]) * invΔx1[i1]^2)))
#        push!(expr,Expr(:(=), Symbol(solname, statename1, statename2), :((y[i1h, i2h, i3, $k] - y[i1h, i2l, i3, $k] - y[i1l, i2h, i3, $k] + y[i1l, i2l, i3, $k]) * invΔx1[i1] * invΔx2[i2])))
#        push!(expr, Expr(:(=), Symbol(solname, statename1, statename3), :((y[i1h, i2, i3h, $k] - y[i1h, i2, i3l, $k] - y[i1l, i2, i3h, $k] + y[i1l, i2, i3l, $k]) * invΔx1[i1] * invΔx3[i3])))
#        push!(expr, Expr(:(=), Symbol(solname, statename2, statename2), :((y[i1, min(i2 + 1, size(y, 2)), i3, $k] + y[i1, max(i2 - 1, 1), i3, $k] - 2 * y[i1, i2, i3, $k]) * invΔx2[i2]^2)))
#        push!(expr, Expr(:(=), Symbol(solname, statename2, statename3), :((y[i1, i2h, i3h, $k] - y[i1, i2l, i3l, $k] - y[i1, i2l, i3h, $k] + y[i1, i2l, i3l, $k]) * invΔx2[i2] * invΔx3[i3])))
#        push!(expr, Expr(:(=), Symbol(solname, statename3, statename3), :((y[i1, i2, min(i3 + 1, size(y, 3)), $k] + y[i1, i2, max(i3 - 1, 1), $k] - 2 * y[i1, i2, i3, $k]) * invΔx3[i3]^2)))
#    end
#    out = Expr(:call, Tsolution, expr...)
#    quote
#        $(Expr(:meta, :inline))
#        i1, i2, i3 = icar[1], icar[2], icar[3]
#        μx1, μx2, μx3 = drift[1], drift[2], drift[3]
#        invΔx1, invΔx2, indvΔx3 = grid.invΔx[1], grid.invΔx[2], grid.invΔx[3]
#        if μx1 >= 0.0
#            i1h = min(i1 + 1, size(y, 1))
#            i1l = i1
#        else
#          i1h = i1
#          i1l = max(i1 - 1, 1)
#        end
#        if μx2 >= 0.0
#            i2h = min(i2 + 1, size(y, 2))
#            i2l = i2
#        else
#          i2h = i2
#          i2l = max(i2 - 1, 1)
#        end
#        if μx3 >= 0.0
#            i3h = min(i3 + 1, size(y, 3))
#            i3l = i3
#        else
#          i3h = i3
#          i3l = max(i3 - 1, 1)
#        end
#        $out
#    end
#end



#========================================================================================

Define function F!(ydot, y) to pass to finiteschemesolve

========================================================================================#


function hjb!(apm, grid::StateGrid{Ngrid, Tstate}, Tsolution, ydot, y, bc) where {Ngrid, Tstate}
    for i in eachindex(grid)
        solution = derive(Tsolution, grid, y, i, bc)
        outi = apm(grid[i], solution)[2]
        #upwind
        solution = derive(Tsolution, grid, y, i, bc, outi)
        outi = apm(grid[i], solution)[1]
        _setindex!(ydot, outi, i)
    end
    return ydot
end


@generated function _setindex!(ydot, outi::NTuple{N, T}, i) where {N, T}
    quote
         $(Expr(:meta, :inline))
         $(Expr(:block, [:(setindex!(ydot, outi[$k], i, $k)) for k in 1:N]...))
    end
end

function create_dictionary(apm, grid::StateGrid{Ngrid, Tstate}, ::Type{Tsolution}, y, bc) where {Ngrid, Tstate, Tsolution}
    i0 = iterate(eachindex(grid))[1]
    state = grid[i0]
    solution = derive(Tsolution, grid, y, i0, bc)
    x = apm(state, solution)[3]
    A = OrderedDict{Symbol, Array{Float64, Ngrid}}(n => Array{Float64}(undef, size(grid)) for n in keys(x))
    for i in eachindex(grid)
        state = grid[i]
        solution = derive(Tsolution, grid, y, i, bc)
        outi = apm(state, solution)[2]
        # upwind
        solution = derive(Tsolution, grid, y, i, bc, outi)
        outi = apm(state, solution)[3]
        for (n, v) in pairs(outi)
            A[n][i] = v
        end
    end
    return A
end




#========================================================================================

Solve the PDE

========================================================================================#
function pdesolve(apm, grid::OrderedDict, y0::OrderedDict; is_algebraic = OrderedDict(k => false for k in keys(y0)), bc = OrderedDict(), kwargs...)
    Tsolution = Type{tuple(keys(y0)...)}
    stategrid = StateGrid(grid)
    l = prod(size(stategrid))
    for e in values(y0)
        if length(e) != l
            throw("The length of initial solution $(length(e)) does not equal the length of the state space $l")
        end
    end
    y0_M = _Matrix(y0)
    bc_m = zero(y0_M)
    if !isempty(bc)
        k = 0
        for yname in keys(y0)
            k += 1
            keys_grid = collect(keys(grid))
            if length(keys_grid) == 1
                bc_m[1, k] = bc[Symbol(yname, keys_grid[1])][1]
                bc_m[end, k] = bc[Symbol(yname, keys_grid[1])][2]
            elseif length(keys_grid) == 2
                bc_m[1, :,  k] = bc[Symbol(yname, keys_grid[1])][1]
                bc_m[end, :, k] = bc[Symbol(yname, keys_grid[1])][2]
                bc_m[:, 1,  k] = bc[Symbol(yname, keys_grid[2])][1]
                bc_m[:, 1,  k] = bc[Symbol(yname, keys_grid[2])][2]
            end
        end
    end
    is_algebraic = OrderedDict(k => fill(is_algebraic[k], size(y0[k])) for k in keys(y0))
    y, distance = finiteschemesolve((ydot, y) -> hjb!(apm, stategrid, Tsolution, ydot, y, bc_m), y0_M; is_algebraic = _Matrix(is_algebraic), kwargs...)
    dy = _Dict(collect(keys(y0)), y)
    try
        a = create_dictionary(apm, stategrid, Tsolution, y, bc_m)
        merge(dy, a)
        return dy, a, distance
    catch
        return dy, nothing, distance
    end
end


# throw("Naming for spaces and solutions lead to ambiguous derivative names. Use different letters for spaces and for solutions")


function _Matrix(y)
    k1 = collect(keys(y))[1]
    if length(y) == 1
        y[k1]
    else
        cat(values(y)..., dims = ndims(y[k1]) + 1)
    end
end
function _Dict(k, y)
    if length(k) == 1
        N = ndims(y)
        OrderedDict{Symbol, Array{Float64, N}}(k[1] => y)
    else
        N = ndims(y) - 1
        OrderedDict{Symbol, Array{Float64, N}}(k[i] => y[(Colon() for _ in 1:N)..., i] for i in 1:length(k))
    end
end
