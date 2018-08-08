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
@generated function derive(::Type{Tsolution}, grid::StateGrid{1, Tstate}, y::AbstractArray, icar, drift = (0.0,)) where {Tsolution, Tstate}
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
        $out
    end
end

# Case with 2 state variables
@generated function derive(::Type{Tsolution}, grid::StateGrid{2, Tstate}, y::AbstractArray, icar, drift = (0.0, 0.0)) where {Tsolution, Tstate}
    N = length(Tsolution.parameters[1])
    statename1 = Tstate[1]
    statename2 = Tstate[2]
    expr = Expr[]
    for k in 1:N
        solname = Tsolution.parameters[1][k]
        push!(expr, Expr(:(=), solname, :(y[i1, i2, $k])))
        push!(expr, Expr(:(=), Symbol(solname, statename1), :((y[i1h, i2, $k] - y[i1l, i2, $k]) * invΔx1[i1])))
        push!(expr, Expr(:(=), Symbol(solname, statename2), :((y[i1, i2h, $k] - y[i1, i2l, $k]) * invΔx2[i2])))
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename1), :((y[min(i1 + 1, size(y, 1)), i2, $k] + y[max(i1 - 1, 1), i2, $k] - 2 * y[i1, i2, $k]) * invΔx1[i1]^2)))
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename2), :((y[i1h, i2h, $k] - y[i1h, i2l, $k] - y[i1l, i2h, $k] + y[i1l, i2l, $k]) * invΔx1[i1] * invΔx2[i2])))
        push!(expr, Expr(:(=), Symbol(solname, statename2, statename2), :((y[i1, min(i2 + 1, size(y, 2)), $k] + y[i1, max(i2 - 1, 1), $k] - 2 * y[i1, i2, $k]) * invΔx2[i2]^2)))
    end
    out = Expr(:tuple, expr...)
    quote
        $(Expr(:meta, :inline))
        i1, i2 = icar[1], icar[2]
        μx1, μx2 = drift[1], drift[2]
        invΔx1, invΔx2 = grid.invΔx[1], grid.invΔx[2]
        if μx1 >= 0.0
            i1h = min(i1 + 1, size(y, 1))
            i1l = i1
        else
          i1h = i1
          i1l = max(i1 - 1, 1)
        end
        if μx2 >= 0.0
            i2h = min(i2 + 1, size(y, 2))
            i2l = i2
        else
          i2h = i2
          i2l = max(i2 - 1, 1)
        end
        $out
    end
end

# Case with 3 state variables
@generated function derive(::Type{Tsolution}, grid::StateGrid{3}, y::AbstractArray, icar, drift = (0.0, 0.0, 0.0)) where {Tsolution}
    N = length(Tsolution.parameters[1])
    statename1 = Tstate[1]
    statename2 = Tstate[2]
    statename3 = Tstate[3]
    expr = Expr[]
    for k in 1:N
        solname = Tsolution.parameters[1][k]
        push!(expr, Expr(:(=), solname, :(y[i1, i2, i3, $k])))
        push!(expr, Expr(:(=), Symbol(solname, statename1), :((y[i1h, i2, i3, $k] - y[i1l, i2, i3, $k]) * invΔx1[i1])))
        push!(expr, Expr(:(=), Symbol(solname, statename2), :((y[i1, i2h, i3, $k] - y[i1, i2l, i3, $k]) * invΔx2[i2])))
        push!(expr, Expr(:(=), Symbol(solname, statename3), :((y[i1, i2, i3h, $k] - y[i1, i2, i3l, $k]) * invΔx3[i3])))
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename1), :((y[min(i1 + 1, size(y, 1)), i2, i3, $k] + y[max(i1 - 1, 1), i2, i3, $k] - 2 * y[i1, i2, i3, $k]) * invΔx1[i1]^2)))
        push!(expr,Expr(:(=), Symbol(solname, statename1, statename2), :((y[i1h, i2h, i3, $k] - y[i1h, i2l, i3, $k] - y[i1l, i2h, i3, $k] + y[i1l, i2l, i3, $k]) * invΔx1[i1] * invΔx2[i2])))
        push!(expr, Expr(:(=), Symbol(solname, statename1, statename3), :((y[i1h, i2, i3h, $k] - y[i1h, i2, i3l, $k] - y[i1l, i2, i3h, $k] + y[i1l, i2, i3l, $k]) * invΔx1[i1] * invΔx3[i3])))
        push!(expr, Expr(:(=), Symbol(solname, statename2, statename2), :((y[i1, min(i2 + 1, size(y, 2)), i3, $k] + y[i1, max(i2 - 1, 1), i3, $k] - 2 * y[i1, i2, i3, $k]) * invΔx2[i2]^2)))
        push!(expr, Expr(:(=), Symbol(solname, statename2, statename3), :((y[i1, i2h, i3h, $k] - y[i1, i2l, i3l, $k] - y[i1, i2l, i3h, $k] + y[i1, i2l, i3l, $k]) * invΔx2[i2] * invΔx3[i3])))
        push!(expr, Expr(:(=), Symbol(solname, statename3, statename3), :((y[i1, i2, min(i3 + 1, size(y, 3)), $k] + y[i1, i2, max(i3 - 1, 1), $k] - 2 * y[i1, i2, i3, $k]) * invΔx3[i3]^2)))
    end
    out = Expr(:call, Tsolution, expr...)
    quote
        $(Expr(:meta, :inline))
        i1, i2, i3 = icar[1], icar[2], icar[3]
        μx1, μx2, μx3 = drift[1], drift[2], drift[3]
        invΔx1, invΔx2, indvΔx3 = grid.invΔx[1], grid.invΔx[2], grid.invΔx[3]
        if μx1 >= 0.0
            i1h = min(i1 + 1, size(y, 1))
            i1l = i1
        else
          i1h = i1
          i1l = max(i1 - 1, 1)
        end
        if μx2 >= 0.0
            i2h = min(i2 + 1, size(y, 2))
            i2l = i2
        else
          i2h = i2
          i2l = max(i2 - 1, 1)
        end
        if μx3 >= 0.0
            i3h = min(i3 + 1, size(y, 3))
            i3l = i3
        else
          i3h = i3
          i3l = max(i3 - 1, 1)
        end
        $out
    end
end



#========================================================================================

Define function F!(ydot, y) to pass to finiteschemesolve

========================================================================================#


function hjb!(apm, grid::StateGrid{Ngrid, Tstate}, Tsolution, ydot, y) where {Ngrid, Tstate}
    for i in eachindex(grid)
        solution = derive(Tsolution, grid, y, i)
        outi = apm(grid[i], solution)
        #upwind
        solution = derive(Tsolution, grid, y, i, get_drift(outi, grid))
        outi = apm(grid[i], solution)
        _setindex!(ydot, outi, i, Tsolution)
    end
    return ydot
end

@generated function get_drift(outi, grid::StateGrid{Ngrid, Tstate}) where {Ngrid, Tstate}
    Expr(:tuple, [Expr(:(.), :outi, QuoteNode(Symbol(:μ, Tstate[k]))) for k in 1:length(Tstate)]...)
end

@generated function _setindex!(ydot, outi::NamedTuple{N, T}, i, ::Type{Tsolution}) where {N, T, Tsolution}
    expr = Expr(:block, [:(setindex!(ydot, $(Expr(:(.), :outi, QuoteNode(Symbol(Tsolution.parameters[1][k], :t)))), i, $k)) for k in 1:length(Tsolution.parameters[1])]...)
    quote
        $(Expr(:meta, :inline))
        $expr
    end
end

function create_dictionary(apm, grid::StateGrid{Ngrid, Tstate}, ::Type{Tsolution}, y) where {Ngrid, Tstate, Tsolution}
    i0 = iterate(eachindex(grid))[1]
    state = grid[i0]
    solution = derive(Tsolution, grid, y, i0)
    x = apm(state, solution)
    A = OrderedDict{Symbol, Array{Float64, Ngrid}}(n => Array{Float64}(undef, size(grid)) for n in keys(x))
    for i in eachindex(grid)
        state = grid[i]
        solution = derive(Tsolution, grid, y, i)
        outi = apm(state, solution)
        # upwind
        solution = derive(Tsolution, grid, y, i, get_drift(outi, grid))
        outi = apm(state, solution)
        for (n, v) in pairs(outi)
            A[n][i] = v
        end
    end
    return A
end




#========================================================================================

Solve the PDE

========================================================================================#
function pdesolve(apm, grid::OrderedDict, y0::OrderedDict; is_algebraic = Dict(k => false for k in keys(y0)), kwargs...)
    Tsolution = Type{tuple(keys(y0)...)}
    stategrid = StateGrid(grid)
    is_algebraic = OrderedDict(k => fill(is_algebraic[k], size(y0[k])) for k in keys(y0))
    y, distance = finiteschemesolve((ydot, y) -> hjb!(apm, stategrid, Tsolution, ydot, y), _concatenate(y0); is_algebraic = _concatenate(is_algebraic), kwargs...)
    a = create_dictionary(apm, stategrid, Tsolution, y)
    y = _deconcatenate(collect(keys(y0)), y)
    if a != nothing
        a = merge(y, a)
    end
    return y, a, distance
end


# throw("Naming for spaces and solutions lead to ambiguous derivative names. Use different letters for spaces and for solutions")

function _concatenate(y)
    k1 = collect(keys(y))[1]
    if length(y) == 1
        y[k1]
    else
        cat(values(y)..., dims = ndims(y[k1]) + 1)
    end
end
function _deconcatenate(k, y)
    if length(k) == 1
        N = ndims(y)
        OrderedDict{Symbol, Array{Float64, N}}(k[1] => y)
    else
        N = ndims(y) - 1
        OrderedDict{Symbol, Array{Float64, N}}(k[i] => y[(Colon() for _ in 1:N)..., i] for i in 1:length(k))
    end
end
