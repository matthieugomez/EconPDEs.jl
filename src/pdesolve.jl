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
        if isa(outi, Number)
            @error "The pde function must returns a tuple of tuples, not a tuple of numbers"
        end
        _setindex!(ydot, outi, i)
    end
    return ydot
end


@generated function _setindex!(ydot::AbstractArray, outi::NTuple{N, T}, i::CartesianIndex) where {N, T}
    quote
         $(Expr(:meta, :inline))
         $(Expr(:block, [:(setindex!(ydot, outi[$k], i, $k)) for k in 1:N]...))
    end
end






#========================================================================================

Solve the stationary solution of the PDE

========================================================================================#
function pdesolve(apm, grid::OrderedDict, y0::OrderedDict; is_algebraic = OrderedDict(k => false for k in keys(y0)), bc = nothing, kwargs...)
    Tsolution = Type{tuple(keys(y0)...)}
    stategrid = StateGrid(grid)
    l = prod(size(stategrid))
    for e in values(y0)
        if length(e) != l
            throw("The length of initial solution $(length(e)) does not equal the length of the state space $l")
        end
    end
    y0_M = _Matrix(y0)
    bc_M = _Matrix_bc(bc, y0_M, y0, grid)
    is_algebraic = OrderedDict(k => fill(is_algebraic[k], size(y0[k])) for k in keys(y0))

    y_M, distance = finiteschemesolve((ydot, y) -> hjb!(apm, stategrid, Tsolution, ydot, y, bc_M), y0_M; is_algebraic = _Matrix(is_algebraic), kwargs...)
    y = _Dict(collect(keys(y0)), y_M)
    try
        a = _Dict_result(apm, stategrid, Tsolution, y_M, bc_M)
        merge(y, a)
        return y, a, distance
    catch
        return y, nothing, distance
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

function _Matrix_bc(bc, y0_M, y0, grid)
    bc_M = zero(y0_M)
    k = 0
    for yname in keys(y0)
        k += 1
        keys_grid = collect(keys(grid))
        if length(keys_grid) == 1
            bc_M[1, k] = bc[Symbol(yname, keys_grid[1])][1]
            bc_M[end, k] = bc[Symbol(yname, keys_grid[1])][2]
        elseif length(keys_grid) == 2
            bc_M[1, :,  k] = bc[Symbol(yname, keys_grid[1])][1]
            bc_M[end, :, k] = bc[Symbol(yname, keys_grid[1])][2]
            bc_M[:, 1,  k] = bc[Symbol(yname, keys_grid[2])][1]
            bc_M[:, end,  k] = bc[Symbol(yname, keys_grid[2])][2]
        end
    end
    return bc_M
end
_Matrix_bc(::Nothing, y0_M, y0, grid) = nothing

function _Dict(k, y_M::AbstractArray)
    if length(k) == 1
        N = ndims(y_M)
        OrderedDict{Symbol, Array{Float64, N}}(k[1] => y_M)
    else
        N = ndims(y_M) - 1
        OrderedDict{Symbol, Array{Float64, N}}(k[i] => y_M[(Colon() for _ in 1:N)..., i] for i in 1:length(k))
    end
end

function _Dict_result(apm, grid::StateGrid{Ngrid, Tstate}, ::Type{Tsolution}, y_M, bc) where {Ngrid, Tstate, Tsolution}
    i0 = iterate(eachindex(grid))[1]
    state = grid[i0]
    solution = derive(Tsolution, grid, y_M, i0, bc)
    x = apm(state, solution)[3]
    A = OrderedDict{Symbol, Array{Float64, Ngrid}}(n => Array{Float64}(undef, size(grid)) for n in keys(x))
    for i in eachindex(grid)
        state = grid[i]
        solution = derive(Tsolution, grid, y_M, i, bc)
        outi = apm(state, solution)[2]
        # upwind
        solution = derive(Tsolution, grid, y_M, i, bc, outi)
        outi = apm(state, solution)[3]
        for (n, v) in pairs(outi)
            A[n][i] = v
        end
    end
    return A
end



#========================================================================================

Solve the PDE on a given time grid

========================================================================================#
function pdesolve(apm, grid::OrderedDict, y0::OrderedDict, τs::AbstractVector; is_algebraic = OrderedDict(k => false for k in keys(y0)), bc = nothing, kwargs...)
    Tsolution = Type{tuple(keys(y0)...)}
    stategrid = StateGrid(grid)
    l = prod(size(stategrid))
    for e in values(y0)
        if length(e) != l
            throw("The length of initial solution $(length(e)) does not equal the length of the state space $l")
        end
    end
    y0_M = _Matrix(y0)
    bc_M = _Matrix_bc(bc, y0_M, y0, grid)
    is_algebraic = OrderedDict(k => fill(is_algebraic[k], size(y0[k])) for k in keys(y0))

    # create storage
    y = _Dict(collect(keys(y0)), (size(y0_M)..., length(τs)))
    a = nothing
    a = _Dict_result((state, grid) -> apm(state, grid, τs[1]), stategrid, Tsolution, y0_M, bc_M, τs)


    y_M = y0_M
    distance = 0.0
    # iterate on time
    for iτ in 1:length(τs)
        apm2 = (state, grid) -> apm(state, grid, τs[iτ])
        _setindex!(y, iτ, y_M)
        if !isa(a, Nothing)
            _setindex!(a, iτ, apm2, stategrid, Tsolution, y_M, bc_M)
        end
        if iτ < length(τs)
            y_M, newdistance = finiteschemesolve((ydot, y) -> hjb!(apm2, stategrid, Tsolution, ydot, y, bc_M), y_M; is_algebraic = _Matrix(is_algebraic), Δ = τs[iτ+1] - τs[iτ], iterations = 1, verbose = false, kwargs...)
        end
    end
    return y, a, distance
end


function _Dict(k, s::Tuple)
    if length(k) == 1
        N = length(s)
        return OrderedDict{Symbol, Array{Float64, N}}(k[1] => Array{Float64, N}(undef, s))
    else
        N = length(s) - 1
        return OrderedDict{Symbol, Array{Float64, N}}(k[i] =>  Array{Float64, N}(undef, s) for i in 1:length(k))
    end
end

function _setindex!(y::OrderedDict, iτ::Integer, y_M::AbstractArray)
    if length(y) == 1
        for k in keys(y)
            y[k][:, iτ] = y_M
        end
    else
        N = ndims(y_M) - 1
        i = 0
        for k in keys(y)
            i += 1
            y[k][:, iτ] = y_M[(Colon() for _ in 1:N)..., i]
        end
    end
end


function _Dict_result(apm, grid::StateGrid{Ngrid, Tstate}, ::Type{Tsolution}, y_M, bc, τs) where {Ngrid, Tstate, Tsolution}
    i0 = iterate(eachindex(grid))[1]
    state = grid[i0]
    solution = derive(Tsolution, grid, y_M, i0, bc)
    x = apm(state, solution)[3]
    return OrderedDict{Symbol, Array{Float64, Ngrid + 1}}(n => Array{Float64}(undef, size(grid)..., length(τs)) for n in keys(x))
end



function _setindex!(a, iτ, apm, grid::StateGrid{Ngrid, Tstate}, ::Type{Tsolution}, y_M, bc) where {Ngrid, Tstate, Tsolution}
    for i in eachindex(grid)
         state = grid[i]
         solution = derive(Tsolution, grid, y_M, i, bc)
         outi = apm(state, solution)[2]
         # upwind
         solution = derive(Tsolution, grid, y_M, i, bc, outi)
         outi = apm(state, solution)[3]
         for (n, v) in pairs(outi)
             a[n][i, iτ] = v
         end
     end
 end


