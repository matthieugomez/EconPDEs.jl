#========================================================================================

PDE solver
1. generates first order and second order derivatibes (using upwinding)
2. call finiteschemesolve on the resulting finite difference scheme

========================================================================================#

#========================================================================================

Type State Grid

========================================================================================#
immutable StateGrid{N}
    x::NTuple{N, Vector{Float64}}
    invΔx::NTuple{N, Vector{Float64}}
    invΔxm::NTuple{N, Vector{Float64}}
    invΔxp::NTuple{N, Vector{Float64}}
    name::NTuple{N, Symbol}
end

function make_Δ(x)
    n = length(x)
    Δxm = zeros(x)
    Δxm[1] = x[2] - x[1]
    for i in 2:n
        Δxm[i] = x[i] - x[i-1]
    end
    Δxp = zeros(x)
    for i in 1:(n-1)
        Δxp[i] = x[i+1] - x[i]
    end
    Δxp[end] = x[n] - x[n-1]
    Δx = (Δxm .+ Δxp) / 2
    return x, 1./Δx, 1./Δxm, 1./Δxp
end

function StateGrid(x)
    names = tuple(keys(x)...)
    StateGrid{length(x)}(
        map(i -> map(v -> make_Δ(v)[i], tuple(values(x)...)), 1:4)...,
        names
    )
end
Base.size(grid::StateGrid) = map(length, grid.x)
Base.eachindex(grid::StateGrid) = CartesianRange(size(grid))
@generated function Base.getindex(grid::StateGrid{N}, ::Type{T}, args::CartesianIndex) where {T, N}
    quote
        $(Expr(:meta, :inline))
        $(Expr(:call, T, [:(getindex(grid.x[$i], args[$i])) for i in 1:N]...))
    end
end
Base.getindex(grid::StateGrid, x::Symbol) = grid.x[find(collect(grid.name) .== x)[1]]

#========================================================================================

Derive

========================================================================================#

# Case with 1 state variable
@generated function derive(::Type{Tsolution}, grid::StateGrid{1}, y::AbstractArray, icar, drift = (0.0,)) where {Tsolution}
    N = div(length(fieldnames(Tsolution)), 3)
    expr = Expr[]
    for k in 1:N
        push!(expr, :(y[i, $k]))
        push!(expr, :(μx >= 0.0 ? (y[min(i + 1, size(y, 1)), $k] - y[i, $k]) * invΔxp[i] : (y[i, $k] - y[max(i - 1, 1), $k]) * invΔxm[i]))
        push!(expr, :(y[min(i + 1, size(y, 1)), $k] * invΔxp[i] * invΔx[i] + y[max(i - 1, 1), $k] * invΔxm[i] * invΔx[i] - 2 * y[i, $k] * invΔxp[i] * invΔxm[i]))
    end
    out = Expr(:call, Tsolution, expr...)
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
@generated function derive(::Type{Tsolution}, grid::StateGrid{2}, y::AbstractArray, icar, drift = (0.0, 0.0)) where {Tsolution}
    N = div(length(fieldnames(Tsolution)), 6)
    expr = Expr[]
    for k in 1:N
        push!(expr, :(y[i1, i2, $k]))
        push!(expr, :((y[i1h, i2, $k] - y[i1l, i2, $k]) * invΔx1[i1]))
        push!(expr, :((y[i1, i2h, $k] - y[i1, i2l, $k]) * invΔx2[i2]))
        push!(expr, :((y[min(i1 + 1, size(y, 1)), i2, $k] + y[max(i1 - 1, 1), i2, $k] - 2 * y[i1, i2, $k]) * invΔx1[i1]^2))
        push!(expr, :((y[i1h, i2h, $k] - y[i1h, i2l, $k] - y[i1l, i2h, $k] + y[i1l, i2l, $k]) * invΔx1[i1] * invΔx2[i2]))
        push!(expr, :((y[i1, min(i2 + 1, size(y, 2)), $k] + y[i1, max(i2 - 1, 1), $k] - 2 * y[i1, i2, $k]) * invΔx2[i2]^2))
    end
    out = Expr(:call, Tsolution, expr...)
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
    N = div(length(fieldnames(Tsolution)), 10)
    expr = Expr[]
    for k in 1:N
        push!(expr, :(y[i1, i2, i3, $k]))
        push!(expr, :((y[i1h, i2, i3, $k] - y[i1l, i2, i3, $k]) * invΔx1[i1]))
        push!(expr, :((y[i1, i2h, i3, $k] - y[i1, i2l, i3, $k]) * invΔx2[i2]))
        push!(expr, :((y[i1, i2, i3h, $k] - y[i1, i2, i3l, $k]) * invΔx3[i3]))
        push!(expr, :((y[min(i1 + 1, size(y, 1)), i2, i3, $k] + y[max(i1 - 1, 1), i2, i3, $k] - 2 * y[i1, i2, i3, $k]) * invΔx1[i1]^2))
        push!(expr, :((y[i1h, i2h, i3, $k] - y[i1h, i2l, i3, $k] - y[i1l, i2h, i3, $k] + y[i1l, i2l, i3, $k]) * invΔx1[i1] * invΔx2[i2]))
        push!(expr, :((y[i1h, i2, i3h, $k] - y[i1h, i2, i3l, $k] - y[i1l, i2, i3h, $k] + y[i1l, i2, i3l, $k]) * invΔx1[i1] * invΔx3[i3]))
        push!(expr, :((y[i1, min(i2 + 1, size(y, 2)), i3, $k] + y[i1, max(i2 - 1, 1), i3, $k] - 2 * y[i1, i2, i3, $k]) * invΔx2[i2]^2))
        push!(expr, :((y[i1, i2h, i3h, $k] - y[i1, i2l, i3l, $k] - y[i1, i2l, i3h, $k] + y[i1, i2l, i3l, $k]) * invΔx2[i2] * invΔx3[i3]))
        push!(expr, :((y[i1, i2, min(i3 + 1, size(y, 3)), $k] + y[i1, i2, max(i3 - 1, 1), $k] - 2 * y[i1, i2, i3, $k]) * invΔx3[i3]^2))
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

##NamedTuple create type in the module
# at parsing time
# 1. creates a type by looking at expression, so something like T_NT_a_b{T1, T2}
# 2. replace @NT by T_NT_ab
_NT(names) =  eval(Expr(:macrocall, Symbol("@NT"), (x for x in names)...))


function hjb!(apm, grid::StateGrid{Ngrid}, ::Type{Tstate}, ::Type{Tsolution}, ydot, y) where {Ngrid, Tstate, Tsolution}
    for i in eachindex(grid)
        state = getindex(grid, Tstate, i)
        solution = derive(Tsolution, grid, y, i)
        drifti = apm(state, solution)[2]
        #upwind
        solution = derive(Tsolution, grid, y, i, drifti)
        outi = apm(state, solution)[1]
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
@inline _setindex!(ydot, outi, i) = setindex!(ydot, outi, i)


function create_dictionary(apm, grid::StateGrid{Ngrid}, ::Type{Tstate}, ::Type{Tsolution}, y) where {Ngrid, Tstate, Tsolution}
    i0 = start(eachindex(grid))
    state = getindex(grid, Tstate, i0)
    solution = derive(Tsolution, grid, y, i0)
    x = apm(state, solution)
    A = nothing
    if length(x) == 3
        names = map(first, x[3])
        A = OrderedDict{Symbol, Array{Float64, Ngrid}}(n => Array{Float64}(size(grid)) for n in names)
        for i in eachindex(grid)
            state = getindex(grid, Tstate, i)
            solution = derive(Tsolution, grid, y, i)
            drifti = apm(state, solution)[2]
            # upwind
            solution = derive(Tsolution, grid, y, i, drifti)
            othersi = apm(state, solution)[3]
            for j in 1:length(names)
                A[names[j]][i] = last(othersi[j])
            end
        end
    end
    return A
end




#========================================================================================

Solve the PDE

========================================================================================#
function pdesolve(apm, grid::OrderedDict, y0::OrderedDict; is_algebraic = Dict(k => false for k in keys(y0)), kwargs...)
    Tstate = _NT(keys(grid))
    Tsolution = _NT(all_symbol(collect(keys(y0)), collect(keys(grid))))
    stategrid = StateGrid(grid)
    is_algebraic = _NT(keys(y0))((fill(is_algebraic[k], size(y0[k])) for k in keys(y0))...)
    y, distance = finiteschemesolve((ydot, y) -> hjb!(apm, stategrid, Tstate, Tsolution, ydot, y), _concatenate(y0); is_algebraic = _concatenate(is_algebraic), kwargs...)
    a = create_dictionary(apm, stategrid, Tstate, Tsolution, y)
    y = _deconcatenate(collect(keys(y0)), y)
    if a != nothing
        a = merge(y, a)
    end
    return y, a, distance
end

function all_symbol(sols, states)
    all_derivatives = vcat((map(x -> Symbol(x...), with_replacement_combinations(states, k)) for k in 0:2)...)
    out = vec([Symbol(a, s) for s in all_derivatives, a in sols])
    if length(out) != length(unique(out))
        throw("Naming for spaces and solutions lead to ambiguous derivative names. Use different letters for spaces and for solutions")
    end
    return out
end

function _concatenate(y)
    k1 = collect(keys(y))[1]
    if length(y) == 1
        y[k1]
    else
        cat(ndims(y[k1]) + 1, values(y)...)
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
