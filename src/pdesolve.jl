#========================================================================================

Type State Grid

========================================================================================#
struct StateGrid{N, V}
    x::NTuple{N, Vector{Float64}}
end

StateGrid(x) = StateGrid{length(x), tuple(keys(x)...)}(tuple(values(x)...))
Base.size(stategrid::StateGrid) = map(length, stategrid.x)
Base.ndims(stategrid::StateGrid{N, V}) where {N, V} = N
Base.eachindex(stategrid::StateGrid) = CartesianIndices(size(stategrid))
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
# 1 state variable
@generated function derive(::Type{Tsolution}, grid::StateGrid{1, Tstate}, y::AbstractArray{T}, icar, bc, drift = (0.0,)) where {Tsolution, Tstate, T}
    N = length(Tsolution.parameters[1])
    statename = Tstate[1]
    expr = Expr[]
    for k in 1:N
        solname = Tsolution.parameters[1][k]
        push!(expr, Expr(:(=), solname, :(y[i, $k])))
        push!(expr, Expr(:(=), Symbol(solname, statename), :((μx >= 0.0) ? ((i < size(y, 1)) ? (y[i+1, $k] - y[i, $k]) / Δxp : convert($T, bc[end, $k])) : ((i > 1) ? (y[i, $k] - y[i-1, $k]) / Δxm : convert($T, bc[1, $k])))))
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

# 2 state variables
@generated function derive(::Type{Tsolution}, grid::StateGrid{2, Tstate}, y::AbstractArray{T}, icar, bc, drift = (0.0, 0.0)) where {Tsolution, Tstate, T}
    N = length(Tsolution.parameters[1])
    statename1 = Tstate[1]
    statename2 = Tstate[2]
    expr = Expr[]
    for k in 1:N
        solname = Tsolution.parameters[1][k]
        push!(expr, Expr(:(=), solname, :(y[i1, i2, $k])))
        push!(expr, Expr(:(=), Symbol(solname, statename1), :((μx1 >= 0.0) ? ((i1 < size(y, 1)) ? (y[i1+1, i2, $k] - y[i1, i2, $k]) / Δx1p : convert($T, bc[end, i2, $k])) : ((i1 > 1) ? (y[i1, i2, $k] - y[i1-1, i2, $k]) / Δx1m : convert($T, bc[1, i2, $k])))))
        push!(expr, Expr(:(=), Symbol(solname, statename2), :((μx2 >= 0.0) ? ((i2 < size(y, 2)) ? (y[i1, i2+1, $k] - y[i1, i2, $k]) / Δx2p : convert($T, bc[i1, end, $k])) : ((i2 > 1) ? (y[i1, i2, $k] - y[i1, i2-1, $k]) / Δx2m : convert($T, bc[i1, 1, $k])))))
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

#========================================================================================

Define function F!(ydot, y)

========================================================================================#

function hjb!(apm, stategrid::StateGrid, Tsolution, ydot_M::AbstractArray, y_M::AbstractArray, bc_M::AbstractArray)
    for i in eachindex(stategrid)
        solution = derive(Tsolution, stategrid, y_M, i, bc_M)
        outi, drifti, = apm(stategrid[i], solution)
        solution = derive(Tsolution, stategrid, y_M, i, bc_M, drifti)
        outi, drifti, = apm(stategrid[i], solution)
        isa(outi, Number) && @error "The pde function must returns a tuple of tuples, not a tuple of numbers"
        _setindex!(ydot_M, outi, i)
    end
    return ydot_M
end

@generated function _setindex!(ydot_M::AbstractArray, outi::NTuple{N, T}, i::CartesianIndex) where {N, T}
    quote
         $(Expr(:meta, :inline))
         $(Expr(:block, [:(setindex!(ydot_M, outi[$k], i, $k)) for k in 1:N]...))
    end
end

# create hjb! that accepts and returns AbstractVector rather than AbstractArrays
function hjb!(apm, stategrid::StateGrid, Tsolution, ydot::AbstractVector, y::AbstractVector, bc_M::AbstractArray, ysize::NTuple)
    y_M = reshape(y, ysize...)
    ydot_M = reshape(ydot, ysize...)
    vec(hjb!(apm, stategrid, Tsolution, ydot_M, y_M, bc_M))
end


#========================================================================================

Solve the PDE on a given time grid

========================================================================================#

function pdesolve(apm, grid::OrderedDict, y0::OrderedDict, τs::AbstractVector; is_algebraic = OrderedDict(k => false for k in keys(y0)), bc = nothing, kwargs...)
    Tsolution = Type{tuple(keys(y0)...)}
    stategrid = StateGrid(grid)
    all(length.(values(y0)) .== prod(size(stategrid))) || throw("The length of initial solution does not equal the length of the state space")
    is_algebraic = OrderedDict(k => fill(is_algebraic[k], size(y0[k])) for k in keys(y0))
    y = OrderedDict{Symbol, Array{Float64, ndims(stategrid) + 1}}(x => Array{Float64}(undef, (size(stategrid)..., length(τs))) for x in keys(y0))

    issorted(reverse(τs)) || throw("The set of times must be decreasing.")

    # convert to Matrix
    y0_M = _Array(y0)
    ysize = size(y0_M)
    is_algebraic_M = _Array(is_algebraic)
    bc_M = _Array_bc(bc, y0_M, y0, grid)

    # prepare dict
    apm_onestep = hasmethod(apm, Tuple{NamedTuple, NamedTuple, Number}) ? (state, grid) -> apm(state, grid, τs[1]) : apm
    a_keys = get_keys(apm_onestep, stategrid, Tsolution, y0_M, bc_M)
    a = nothing
    if a_keys !== nothing
       a = OrderedDict{Symbol, Array{Float64, ndims(stategrid) + 1}}(a_key => Array{Float64}(undef, size(stategrid)..., length(τs)) for a_key in a_keys)
    end

    # create sparsity
    J0c = sparsity_jac(stategrid, y0)

    # iterate on time
    y_M = y0_M
    distance = 0.0
    for iτ in 1:length(τs)
        a_keys !== nothing && _setindex!(a, iτ, localize(apm, τs[iτ]), stategrid, Tsolution, y_M, bc_M)
        _setindex!(y, iτ, y_M)
        if iτ < length(τs)
            y_M, newdistance = implicit_timestep((ydot, y) -> hjb!(localize(apm, τs[iτ]), stategrid, Tsolution, ydot, y, bc_M, ysize), vec(y_M), τs[iτ] - τs[iτ+1]; is_algebraic = vec(is_algebraic_M), verbose = false, J0c = J0c, kwargs...)
            y_M = reshape(y_M, ysize...)
        end
    end
    return y, a, distance
end

_Array(y0) = cat(collect.(values(y0))...; dims = ndims(first(values(y0))) + 1)

function _Array_bc(bc, y0_M, y0, grid)
    bc_M = zero(y0_M)
    k = 0
    for yname in keys(y0)
        k += 1
        keys_grid = collect(keys(grid))
        if length(keys_grid) == 1
            bc_M[1, k], bc_M[end, k] = bc[Symbol(yname, keys_grid[1])]
        elseif length(keys_grid) == 2
            bc_M[1, :,  k],  bc_M[end, :, k] = bc[Symbol(yname, keys_grid[1])]
            bc_M[:, 1,  k], bc_M[:, end,  k] = bc[Symbol(yname, keys_grid[2])]
        end
    end
    return bc_M
end
_Array_bc(::Nothing, y0_M, y0, grid) = zero(y0_M)


function get_keys(apm, stategrid::StateGrid, Tsolution, y_M::AbstractArray, bc_M::AbstractArray)
    i0 = first(eachindex(stategrid))
    solution = derive(Tsolution, stategrid, y_M, i0, bc_M)
    result = apm(stategrid[i0], solution)
    return (length(result) == 3) ? keys(result[3]) : nothing
end

function localize(apm, τ::Number)
    if hasmethod(apm, Tuple{NamedTuple, NamedTuple, Number})
        (state, grid) -> apm(state, grid, τ)
    elseif hasmethod(apm, Tuple{NamedTuple, NamedTuple})
        apm
    else
        throw("The function encoding the PDE must accept NamedTuples for arguments")
    end
end

function sparsity_jac(stategrid::StateGrid, y0::OrderedDict)
    s = size(stategrid)
    l = prod(s)
    t = (ndims(stategrid), length(y0) > 1)
    if t == (1, 0)
        J = Tridiagonal(ones(l - 1), ones(l), ones(l -1))
        return J, matrix_colors(J)
    elseif t == (2, 0)
        J = BandedBlockBandedMatrix(Ones(l, l), (fill(s[1], s[2]), fill(s[1], s[2])), (1, 1), (1, 1))
        return sparse(J), matrix_colors(J)
    elseif t == (1, 1)
        J = BandedBlockBandedMatrix(Ones(l * length(y0), l * length(y0)), (fill(l, length(y0)) ,fill(l, length(y0))), (length(y0) - 1, length(y0) - 1), (1, 1))
        return sparse(J), matrix_colors(J)
    elseif t == (2, 1)
        J = BandedBlockBandedMatrix(Ones(l * length(y0), l * length(y0)), (repeat(fill(s[1], s[2]), outer = length(y0)), repeat(fill(s[1], s[2]), outer = length(y0))), (s[2] * length(y0) - 1, s[2] * length(y0) - 1), (1, 1))
        return sparse(J), matrix_colors(J)
    else
        return nothing, nothing
    end
end

function _setindex!(y::OrderedDict, iτ::Integer, y_M::AbstractArray)
    N = ndims(y_M) - 1
    i = 0
    for v in values(y)
        i += 1
        v[(Colon() for _ in 1:N)..., iτ] = y_M[(Colon() for _ in 1:N)..., i]
    end
end

function _setindex!(a::OrderedDict, iτ::Integer, apm, stategrid::StateGrid, Tsolution, y_M::AbstractArray, bc_M::AbstractArray)
    for i in eachindex(stategrid)
         state = stategrid[i]
         solution = derive(Tsolution, stategrid, y_M, i, bc_M)
         # upwind
         solution = derive(Tsolution, stategrid, y_M, i, bc_M, apm(state, solution)[2])
         outi = apm(state, solution)[3]
         for (k, v) in zip(values(a), values(outi))
            k[i, iτ] = v
         end
     end
 end

 #========================================================================================

 Solve the stationary solution of the PDE

========================================================================================#

 function pdesolve(apm, grid::OrderedDict, y0::OrderedDict; is_algebraic = OrderedDict(k => false for k in keys(y0)), bc = nothing, kwargs...)
     Tsolution = Type{tuple(keys(y0)...)}
     stategrid = StateGrid(grid)
     all(length.(values(y0)) .== prod(size(stategrid))) || throw("The length of initial solution does not equal the length of the state space")
     is_algebraic = OrderedDict(k => fill(is_algebraic[k], size(y0[k])) for k in keys(y0))
     y = OrderedDict{Symbol, Array{Float64, ndims(stategrid)}}(x =>  Array{Float64}(undef, size(stategrid)) for x in keys(y0))

     # Convert to Matrix
     y0_M = _Array(y0)
     ysize = size(y0_M)
     is_algebraic_M = _Array(is_algebraic)
     bc_M = _Array_bc(bc, y0_M, y0, grid)

     # prepare dict
     a_keys = get_keys(apm, stategrid, Tsolution, y0_M, bc_M)
     a = nothing
     if a_keys !== nothing
        a = OrderedDict{Symbol, Array{Float64, ndims(stategrid)}}(a_key => Array{Float64}(undef, size(stategrid)) for a_key in a_keys)
     end

     # create sparsity
     J0c = sparsity_jac(stategrid, y0)

     y_M, distance = finiteschemesolve((ydot, y) -> hjb!(apm, stategrid, Tsolution, ydot, y, bc_M, ysize), vec(y0_M); is_algebraic = vec(is_algebraic_M),  J0c = J0c, kwargs... )
     y_M = reshape(y_M, ysize...)
     _setindex!(y, 1, y_M)
     if a_keys !== nothing
         _setindex!(a, 1, apm, stategrid, Tsolution, y_M, bc_M)
         a = merge(y, a)
     end
     return y, a, distance
 end


