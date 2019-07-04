
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
    is_algebraic = OrderedDict(k => fill(is_algebraic[k], size(y0[k])) for k in keys(y0))

    y0_M = _Matrix(y0)
    bc_M = _Matrix_bc(bc, y0_M, y0, grid)
    is_algebraic_M = _Matrix(is_algebraic)
    ysize = size(y0_M)
    function F!(ydot, y)
        y = reshape(y, ysize...)
        ydot = reshape(ydot, ysize...)
        vec(hjb!(apm, stategrid, Tsolution, ydot, y, bc_M))
    end
    y_M, distance = finiteschemesolve(F!, vec(y0_M); is_algebraic = vec(is_algebraic_M),  J0c = sparsity_jac(stategrid, y0), kwargs... )
    y_M = reshape(y_M, ysize...)
    y = _Dict(collect(keys(y0)), y_M)
    try
        a = _Dict_result(apm, stategrid, Tsolution, y_M, bc_M)
        merge(y, a)
        return y, a, distance
    catch
        return y, nothing, distance
    end
end

#========================================================================================

Sparsity pattern

========================================================================================#

function sparsity_jac(stategrid, y0)
    s = size(stategrid)
    l = prod(s)
    t = (ndims(stategrid), length(y0) > 1)
    if t == (1, 0)
        J = Tridiagonal(ones(l - 1), ones(l), ones(l -1))
    elseif t == (2, 0)
        J = BandedBlockBandedMatrix(Ones(l, l), (fill(s[1], s[2]), fill(s[1], s[2])), (1, 1), (1, 1))
    elseif t == (1, 1)
        J = BandedBlockBandedMatrix(Ones(l * length(y0), l * length(y0)), (fill(l, length(y0)) ,fill(l, length(y0))), (length(y0) - 1, length(y0) - 1), (1, 1))
    elseif t == (2, 1)
        J = BandedBlockBandedMatrix(Ones(l * length(y0), l * length(y0)), (repeat(fill(s[1], s[2]), outer = length(y0)), repeat(fill(s[1], s[2]), outer = length(y0))), (s[2] * length(y0) - 1, s[2] * length(y0) - 1), (1, 1))
    else
        J = nothing
    end
    return (J === nothing) ? (nothing, nothing) : (sparse(J), matrix_colors(J))
end

#========================================================================================

Dict to Matrix

========================================================================================#

# throw("Naming for spaces and solutions lead to ambiguous derivative names. Use different letters for spaces and for solutions")
function _Matrix(y)
    k1 = collect(keys(y))[1]
    if length(y) == 1
        collect(y[k1])
    else
        cat(collect.(values(y))..., dims = ndims(y[k1]) + 1)
    end
end

function _Matrix_bc(::Nothing, y0_M, y0, grid)
    zero(y0_M)
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