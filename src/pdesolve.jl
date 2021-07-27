"""
    pdesolve(f, grid, yend [, τs]; is_algebraic = OrderedDict(k => false for k in keys(yend)), bc = nothing)
    
    solves a system of pdes with an arbitrary number of functions and up to state variable. Denote  F the number of functions to solve for and S the number of state variables.

    ### Arguments
    *  f: is a function (state, y) -> (out, μ) 
    where 
       state is a tuple of real numbers of  length S (for the state values)
       y is a tuple of real numbers of length F * 3 with the functions as well as their first and second derivatives at the state
       out is a tuple of length F with the value of the time derivatives of the functions at the state 
       μ is a tuple of length F with the value for the drift of the state variable at the state
    * grid is an OrderedDict that, for each state, associates an AbstractVector (for the grid)
    * yend is an OrderedDict that, for each function, associates an initial guess (or a terminal value to start the backward iteration with)
    * τs (optional) is a time grid on which to solve the pde on. In this case, yend corresponds to the solution at time τs[end]
"""
function pdesolve(apm, grid::OrderedDict, yend::OrderedDict, τs::Union{Nothing, AbstractVector} = nothing; is_algebraic = OrderedDict(k => false for k in keys(yend)), bc = nothing, kwargs...)
    Tsolution = Type{tuple(keys(yend)...)}
    stategrid = StateGrid(grid)
    all(length.(values(yend)) .== prod(size(stategrid))) || throw("The length of initial guess (e.g. terminal value) does not equal the length of the state space")
    is_algebraic = OrderedDict(k => fill(is_algebraic[k], size(yend[k])) for k in keys(yend))
    if τs isa AbstractVector
        issorted(τs) || throw("The set of times must be increasing.")
        ys = [OrderedDict{Symbol, Array{Float64, ndims(stategrid)}}(x =>  Array{Float64}(undef, size(stategrid)) for x in keys(yend)) for t in τs]
    else
        y = OrderedDict{Symbol, Array{Float64, ndims(stategrid)}}(x =>  Array{Float64}(undef, size(stategrid)) for x in keys(yend))
    end

    # convert to Matrix
    yend_M = _Array(yend)
    ysize = size(yend_M)
    is_algebraic_M = _Array(is_algebraic)
    bc_M = _Array_bc(bc, yend_M, yend, grid)

    # prepare dict
    if τs isa AbstractVector
        apm_onestep = hasmethod(apm, Tuple{NamedTuple, NamedTuple, Number}) ? (state, grid) -> apm(state, grid, τs[end]) : apm
        a_keys = get_keys(apm_onestep, stategrid, Tsolution, yend_M, bc_M)
        as = nothing
        if a_keys !== nothing
           as = [OrderedDict{Symbol, Array{Float64, ndims(stategrid)}}(a_key => Array{Float64}(undef, size(stategrid)) for a_key in a_keys) for a_key in a_keys]
        end
    else
        a_keys = get_keys(apm, stategrid, Tsolution, yend_M, bc_M)
        a = nothing
        if a_keys !== nothing
            a = OrderedDict{Symbol, Array{Float64, ndims(stategrid)}}(a_key => Array{Float64}(undef, size(stategrid)) for a_key in a_keys)
        end
    end

    # create sparsity
    J0c = sparsity_jac(stategrid, yend)

    # iterate on time
    if τs isa AbstractVector
        y_M = yend_M
        distances = zeros(length(τs))
        for iτ in length(τs):(-1):1
            a_keys !== nothing && _setindex!(as[iτ], localize(apm, τs[iτ]), stategrid, Tsolution, y_M, bc_M)
            _setindex!(ys[iτ], y_M)
            if iτ > 1
                y_M, distances[iτ] = implicit_timestep((ydot, y) -> hjb!(localize(apm, τs[iτ]), stategrid, Tsolution, ydot, y, bc_M, ysize), vec(y_M), τs[iτ] - τs[iτ-1]; is_algebraic = vec(is_algebraic_M), verbose = false, J0c = J0c, kwargs...)
                y_M = reshape(y_M, ysize...)
            end
        end
        return ys, as, distances
    else
        y_M, distance = finiteschemesolve((ydot, y) -> hjb!(apm, stategrid, Tsolution, ydot, y, bc_M, ysize), vec(yend_M); is_algebraic = vec(is_algebraic_M),  J0c = J0c, kwargs... )
        y_M = reshape(y_M, ysize...)
        _setindex!(y, y_M)
        if a_keys !== nothing
            _setindex!(a, apm, stategrid, Tsolution, y_M, bc_M)
            a = merge(y, a)
        end
        return y, a, distance
    end
end




_Array(yend) = cat(collect.(values(yend))...; dims = ndims(first(values(yend))) + 1)

function _Array_bc(bc, yend_M, yend, grid)
    bc_M = zero(yend_M)
    k = 0
    for yname in keys(yend)
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
_Array_bc(::Nothing, yend_M, yend, grid) = zero(yend_M)


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

function sparsity_jac(stategrid::StateGrid, yend::OrderedDict)
    s = size(stategrid)
    l = prod(s)
    t = (ndims(stategrid), length(yend) > 1)
    if t == (1, 0)
        J = Tridiagonal(ones(l - 1), ones(l), ones(l -1))
        return J, matrix_colors(J)
    elseif t == (2, 0)
        J = BandedBlockBandedMatrix(Ones(l, l), fill(s[1], s[2]), fill(s[1], s[2]), (1, 1), (1, 1))
        return sparse(J), matrix_colors(J)
    elseif t == (1, 1)
        J = BandedBlockBandedMatrix(Ones(l * length(yend), l * length(yend)), fill(l, length(yend)) ,fill(l, length(yend)), (length(yend) - 1, length(yend) - 1), (1, 1))
        return sparse(J), matrix_colors(J)
    elseif t == (2, 1)
        J = BandedBlockBandedMatrix(Ones(l * length(yend), l * length(yend)), repeat(fill(s[1], s[2]), outer = length(yend)), repeat(fill(s[1], s[2]), outer = length(yend)), (s[2] * length(yend) - 1, s[2] * length(yend) - 1), (1, 1))
        return sparse(J), matrix_colors(J)
    else
        return nothing, nothing
    end
end

function _setindex!(y::OrderedDict, y_M::AbstractArray)
    N = ndims(y_M) - 1
    i = 0
    for v in values(y)
        i += 1
        v[(Colon() for _ in 1:N)...] = y_M[(Colon() for _ in 1:N)..., i]
    end
end

function _setindex!(a::OrderedDict, apm, stategrid::StateGrid, Tsolution, y_M::AbstractArray, bc_M::AbstractArray)
    for i in eachindex(stategrid)
         state = stategrid[i]
         solution = derive(Tsolution, stategrid, y_M, i, bc_M)
         # upwind
         solution = derive(Tsolution, stategrid, y_M, i, bc_M, apm(state, solution)[2])
         outi = apm(state, solution)[3]
         for (k, v) in zip(values(a), values(outi))
            k[i] = v
         end
     end
 end

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


