"""
    pdesolve(f, grid, yend [, τs]; is_algebraic = OrderedDict(k => false for k in keys(yend)), bc = nothing)
    
    solves a system of pdes with an arbitrary number of functions and up to state variable. Denote  F the number of functions to solve for and S the number of state variables.

    ### Arguments
    *  f: is a function (state, y) -> out
    where: 
       state is a tuple of real numbers of  length S (for the state values)
       y is a tuple of real numbers of length F * 3 with the functions as well as their first and second derivatives at the state
       out is a named tuple of length F with the value of the time derivatives of the functions at the state 
    * grid is an NamedTuple that, for each state, associates an AbstractVector (for the grid)
    * yend is an NamedTuple that gives the terminal value of each function in the system of PDEs
    * τs (optional) is a time grid on which to solve the pde on. In this case, yend corresponds to the solution at time τs[end]
"""
function pdesolve(apm, @nospecialize(grid), @nospecialize(yend), τs::Union{Nothing, AbstractVector} = nothing; is_algebraic = OrderedDict(k => false for k in keys(yend)), bc = nothing, verbose = true, kwargs...)
    stategrid = StateGrid(NamedTuple(grid))
    S = size(stategrid)
    all(size(v) == S for v in values(yend)) || throw(ArgumentError("The length of initial guess (e.g. terminal value) does not equal the length of the state space"))
    all(keys(is_algebraic) .== keys(yend)) || throw(ArgumentError("the terminal guess yend and the is_algebric keyword argument must have the same names"))
    is_algebraic = OrderedDict(first(p) => fill(last(p), S) for p in pairs(is_algebraic))
    if τs isa AbstractVector
        issorted(τs) || throw(ArgumentError("The set of times must be increasing."))
        ys = [OrderedDict(first(p) => collect(last(p)) for p in pairs(yend)) for t in τs]
    else
        y = OrderedDict(first(p) => collect(last(p)) for p in pairs(yend))
    end
    # convert to Matrix
    yend_M = catlast(values(yend))
    is_algebraic_M = catlast(values(is_algebraic))
    bc_M = _Array_bc(bc, yend, grid)

    # create sparsity
    J0c = sparsity_jac(stategrid, yend)
    Tsolution = Type{tuple(keys(yend)...)}

    # iterate on time
    if τs isa AbstractVector
        apm_onestep = hasmethod(apm, Tuple{NamedTuple, NamedTuple, Number}) ? (state, grid) -> apm(state, grid, τs[end]) : apm
        a = get_a(apm_onestep, stategrid, Tsolution, yend_M, bc_M)
        as = (a === nothing) ? nothing : [deepcopy(a) for τ in τs]
        y_M = yend_M
        residual_norms = zeros(length(τs))
        if verbose
            @printf "    Time Residual\n"
            @printf "-------- --------\n"
        end
        for iτ in length(τs):(-1):1
            _setindex!(ys[iτ], y_M)
            apm_onestep = hasmethod(apm, Tuple{NamedTuple, NamedTuple, Number}) ? (state, grid) -> apm(state, grid, τs[iτ]) : apm
            if a !== nothing
                _setindex!(as[iτ], apm_onestep, stategrid, Tsolution, y_M, bc_M)
                as[iτ] = merge(ys[iτ], as[iτ])
            end
            if iτ > 1
                y_M, residual_norms[iτ] = implicit_timestep((ydot, y) -> hjb!(apm_onestep, stategrid, Tsolution, ydot, y, bc_M, size(yend_M)), vec(y_M), τs[iτ] - τs[iτ-1]; is_algebraic = vec(is_algebraic_M), verbose = false, J0c = J0c, kwargs...)
                if verbose
                    @printf "%8g   %8.4e\n" τs[iτ-1] residual_norms[iτ]
                end
                y_M = reshape(y_M, size(yend_M)...)
            end
        end
        return EconPDEResult(ys, residual_norms, as)
    else
        a = get_a(apm, stategrid, Tsolution, yend_M, bc_M)
        y_M, residual_norm = finiteschemesolve((ydot, y) -> hjb!(apm, stategrid, Tsolution, ydot, y, bc_M, size(yend_M)), vec(yend_M); is_algebraic = vec(is_algebraic_M),  J0c = J0c, verbose = verbose, kwargs... )
        y_M = reshape(y_M, size(yend_M)...)
        _setindex!(y, y_M)
        if a !== nothing
            _setindex!(a, apm, stategrid, Tsolution, y_M, bc_M)
            a = merge(y, a)
        end
        return EconPDEResult(y, residual_norm, a)
    end
end




catlast(iter) = cat(iter...; dims = ndims(first(iter)) + 1)

_Array_bc(::Nothing, yend, grid) = zeros(size(first(values(yend)))..., length(yend))
function _Array_bc(bc, yend, grid)
    keys_grid = collect(keys(grid))
    bc_M = _Array_bc(nothing, yend, grid)
    for (k, yname) in enumerate(keys(yend))
        if length(keys_grid) == 1
            bc_M[1, k], bc_M[end, k] = bc[Symbol(yname, keys_grid[1])]
        elseif length(keys_grid) == 2
            bc_M[1, :,  k],  bc_M[end, :, k] = bc[Symbol(yname, keys_grid[1])]
            bc_M[:, 1,  k], bc_M[:, end,  k] = bc[Symbol(yname, keys_grid[2])]
        end
    end
    return bc_M
end


function get_a(apm, stategrid::StateGrid, Tsolution, y_M::AbstractArray, bc_M::AbstractArray)
    i0 = first(eachindex(stategrid))
    derivatives = differentiate(Tsolution, stategrid, y_M, i0, bc_M)
    result = apm(stategrid[i0], derivatives)
    if length(result) == 1
        return nothing
    else
        return OrderedDict(a_key => Array{Float64}(undef, size(stategrid)) for a_key in keys(result[2]))
    end
end

function sparsity_jac(stategrid::StateGrid, @nospecialize(yend))
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

function _setindex!(@nospecialize(y), y_M::AbstractArray)
    for (i, v) in enumerate(values(y))
        v[:] = selectdim(y_M, ndims(y_M), i)
    end
end

function _setindex!(@nospecialize(a), apm, stategrid::StateGrid, Tsolution, y_M::AbstractArray, bc_M::AbstractArray)
    for i in eachindex(stategrid)
         solution = differentiate(Tsolution, stategrid, y_M, i, bc_M)
         outi = apm(stategrid[i], solution)[2]
         for (k, v) in zip(values(a), values(outi))
            k[i] = v
         end
     end
 end


 # create hjb! that accepts and returns AbstractVector rather than AbstractArrays
 function hjb!(apm, stategrid::StateGrid, Tsolution, ydot::AbstractVector, y::AbstractVector, bc_M::AbstractArray, ysize::NTuple)
     y_M = reshape(y, ysize...)
     ydot_M = reshape(ydot, ysize...)
     vec(hjb!(apm, stategrid, Tsolution, ydot_M, y_M, bc_M))
 end

 function hjb!(apm, stategrid::StateGrid, Tsolution, ydot_M::AbstractArray, y_M::AbstractArray, bc_M::AbstractArray)
     Tt = [Symbol(v, :t) for v in Tsolution.parameters[1]]
     for i in eachindex(stategrid)
         solution = differentiate(Tsolution, stategrid, y_M, i, bc_M)
         outi = apm(stategrid[i], solution)
         if isa(outi[1], Number)
            _setindex!(ydot_M, Tsolution, outi, i)
        else
            _setindex!(ydot_M, Tsolution, outi[1], i)
        end
     end
     return ydot_M
 end

 @generated function _setindex!(ydot_M::AbstractArray, ::Type{Tsolution}, outi::NamedTuple, i::CartesianIndex) where {Tsolution}
     N = length(Tsolution.parameters[1])
     quote
          $(Expr(:meta, :inline))
          $(Expr(:block, [Expr(:call, :setindex!, :ydot_M, Expr(:call, :getproperty, :outi, Meta.quot(Symbol(Tsolution.parameters[1][k], :t))), :i, k) for k in 1:N]...))
     end
 end


