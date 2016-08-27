##############################################################################
##
## Reflecting Array (i.e. 0 is 1 and N+1 is N)
##
##############################################################################

if !isdefined(:ReflectingArray)
    type ReflectingArray{T, N}
        A::Array{T, N}
    end
end
Base.size(y::ReflectingArray, args...) = size(y.A, args...)
Base.eltype(y::ReflectingArray) = eltype(y.A)
Base.eachindex(y) = eachindex(y.A)
@generated function Base.getindex{N, T}(A::ReflectingArray{T, N}, args...)
    expr = tuple([_helper(args[i], i) for i in 1:N]...)
    Expr(:call, :getindex, :(A.A), expr...)
end
_helper(x::Type{Int64}, i) = :(min(max(args[$i], 1), size(A, $i)))
_helper(x::Type{UnitRange{Int64}}, i) = :(:)
_helper(x::Type{Colon}, i) = :(:)


##############################################################################
##
## State Space Grid
##
##############################################################################

type StateSpace{N}
    a::NTuple{N, Vector{Float64}}
    inva::NTuple{N, Float64}
    n::NTuple{N, Int}
end

function StateSpace(args...)
    inva = map(x -> 1 / (x[2] - x[1]), args)
    n = map(length, args)
    StateSpace(args, inva, n)
end
Base.eachindex(grid::StateSpace) = CartesianRange(grid.n)
@generated function Base.getindex{N}(grid::StateSpace{N}, args...)
    Expr(:call, :tuple, [:(getindex(grid.a[$i], args[$i])) for i in 1:N]...)
end

function convert(::Type{Tuple{N}{Vector{Float64}}}, x::StateSpace{N})
    A = Array{Float64, N}[zeros(grid.n) for i in 1:N]
    for ituple in eachindex(x)
        for j in 1:N
            A[j][ituple...] = x.a[j][ituple[j]]
        end
    end
    return tuple(A...)
end



##############################################################################
##
## HJB Solver
##
## solves y such that 0 = max \{ f(c, y), + E[dy] \}
##
## F!(y, ydot) is a function that updates in place ydot to max \{ f(c, y), + E[dy] \}
##
##############################################################################


function F_step!(F!, y, ydot, ypost, invΔ)
    for i in 1:length(ydot)
        ydot[i] = ydot[i] + (ypost[i] - y[i]) * invΔ
    end
    return ydot
end

function iterate_backward(X, ypost, invΔ; iterations = 100, method = :newton, maxdist = 1e-9, verbose = false) 
    y0 = deepcopy(ypost)
    out = nlsolve((y, ydot) -> F_step!(F!, y, ydot, ypost, invΔ), y0, method = method, autodiff = true, show_trace = verbose, ftol = maxdist, iterations = iterations)
    return out.zero, out.residual_norm
end

function solve_backward(F!, y0; iterations = 100, method = :newton, maxdist = 1e-9, verbose = true, autodiff = true, kwargs...)
    ypost = deepcopy(y0)
    y = deepcopy(y0)
    ydot = deepcopy(y0)
    distance = Inf
    olddistance = Inf
    oldolddistance = Inf
    invΔ = 10
    iter = 0
    while (iter <= iterations) & (invΔ <= 1e6)
        iter += 1
        y, distance = iterate_backward(F!, ypost, invΔ)
        distance = maxabs(F!(y, ydot))
        if show_trace
            @show iter, distance
        end
        if (distance >=  olddistance)
            invΔ = 10 * invΔ
            copy!(v, ypost)
        else
            invΔ = max(1e-3, invΔ / 10)
        end
        if (distance <= maxdist) 
            break
        else
            ypost, y = y, ypost
            olddistance, oldolddistance = distance, olddistance
        end
    end
    return y, distance
end


##############################################################################
##
## Model Solver
##
##############################################################################

function F_fd!(byp::::AbstractAssetPricingContinuousTimeModel, grid::StateSpace, y::Vector, ydot::Vector)
    byp, grid = X
    n = n_functions(byp)
    fy = ReflectingArray(reshape(y, grid.n..., n))
    fydot = reshape(ydot, grid.n..., n)
    for ituple in eachindex(grid)
        gridi = derive(byp, ituple, grid, fy)
        out, μ, others = pde(byp, gridi)
        gridi = derive(byp, ituple, grid, fy, μ)
        out, μ, others = pde(byp, gridi)
        for k in 1:n
            fydot[ituple, k] = out[k]
        end
    end
    return ydot
end



function solve(byp::AbstractAssetPricingContinuousTimeModel, grid::StateSpace)
    F!(y, ydot) = F_Fd!(byp, grid, y, ydot)
    solve_backward((y, ydot) -> F_Fd!(byp, grid, y, ydot), y0)
    y, distance
end