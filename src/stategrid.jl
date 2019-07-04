#========================================================================================

Type State Grid

========================================================================================#
struct StateGrid{N, V}
    x::NTuple{N, Vector{Float64}}
    invΔx::NTuple{N, Vector{Float64}}
    invΔxm::NTuple{N, Vector{Float64}}
    invΔxp::NTuple{N, Vector{Float64}}
end

function make_Δ(x::AbstractVector)
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
Base.ndims(grid::StateGrid) = length(grid.x)
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


