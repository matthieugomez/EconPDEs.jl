#========================================================================================

Type State Grid

========================================================================================#
struct StateGrid{N, V}
    x::NTuple{N, Vector{Float64}}
end

function StateGrid(x)
    StateGrid{length(x), tuple(keys(x)...)}(tuple(values(x)...))
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


