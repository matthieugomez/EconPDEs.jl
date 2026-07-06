"""
StateGrid(x::NamedTuple)

Let x = (k1 = v1, k2 = v2, ...., kN = vN) be a NamedTuple of AbstractVectors for state space,
StateGrid(x) returns an AbstractArray M such that
M[i1, i2..., iN] = (k1 = v1[i1], k2 = v2[i2], ...., kN = vN[iN]
"""
struct StateGrid{T, N, C <: NamedTuple} <: AbstractArray{T, N}
    x::C
end

function StateGrid(x::NamedTuple{Names, V}) where {Names, V <: Tuple}
    N = length(Names)
    N > 0 || throw(ArgumentError("state grid must contain at least one state variable"))
    all(v -> v isa AbstractVector, values(x)) || throw(ArgumentError("state grid entries must be vectors"))
    all(v -> length(v) >= 2, values(x)) || throw(ArgumentError("state grid entries must have at least 2 points"))
    T = promote_type(map(eltype, values(x))...)
    StateGrid{T, N, typeof(x)}(x)
end
function Base.eltype(stategrid::StateGrid{T, N, <: NamedTuple{Names, V}}) where {T, N, Names, V}
    NamedTuple{Names, NTuple{N, T}}
end
Base.ndims(stategrid::StateGrid{T, N}) where {T, N} = N
Base.size(stategrid::StateGrid{T, N}) where {T, N} = ntuple(i -> length(stategrid.x[i]), N)
Base.eachindex(stategrid::StateGrid) = CartesianIndices(size(stategrid))
function Base.getindex(stategrid::StateGrid, args::CartesianIndex)
    eltype(stategrid)(ntuple(i -> stategrid.x[i][args[i]], ndims(stategrid)))
end
Base.getindex(stategrid::StateGrid, x::Symbol) = stategrid.x[x]
