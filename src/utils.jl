"""
StrateGrid(x::NamedTuple)

Let x = (k1 = v1, k2 = v2, ...., kN = vN) be a NamedTuple of AbstractVectors for state space, 
StateGrid(x) returns an AbstractArray M such that 
M[i1, i2..., iN] = (k1 = v1[i1], k2 = v2[i2], ...., kN = vN[iN]
"""
struct StateGrid{T, N, C <: NamedTuple} <: AbstractArray{T, N}
    x::C
end
function StateGrid(x::NamedTuple{Names, <: NTuple{N, <: AbstractVector{T}}}) where {Names, N, T}
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

#========================================================================================

Derive

========================================================================================#

# Helpers for building index expressions at compile time.
# Called inside the @generated function body to construct Expr nodes.

# Build y[i1, ..., iN, k] with optional substitutions at specific dimensions.
function _y_ref(N::Int, k::Int, subs::Pair{Int}...)
    idx = Any[Symbol("i", d) for d in 1:N]
    for (d, v) in subs
        idx[d] = v
    end
    push!(idx, k)
    Expr(:ref, :y, idx...)
end

# Build bc[i1, ..., (1|end), ..., iN, k] with dimension d fixed to boundary.
function _bc_ref(N::Int, k::Int, d::Int, boundary)
    idx = Any[Symbol("i", dim) for dim in 1:N]
    idx[d] = boundary  # 1 or :end
    push!(idx, k)
    Expr(:ref, :bc, idx...)
end

@generated function differentiate(::Type{Tsolution}, grid::StateGrid{T1, Ndim, <: NamedTuple{Names}}, y::AbstractArray{T}, icar, bc) where {Tsolution, T1, Ndim, Names, T}
    solnames = Tsolution.parameters[1]
    statenames = Names

    # Preamble: extract indices and compute grid spacings
    preamble = Expr[]
    for d in 1:Ndim
        id = Symbol("i", d)
        gd = Symbol("grid", d)
        Δm = Symbol("Δx", d, "m")
        Δp = Symbol("Δx", d, "p")
        Δ  = Symbol("Δx", d)
        push!(preamble, :($id = icar[$d]))
        push!(preamble, :($gd = grid.x[$d]))
        push!(preamble, :(@inbounds $Δm = $gd[max($id, 2)] - $gd[max($id - 1, 1)]))
        push!(preamble, :(@inbounds $Δp = $gd[min($id + 1, size(y, $d))] - $gd[min($id, size(y, $d) - 1)]))
        push!(preamble, :($Δ = ($Δm + $Δp) / 2))
    end

    # Build derivative expressions for each solution variable
    expr = Expr[]
    for k in 1:length(solnames)
        solname = solnames[k]

        # Value
        push!(expr, Expr(:(=), solname, _y_ref(Ndim, k)))

        # First derivatives (upwind and downwind) per dimension
        for d in 1:Ndim
            sn = statenames[d]
            id = Symbol("i", d)
            Δp = Symbol("Δx", d, "p")
            Δm = Symbol("Δx", d, "m")
            y_base = _y_ref(Ndim, k)
            y_fwd  = _y_ref(Ndim, k, d => :($id + 1))
            y_bwd  = _y_ref(Ndim, k, d => :($id - 1))
            bc_hi  = _bc_ref(Ndim, k, d, :end)
            bc_lo  = _bc_ref(Ndim, k, d, 1)
            push!(expr, Expr(:(=), Symbol(solname, sn, :_up),
                :(($id < size(y, $d)) ? ($y_fwd - $y_base) / $Δp : convert($T, $bc_hi))))
            push!(expr, Expr(:(=), Symbol(solname, sn, :_down),
                :(($id > 1) ? ($y_base - $y_bwd) / $Δm : convert($T, $bc_lo))))
        end

        # Second derivatives per dimension (with boundary conditions)
        for d in 1:Ndim
            sn = statenames[d]
            id = Symbol("i", d)
            Δp = Symbol("Δx", d, "p")
            Δm = Symbol("Δx", d, "m")
            Δ  = Symbol("Δx", d)
            y_fwd  = _y_ref(Ndim, k, d => :($id + 1))
            y_base = _y_ref(Ndim, k)
            y_bwd  = _y_ref(Ndim, k, d => :($id - 1))
            # left boundary (id == 1)
            y_at_2  = _y_ref(Ndim, k, d => 2)
            y_at_1  = _y_ref(Ndim, k, d => 1)
            bc_lo   = _bc_ref(Ndim, k, d, 1)
            # right boundary (id == end)
            y_at_end   = _y_ref(Ndim, k, d => :(size(y, $d)))
            y_at_endm1 = _y_ref(Ndim, k, d => :(size(y, $d) - 1))
            bc_hi      = _bc_ref(Ndim, k, d, :end)
            interior = :($y_fwd / ($Δp * $Δ) + $y_bwd / ($Δm * $Δ) - 2 * $y_base / ($Δp * $Δm))
            left_bc  = :($y_at_2 / ($Δp * $Δ) + ($y_at_1 - $bc_lo * $Δm) / ($Δm * $Δ) - 2 * $y_at_1 / ($Δp * $Δm))
            right_bc = :(($y_at_end + $bc_hi * $Δp) / ($Δp * $Δ) + $y_at_endm1 / ($Δm * $Δ) - 2 * $y_at_end / ($Δp * $Δm))
            push!(expr, Expr(:(=), Symbol(solname, sn, sn),
                :((1 < $id < size(y, $d)) ? $interior : (($id == 1) ? $left_bc : $right_bc))))
        end

        # Cross derivatives for each pair of dimensions
        for d1 in 1:Ndim, d2 in (d1+1):Ndim
            sn1 = statenames[d1]
            sn2 = statenames[d2]
            id1 = Symbol("i", d1)
            id2 = Symbol("i", d2)
            Δ1  = Symbol("Δx", d1)
            Δ2  = Symbol("Δx", d2)
            y_pp = _y_ref(Ndim, k, d1 => :(min($id1 + 1, size(y, $d1))), d2 => :(min($id2 + 1, size(y, $d2))))
            y_pm = _y_ref(Ndim, k, d1 => :(min($id1 + 1, size(y, $d1))), d2 => :(max($id2 - 1, 1)))
            y_mp = _y_ref(Ndim, k, d1 => :(max($id1 - 1, 1)),            d2 => :(min($id2 + 1, size(y, $d2))))
            y_mm = _y_ref(Ndim, k, d1 => :(max($id1 - 1, 1)),            d2 => :(max($id2 - 1, 1)))
            push!(expr, Expr(:(=), Symbol(solname, sn1, sn2),
                :(($y_pp - $y_pm - $y_mp + $y_mm) / (4 * $Δ1 * $Δ2))))
        end
    end

    quote
        $(Expr(:meta, :inline))
        $(preamble...)
        @inbounds $(Expr(:tuple, expr...))
    end
end

#========================================================================================

matrix_colors

========================================================================================#

# from ArraysInterface (could just import it but might be big import and changing all the time)
matrix_colors(A::Tridiagonal) = _cycle(1:3, size(A, 2))
function matrix_colors(A::BandedBlockBandedMatrix)
    l, u = blockbandwidths(A)
    lambda, mu = subblockbandwidths(A)
    blockwidth = l + u + 1
    subblockwidth = lambda + mu + 1
    nblock = blocksize(A, 2)
    cols = blocklengths(axes(A, 2))
    blockcolors = _cycle(1:blockwidth, nblock)
    # the reserved number of colors of a block is the min of subblockwidth and the largest length of columns of blocks with the same block color
    ncolors = [
        min(subblockwidth, maximum(cols[i:blockwidth:nblock]))
        for i = 1:min(blockwidth, nblock)
    ]
    endinds = cumsum(ncolors)
    startinds = [endinds[i] - ncolors[i] + 1 for i = 1:min(blockwidth, nblock)]
    colors = [
        _cycle(startinds[blockcolors[i]]:endinds[blockcolors[i]], cols[i])
        for i = 1:nblock
    ]
    return reduce(vcat, colors)
end
_cycle(repetend, len) = repeat(repetend, div(len, length(repetend)) + 1)[1:len]
