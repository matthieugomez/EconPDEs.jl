#========================================================================================

Jacobian sparsity pattern and column coloring.

The finite-difference residual at a grid point only depends on the unknowns on the ±1
box around that point, so the Jacobian pattern is known in closed form
(`local_stencil_jacobian`), as is an optimal distance-2 coloring of its columns
(`stencil_colors`). The coloring lets FiniteDiff compute the whole sparse Jacobian in
one function evaluation per color.

========================================================================================#

function sparse_jacobian(stategrid::StateGrid, @nospecialize(guess))
    s = size(stategrid)
    if 1 <= ndims(stategrid) <= 3
        return local_stencil_jacobian(s, length(guess))
    else
        return nothing
    end
end

function local_stencil_jacobian(s::NTuple{N, Int}, F::Int) where {N}
    l = prod(s)
    nnz_max = l * F^2 * 3^N
    rows = Vector{Int}(undef, nnz_max)
    cols = Vector{Int}(undef, nnz_max)
    k = 0
    state_indices = CartesianIndices(s)
    for fout in 1:F, ci in state_indices
        row = LinearIndices(s)[ci] + (fout - 1) * l
        stencil_ranges = ntuple(d -> max(ci[d] - 1, 1):min(ci[d] + 1, s[d]), N)
        for fin in 1:F, neighbor in CartesianIndices(stencil_ranges)
            k += 1
            rows[k] = row
            cols[k] = LinearIndices(s)[neighbor] + (fin - 1) * l
        end
    end
    return sparse(rows[1:k], cols[1:k], ones(k), l * F, l * F)
end

# Closed-form distance-2 coloring of the pattern built by `local_stencil_jacobian` —
# orders of magnitude faster than greedy coloring of the assembled matrix, with the
# same (provably optimal) number of colors. Keep the two functions in sync: the
# argument below relies on the stencil being the full ±1 box coupling all unknowns.
#
# Two columns (grid point c, unknown f) and (c′, f′) share a nonzero row iff c and c′
# are within Chebyshev distance 2: every row couples all unknowns on the ±1 box around
# its grid point, so f plays no role, and a witness point within distance 1 of both c
# and c′ always exists on the grid (boundary clipping never removes it). Coloring by
# coordinate mod 3 in each dimension, crossed with the unknown index, is therefore
# valid: two same-colored columns differ by a nonzero multiple of 3 in some dimension.
# It is also optimal: any 3×…×3 sub-box crossed with the F unknowns is a clique of
# exactly F·∏ min(3, s_d) columns. Mixed-radix encoding with radix min(3, s_d) keeps
# the colors contiguous 1:ncolors, so FiniteDiff wastes no function evaluation on an
# empty color class when some dimension has fewer than 3 points.
function stencil_colors(s::NTuple{N, Int}, F::Int) where {N}
    l = prod(s)
    radices = map(sd -> min(3, sd), s)
    colors = Vector{Int}(undef, l * F)
    for (i, ci) in enumerate(CartesianIndices(s))
        base = 0
        stride = 1
        for d in 1:N
            base += ((ci[d] - 1) % 3) * stride
            stride *= radices[d]
        end
        colors[i] = base + 1
    end
    ncolors_grid = prod(radices)
    for f in 2:F, j in 1:l
        colors[(f - 1) * l + j] = colors[j] + (f - 1) * ncolors_grid
    end
    return colors
end

# Local replacement for the old SparseDiffTools/ArrayInterface `matrix_colors`
# helper: a greedy column coloring of the sparse Jacobian pattern.
function matrix_colors(A::SparseMatrixCSC)
    rows = rowvals(A)
    row_columns = [Int[] for _ in 1:size(A, 1)]
    for j in 1:size(A, 2)
        for ptr in nzrange(A, j)
            push!(row_columns[rows[ptr]], j)
        end
    end

    colors = zeros(Int, size(A, 2))
    forbidden = zeros(Int, size(A, 2))
    maxcolor = 0
    for j in 1:size(A, 2)
        for ptr in nzrange(A, j)
            for k in row_columns[rows[ptr]]
                color = colors[k]
                if color > 0
                    forbidden[color] = j
                end
            end
        end
        color = 1
        while color <= maxcolor && forbidden[color] == j
            color += 1
        end
        colors[j] = color
        maxcolor = max(maxcolor, color)
    end
    return colors
end
