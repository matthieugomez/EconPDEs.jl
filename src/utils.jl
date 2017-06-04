
#========================================================================================

Stationary Distribution

========================================================================================#

# Let's prove that for A matrix with off diagonal non negative + sum of columns = 0, there exists g non negative vector such that Ag = 0
# Proof: A  has off-diagonal element non negative, therefore can be written as r(I - B) with B non negative matrix.  Sum of A columns are 0, implies sum of B columns are 1. We can apply math of markov chain: there exists g non negative vector such that Bg = I, i.e. Ag = 0. Note that -A is actually said to be a "singular M matrix"

# Let's prove that for for A matrix with off diagonal non negative + sum of columns = 0, there exists g non negative vector such that  (δI - A) g = δψ.  where ψ is psotivie elementwise
# Proof We know from paragraph above that zero is the highest eigenvalue for A thereofre δI - A is invertible and it is a Z matrix therefore it is a M matrix. Therefore we know that (δI - a)^{-1} has all elements positive  and therefore ψ has all element non negative implies g has all element non negative.
# Note that this also proof the case without death rate by using limit argument.

# Note that this set of conditions are the same as barles souganadas conditions for convergence of finite scheme. Scheme monoticity ensures off diagonal elements are psotivie or equal tto zero. Scheme consistency implies sum of A columbs are 0.



#now there are still two issues
#1. Does not satisfy walras law. Or mathematically does not satisfy IPP ∑ μ.g = ∑ a.Ag. In the case drift is positive, there is a remaning term μ_NdG(a_N) To satisfy it, do amax super super high (intuitively, x high enough so that cutting behavior at the top does not matter for aggregate as g(x)x -> 0)
#2. A g can be negative when updating forward. Use implicit scheme

# Case with 1 state variable
function clean!(density)
    scale!(density, 1/sum(density))
end

function computeA(grid, μ::Vector{T}, σ2::Vector{T}) where {T <: Number}
    x, invΔx, invΔxm, invΔxp = make_Δ(grid)
    n = length(x)
    A = zeros(n, n)
    for i in 1:n
        # ensures that sum of columns is always 0
        if μ[i] >= 0
            A[min(i + 1, size(A, 1)), i] += μ[i] * invΔxp[i]
            A[i, i] -= μ[i] * invΔxp[i]
        else
            A[i, i] += μ[i] * invΔxm[i] 
            A[max(i - 1, 1), i] -= μ[i] * invΔxm[i] 
        end
        A[max(i - 1, 1), i] += 0.5 * σ2[i] * invΔx[i] * invΔxm[i] 
        A[i, i] -= 0.5 * σ2[i] * 2 * invΔxm[i] * invΔxp[i]
        A[min(i + 1, size(A, 1)), i] += 0.5 * σ2[i] * invΔx[i] * invΔxp[i]
    end
    return A
end


function stationary_distribution(grid, μ::Vector{T}, σ2::Vector{T}) where {T <: Number}
    A = computeA(grid, μ, σ2)
    n = size(A, 2)
    for j in 1:n
        A[1, j] = 1.0
    end
    b = vcat(1.0, zeros(n - 1))
    density = A \ b
    @assert all(density .> 0)
    @assert sum(density) ≈ 1.0
    clean!(density) 
end

function stationary_distribution(grid, μ::Vector{T}, σ2::Vector{T}, δ, ψ) where {T <: Number}
    A = computeA(grid, μ, σ2)
    density = (δ .* eye(A) .- A) \ (δ .* ψ)
    @assert all(density .> 0)
    @assert sum(density) ≈ 1.0
    clean!(density)
end




# Case with 2 state variables
function stationary_distribution(grid, μ::Vector{Array{T, 2}}, σ2::Vector{Array{T, 2}}) where {T}
    x1, invΔx1, invΔx1m, invΔx1p = make_Δ(grid[1])
    x2, invΔx2, invΔx2m, invΔx2p = make_Δ(grid[2])
    for i2 in 1:n2
        for i1 in 1:n1
            if μ[1][i1, i2] >= 0
                i1h = min(i1 + 1, size(A, 1))
                i1l = i1
                invΔx1i = invΔx1p[i1]
            else
               i1h = i1
               i1l = max(i1 - 1, 1)
               invΔx1i = invΔx1m[i1]
            end
            A[i1h, i2, i1, i2] += μ[1][i1, i2] * invΔx1i
            A[i1l, i2, i1, i2] -= μ[1][i1, i2] * invΔx1i
            A[max(i1 - 1, 1), i2, i1, i2] += 0.5 * σ2[1][i1, i2] * invΔx1[i1] * invΔx1m[i1] 
            A[i1, i2, i1, i2] -= 0.5 * σ[1][i1, i2] * 2 * invΔx1m[i1] * invΔx1p[i1]
            A[min(i1 + 1, size(A, 1)),i2, i1, i2] += 0.5 * σ2[1][i1, i2] * invΔx1[i1] * invΔx1p[i1]
            if μ[2][i1, i2] >= 0
                i2h = min(i2 + 1, size(A, 2))
                i2l = i2
                invΔx2i = invΔx2p[i1]
            else
               i2h = i2
               i2l = max(i2 - 1, 1)
               invΔx2i = invΔx2m[i2]
            end
            A[i1, i2h, i1, i2] += μ[2][i1, i2] * invΔx2i
            A[i1, i2l, i1, i2] -= μ[2][i1, i2] * invΔx2i
            A[i1, max(i2 - 1, 1), i1, i2] += 0.5 * σ2[2][i1, i2] * invΔx2[i2] * invΔx2m[i2] 
            A[i1, i2, i1, i2] -= 0.5 * σ[2][i1, i2] * 2 * invΔx2m[i2] * invΔx2p[i2]
            A[i1,min(i2 + 1, size(A, 2)), i1, i2] += 0.5 * σ2[2][i1, i2] * invΔx2[i2] * invΔx2p[i2]

            #this way only adds bad sign in stuff already filled (at least when sigma2[3] ⫺ 0). works as long as this expression smaller than the drift
            A[i1h, i2h, i1, i2] += σ2[3][i1, i2] * invΔx1i * invΔx2i
            A[i1l, i2h, i1, i2] -= σ2[3][i1, i2] * invΔx1i * invΔx2i
            A[i1h, i2l, i1, i2] -= σ2[3][i1, i2] * invΔx1i * invΔx2i
            A[i1l, i2l, i1, i2] += σ2[3][i1, i2] * invΔx1i * invΔx2i
        end
    end
    A = reshape(A.A, (n1 * n2, n1 * n2))
    for j in 1:size(A, 2)
        A[1, j] = 1.0
    end
    b = vcat(1.0, zeros(size(A, 2) - 1))
    density = A \ b
    density = abs.(density) ./ sum(abs, density)  
    return reshape(density, (n1, n2))
end




function stationary_distribution(grid::OrderedDict, a::OrderedDict)
    k = collect(keys(grid))
    if length(k) == 1
        stationary_distribution(collect(values(grid))[1], a[Symbol(:μ, k[1])],  a[Symbol(:σ, k[1])].^2)
    elseif length(k) == 2
        stationary_distribution(collect(values(grid)), [a[Symbol(:μ, k[1])], a[Symbol(:μ, k[2])]],  [a[Symbol(:σ, k[1])].^2, a[Symbol(:σ, k[2])].^2, a[Symbol(:σ, k[1], k[2])]])
    end
end
#========================================================================================

Simulate

========================================================================================#

function simulate(grid, a, shocks::OrderedDict; dt = 1 / 12, x0 = nothing)
    grid = StateGrid(grid)
    T = size(shocks[first(keys(shocks))], 1)
    I = size(shocks[first(keys(shocks))], 2)
    if x0 == nothing
        i0 = rand(Categorical(vec(stationary_distribution(grid, a))), I)
        if N == 1
            x0 = OrderedDict(grid.name[1] => grid.x[1][i0])
        elseif N == 2
            i10 = mod(i0 - 1, length(grid.x[1])) + 1
            i20 = div(i0 - 1, length(grid.x[1])) + 1
            x0 = OrderedDict(grid.name[1] => grid.x[1][i10], grid.name[2] => grid.x[2][i20])
        end
    end
    # interpolate all functions
    ai = OrderedDict([Pair(k => interpolate(grid.x, a[k], Gridded(Linear()))) for k in keys(a)])
    aT = OrderedDict([Pair(k => zeros(T, I)) for k in keys(a)])
    aT[:id] = zeros(T, I)
    aT[:t] = zeros(T, I)
    for k in keys(shocks)
        aT[k] = zeros(T, I)
    end
    sqrtdt = sqrt(dt)
    for id in 1:I
        xt = tuple([x0[grid.name[i]][id] for i in 1:N]...)
        for t in 1:T
            for k in keys(a)
                aT[k][t, id] = ai[k][xt...]
            end
            aT[:id][t, id] = id
            aT[:t][t, id] = t
            for k in keys(shocks)
                aT[k][t, id] = shocks[k][t, id]
            end
            xt = tuple([xt[i] + _update_state(xt, grid.name[i], shocks, ai, t, id, dt, sqrtdt) for i in 1:N]...)
        end
    end
    return aT
end

function _update_state(xt, name, shocks, ai, t, id, dt, sqrtdt)
    out = ai[Symbol(:μ, name)][xt...] * dt
    if length(keys(shocks)) == 1
        out += ai[Symbol(:σ, name)][xt...] * shocks[first(keys(shocks))][t, id] * sqrtdt
    else
        for k in keys(shocks)
            out += ai[Symbol(:σ, name, :_, k)][xt...] * shocks[k][t, id] * sqrtdt
        end
    end
    return out
end
