δ = 1/30
n = 1000
σ = 0.1 * ones(n)
μ = 0.0027 * ones(n)
grid = collect(linspace(-2, 8, n))
ψ = zeros(n)
ψ[round(Int, n/4)] = 1