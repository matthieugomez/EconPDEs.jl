using Distributions

struct BansalYaronModel
    # consumption process parameters
    μbar::Float64 
    νc::Float64
    κμ::Float64 
    κσ::Float64 
    νμ::Float64 
    νσ::Float64 

    # utility parameters
    ρ::Float64  
    γ::Float64 
    ψ::Float64
end

function BansalYaronModel(;μbar = 0.018, νc = 0.027, κμ = 0.252, κσ = 0.156, νμ = 0.0143, νσ = 0.131, ρ = 0.024, γ = 7.5, ψ = 1.5)
    BansalYaronModel(μbar, νc, κμ, κσ, νμ, νσ, ρ, γ, ψ)
end
# Bansal, Kiku, Yaron (2009) calibration
#m = BansalYaronModel(μbar = 0.018, νc = 0.025, κμ = 0.3, κσ = 0.012, νμ = 0.0114, νσ = 0.189, ρ = 0.0132, γ = 7.5, ψ = 1.5)

function initialize_state(m::BansalYaronModel; μn = 30, σn = 30)
    μbar = m.μbar ; νc = m.νc ; κμ = m.κμ ; κσ = m.κσ ; νμ = m.νμ ; νσ = m.νσ ; ρ = m.ρ ; γ = m.γ ; ψ = m.ψ

    σ = sqrt(νμ^2 / (2 * κμ))
    μmin = quantile(Normal(μbar, σ), 0.001)
    μmax = quantile(Normal(μbar, σ), 0.999)
    μs = collect(linspace(μmin, μmax, μn))

    σ = sqrt(νσ^2 / (2 * κσ))
    σmin = max(0.01, quantile(Normal(1.0, σ), 0.001))
    σmax = quantile(Normal(1.0, σ), 0.999)
    σs = collect(linspace(σmin, σmax, σn))

    OrderedDict(:μ => μs, :σ => σs)
end

function initialize_y(m::BansalYaronModel, state)
    OrderedDict(:p => fill(1.0, length(state[:μ]), length(state[:σ])))
end

function (m::BansalYaronModel)(state, y)
    μbar = m.μbar ; νc = m.νc ; κμ = m.κμ ; κσ = m.κσ ; νμ = m.νμ ; νσ = m.νσ ; ρ = m.ρ ; γ = m.γ ; ψ = m.ψ
    μ, σ = state.μ, state.σ
    p, pμ, pσ, pμμ, pμσ, pσσ = y.p, y.pμ, y.pσ, y.pμμ, y.pμσ, y.pσσ

    # drift and volatility of c, μ, σ, p
    μc = μ
    σc_Zc = νc * sqrt(σ)
    μμ = κμ * (μbar - μ)
    σμ_Zμ = νμ * sqrt(σ)
    μσ = κσ * (1 - σ)
    σσ_Zσ = νσ 
    σp_Zμ = pμ / p * σμ_Zμ
    σp_Zσ = pσ / p * σσ_Zσ
    σp2 = σp_Zμ^2 + σp_Zσ^2
    μp = pμ / p * μμ + pσ / p * μσ + 0.5 * pμμ / p * σμ_Zμ^2 + 0.5 * pσσ / p * σσ_Zσ^2

    # Market price of risk κ
    κ_Zc = γ * σc_Zc
    κ_Zμ = - (1 - γ * ψ) / (ψ - 1) * σp_Zμ
    κ_Zσ = - (1 - γ * ψ) / (ψ - 1) * σp_Zσ
    κσC = κ_Zc * σc_Zc
    κσp = κ_Zμ * σp_Zμ + κ_Zσ * σp_Zσ
    κ2 = κ_Zc^2 + κ_Zμ^2 + κ_Zσ^2

    # Risk free rate r
    r = ρ + 1 / ψ * (μc - (1 + ψ)/ (2 * γ) * κ2 - (1 - ψ * γ) / (γ * (ψ - 1)) * κσp + (1 - γ * ψ) / (2 * γ * (ψ - 1)) * σp2)

    # PDE
    pt = p * (1 / p + μc + μp - r - κσC - κσp)
    # dp = p * (1 / p - ρ + (1 - 1 / ψ) * (μc - 0.5 * γ * σc_Zc^2) + μp + 0.5 * (1 / ψ - γ) / (1 - 1 / ψ) * σp2)

    return pt, (μμ, μσ), tuple(:p => p, :μμ => μμ, :σμ_Zμ => σμ_Zμ, :σμ_Zσ => 0.0, :μσ => μσ, :σσ_Zμ => 0.0, :σσ_Zσ => σσ_Zσ, :μ => μ, :σ => σ, :σμ2 => σμ_Zμ^2, :σσ2 => σσ_Zσ^2, :σμσσ => 0.0)
end

## Long Run Risk Models
### Bansal Yaron (2004)
#m = BansalYaronModel()
#state = initialize_state(m)
#y0 = initialize_y(m, state)
#result, distance = pde_solve(m, state, y0)

### Bansal, Kiku, Yaron (2009) calibration
#m = BansalYaronModel(μbar = 0.018, νc = 0.025, κμ = 0.3, κσ = 0.012, νμ = 0.0114, νσ = 0.189, ρ = 0.0132, γ = 7.5, ψ = 1.5#)
#state = initialize_state(m)
#y0 = initialize_y(m, state)
#result, distance = pde_solve(m, state, y0)