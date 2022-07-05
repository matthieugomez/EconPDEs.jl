using EconPDEs, Distributions

# Wachter (2013) Can time‐varying risk of rare disasters explain aggregate stock market volatility?
Base.@kwdef struct WachterModel{T<: Distribution}
    # consumption process parameters
    μ::Float64 = 0.025
    σ::Float64 = 0.02

    # idiosyncratic
    λbar::Float64 = 0.0355
    κλ::Float64 = 0.08
    νλ::Float64 = 0.067
    ZDistribution::T = Normal(-0.4, 0.25) # from Ian Martin higher order cumulants paper
    
    # utility parameters
    ρ::Float64 = 0.012
    γ::Float64 = 3.0
    ψ::Float64 = 1.1
    ϕ::Float64 = 2.6
end

function initialize_stategrid(m::WachterModel; λn = 30)
  OrderedDict(:λ => range(0.0, 0.1, length = λn))
end

function initialize_y(m::WachterModel, stategrid)
    λn = length(stategrid[:λ])
    OrderedDict(:p => ones(λn))
end


function (m::WachterModel)(state::NamedTuple, y::NamedTuple)
    (; μ, σ, λbar, κλ, νλ, ZDistribution, ρ, γ, ψ, ϕ) = m
    (; λ) = state
    (; p, pλ_up, pλ_down, pλλ) = y

    # Drift and volatility of λ, p
    μλ = κλ * (λbar - λ)
    σλ = νλ * sqrt(λ)
    # upwinding
    pλ = (μλ >= 0) ? pλ_up : pλ_down
    μp = pλ / p * μλ + 0.5 * pλλ / p * σλ^2 
    σp_Zλ = pλ / p * σλ

    # Market price of risk 
    κ_Zc = γ * σ 
    κ_Zλ = (γ * ψ - 1) / (ψ - 1) * σp_Zλ
    η =  λ * (mgf(ZDistribution, 1) - 1 + mgf(ZDistribution, -γ) - 1 - (mgf(ZDistribution, 1 - γ) - 1))

    # Interest rate r
    r = ρ + μ / ψ - (1 + 1 / ψ) / 2 * γ * σ^2 - λ * (mgf(ZDistribution, -γ) - 1) + (1 / ψ - γ) / (1 - γ) * λ * (mgf(ZDistribution, 1 - γ) - 1) - (γ * ψ - 1) / (2 * (ψ - 1)) * σp_Zλ^2

    # Market Pricing
    pt = - p * (1 / p  + μ + μp + λ * (mgf(ZDistribution, 1) - 1) - r - κ_Zc * σ - κ_Zλ * σp_Zλ - η)
    
    return (; pt)
end

m = WachterModel()
stategrid = initialize_stategrid(m)
yend = initialize_y(m, stategrid)
@assert result.residual_norm <= 1e-5