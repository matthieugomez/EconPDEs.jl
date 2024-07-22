# Haddad WP "Concentrated Ownership and Equilibrium Asset Prices"

using EconPDEs, Distributions


Base.@kwdef struct HaddadModel
  # consumption process parameters
  μbar::Float64 = 0.018
  vbar::Float64 = 0.00052
  κμ::Float64 = 0.3
  νμ::Float64 = 0.456
  κv::Float64 = 0.012
  νv::Float64 = 0.00472

  # active capital
  αbar::Float64 = 1.2
  λ::Float64 = 0.018

  # utility parameters
  ρ::Float64 = 0.0132
  γ::Float64 = 10.0
  ψ::Float64 = 1.5
end

function (m::HaddadModel)(state::NamedTuple, y::NamedTuple)
  (; μbar, vbar, κμ, νμ, κv, νv, αbar, λ, ρ, γ, ψ) = m
  (; μ, v) = state
  (; p, pμ_up, pμ_down, pv_up, pv_down, pμμ, pμv, pvv) = y

  # drift and volatility of μ, ν, p
  μc = μ
  σc = sqrt(v)
  μμ = κμ * (μbar - μ)
  σμ = νμ * sqrt(v)
  μv = κv * (vbar - v)

  pμ = (μμ >= 0) ? pμ_up : pμ_down
  pv = (νμ >= 0) ? pv_up : pv_down 
  σv = νv * sqrt(v) 
  σp_Zμ = pμ / p * σμ
  σp_Zv = pv / p * σv
  σp2 = σp_Zμ^2 + σp_Zv^2
  μp = pμ / p * μμ + pv / p * μv + 0.5 * pμμ / p * σμ^2 + 0.5 * pvv / p * σv^2

  # Market Price of Risk κ
  αstar = αbar - sqrt(λ * αbar / (0.5 * γ * (σc^2 + σp_Zμ^2 + σp_Zv^2)))
  αstar = clamp(αstar, 0.0, 1.0)
  κ_Zc = αstar * γ * σc
  κ_Zμ = (αstar - (1 / γ - 1) / (ψ - 1)) * γ * σp_Zμ
  κ_Zv = (αstar - (1 / γ - 1) / (ψ - 1)) * γ * σp_Zv
  κ2 = κ_Zc^2 + κ_Zμ^2 + κ_Zv^2

  # Interest rate r
  #r = ρ + μc / ψ - (1 - 1 / ψ) / (2 * γ) * κ2 - (1 / ψ - γ) / (2 * γ * (ψ - 1)) * σp2 - (1 - γ) / (γ * ψ) * (κμ * σμ + κv * σv) + 1 / ψ  * (κ_Zc * σc + κμ * σμ + κv * σv)
  #out = p * (1 / p + μc + μp - r - κ_Zc * σc - κ_Zμ * σp_Zμ - κ_Zv * σp_Zv)
  pt = - p * (1 / p - ρ + (1 - 1 / ψ) * (μc - 0.5 * γ * σc^2 * (1 - (αstar - 1)^2)) + μp + (0.5 * (1 / ψ - γ) / (1 - 1 / ψ) + 0.5 * γ * (1 - 1 / ψ) * (αstar - 1)^2) * σp2)
  return (; pt)
end

m = HaddadModel()
# initialize grid of μ and v
σ = sqrt(m.νμ^2 * m.vbar / (2 * m.κμ))
μmin = quantile(Normal(m.μbar, σ), 0.001)
μmax = quantile(Normal(m.μbar, σ), 0.999)
μs = range(μmin,  μmax, length = 30)
α = 2 * m.κv * m.vbar / m.νv^2
β = m.νv^2 / (2 * m.κv)
vmin = quantile(Gamma(α, β), 0.001)
vmax = quantile(Gamma(α, β), 0.999)
vs = range(vmin,  vmax, length = 30)
stategrid = OrderedDict(:μ => μs, :v => vs)

# guess for price-dividend ratio
yend =   OrderedDict(:p => ones(length(stategrid[:μ]), length(stategrid[:v])))

# solve
result = pdesolve(m, stategrid, yend)
@assert result.residual_norm <= 1e-5