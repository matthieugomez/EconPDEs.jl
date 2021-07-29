using EconPDEs, Distributions

Base.@kwdef mutable struct DiTellaModel
  # Utility Function
  γ::Float64 = 5.0 
  ψ::Float64 = 1.5
  ρ::Float64 = 0.05
  τ::Float64 = 0.4

  # Technology
  A::Float64 = 200.0
  σ::Float64 = 0.03

  # MoralHazard
  ϕ::Float64 = 0.2

  # Idiosyncratic
  νbar::Float64 = 0.24
  κν::Float64 = 0.22
  σνbar::Float64 = -0.13
end


function (m::DiTellaModel)(state::NamedTuple, y::NamedTuple)
  (; γ, ψ, ρ, τ, A, σ, ϕ, νbar, κν, σνbar) = m  
  (; x, ν) = state
  (; pA, pAx_up, pAx_down, pAν_up, pAν_down, pAxx, pAxν, pAνν, pB, pBx_up, pBx_down, pBν_up, pBν_down, pBxx, pBxν, pBνν, p, px_up, px_down, pν_up, pν_down, pxx, pxν, pνν) = y

  # drift and volatility of state variable ν
  g = p / (2 * A)
  i = A * g^2
  μν = κν * (νbar - ν)
  σν = σνbar * sqrt(ν)
  pAν = (μν >= 0) ? pAν_up : pAν_down
  pBν = (μν >= 0) ? pBν_up : pBν_down
  pν = (μν >= 0) ? pν_up : pν_down

 # Market price of risk κ
  σX_up = x * (1 - x) * (1 - γ) / (γ * (ψ - 1)) * (pAν / pA - pBν / pB) * σν / (1 - x * (1 - x) * (1 - γ) / (γ * (ψ - 1)) * (pAx_up / pA - pBx_up / pB))
  σX_down = x * (1 - x) * (1 - γ) / (γ * (ψ - 1)) * (pAν / pA - pBν / pB) * σν / (1 - x * (1 - x) * (1 - γ) / (γ * (ψ - 1)) * (pAx_down / pA - pBx_down / pB))
  σpA_up = pAx_up / pA * σX_up + pAν / pA * σν
  σpA_down = pAx_down / pA * σX_down + pAν / pA * σν
  σpB_up = pBx_up / pB * σX_up + pBν / pB * σν
  σpB_down = pBx_down / pB * σX_down + pBν / pB * σν
  σp_up = px_up / p * σX_up + pν / p * σν
  σp_down = px_down / p * σX_down + pν / p * σν
  κ_up = (σp_up + σ - (1 - γ) / (γ * (ψ - 1)) * (x * σpA_up + (1 - x) * σpB_up)) / (1 / γ)
  κ_down = (σp_down + σ - (1 - γ) / (γ * (ψ - 1)) * (x * σpA_down + (1 - x) * σpB_down)) / (1 / γ)
  κν = γ * ϕ * ν / x
  σA_up = κ_up / γ + (1 - γ) / (γ * (ψ - 1)) * σpA_up
  σA_down = κ_down / γ + (1 - γ) / (γ * (ψ - 1)) * σpA_down
  νA = κν / γ
  σB_up = κ_up / γ + (1 - γ) / (γ * (ψ - 1)) * σpB_up
  σB_down = κ_down / γ + (1 - γ) / (γ * (ψ - 1)) * σpB_down

  # Interest rate r
  μX_up = x * (1 - x) * ((σA_up * κ_up + νA * κν - 1 / pA - τ) - (σB_up * κ_up -  1 / pB + τ * x / (1 - x)) - (σA_up - σB_up) * (σ + σp_up))
  μX_down = x * (1 - x) * ((σA_down * κ_down + νA * κν - 1 / pA - τ) - (σB_down * κ_down -  1 / pB + τ * x / (1 - x)) - (σA_down - σB_down) * (σ + σp_down))

  pAx = (μX_up >= 0) ? pAx_up : pAx_down
  pBx = (μX_up >= 0) ? pBx_up : pBx_down
  px =  (μX_up >= 0) ? px_up : px_down
  σpA = (μX_up >= 0) ? σpA_up : σpA_down
  σpB = (μX_up >= 0) ? σpB_up : σpB_down
  σp =  (μX_up >= 0) ? σp_up : σp_down
  σA =  (μX_up >= 0) ? σA_up : σA_down
  σB =  (μX_up >= 0) ? σB_up : σB_down
  κ =   (μX_up >= 0) ? κ_up : κ_down
  σX =  (μX_up >= 0) ? σX_up : σX_down
  μX =  (μX_up >= 0) ? μX_up : μX_down

  μpA = pAx / pA * μX + pAν / pA * μν + 0.5 * pAxx / pA * σX^2 + 0.5 * pAνν / pA * σν^2 + pAxν / pA * σX * σν
  μpB = pBx / pB * μX + pBν / pB * μν + 0.5 * pBxx / pB * σX^2 + 0.5 * pBνν / pB * σν^2 + pBxν / pB * σX * σν
  μp = px / p * μX + pν / p * μν + 0.5 * pxx / p * σX^2 + 0.5 * pνν / p * σν^2 + pxν / p * σX * σν
  r = (1 - i) / p + g + μp + σ * σp - κ * (σ + σp) - γ / x * (ϕ * ν)^2

  # Market Pricing
  pAt = - pA * (1 / pA  + (ψ - 1) * τ / (1 - γ) * ((pA / pB)^((1 - γ) / (1 - ψ)) - 1) - ψ * ρ + (ψ - 1) * (r + κ * σA + κν * νA) + μpA - (ψ - 1) * γ / 2 * (σA^2 + νA^2) + (2 - ψ - γ) / (2 * (ψ - 1)) * σpA^2 + (1 - γ) * σpA * σA)
  pBt = - pB * (1 / pB - ψ * ρ + (ψ - 1) * (r + κ * σB) + μpB - (ψ - 1) * γ / 2 * σB^2 + (2 - ψ - γ) / (2 * (ψ - 1)) * σpB^2 + (1 - γ) * σpB * σB)
  # algebraic constraint
  pt = p * ((1 - i) / p - x / pA - (1 - x) / pB)

  return (; pAt, pBt, pt)
end

m = DiTellaModel()
xn = 80
νn = 10
xs = range(0.01, 0.99, length = xn)
distribution = Gamma(2 * m.κν * m.νbar / m.σνbar^2, m.σνbar^2 / (2 * m.κν))
νs = range(quantile(distribution, 0.001), quantile(distribution, 0.999), length = νn)
stategrid = OrderedDict(:x => xs, :ν => νs)
yend = OrderedDict(:pA => ones(xn, νn), :pB => ones(xn, νn), :p => ones(xn, νn))
y, residual_norm = pdesolve(m, stategrid, yend; is_algebraic = OrderedDict(:pA => false, :pB => false, :p => true))