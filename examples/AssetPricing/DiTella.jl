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

  pAx, pBx, px = pAx_up, pBx_up, px_up
  iter = 0
  @label start
  σX = x * (1 - x) * (1 - γ) / (γ * (ψ - 1)) * (pAν / pA - pBν / pB) * σν / (1 - x * (1 - x) * (1 - γ) / (γ * (ψ - 1)) * (pAx / pA - pBx / pB))
  σpA = pAx / pA * σX + pAν / pA * σν
  σpB = pBx / pB * σX + pBν / pB * σν
  σp = px / p * σX + pν / p * σν
  κ = (σp + σ - (1 - γ) / (γ * (ψ - 1)) * (x * σpA + (1 - x) * σpB)) / (1 / γ)
  κν = γ * ϕ * ν / x
  σA = κ / γ + (1 - γ) / (γ * (ψ - 1)) * σpA
  νA = κν / γ
  σB = κ / γ + (1 - γ) / (γ * (ψ - 1)) * σpB
  # Interest rate r
  μX = x * (1 - x) * ((σA * κ + νA * κν - 1 / pA - τ) - (σB * κ -  1 / pB + τ * x / (1 - x)) - (σA - σB) * (σ + σp))
  if (iter == 0) & (μX <= 0)
    iter += 1
    pAx, pBx, px = pAx_down, pBx_down, px_down
    @goto start
  end

  μpA = pAx / pA * μX + pAν / pA * μν + 0.5 * pAxx / pA * σX^2 + 0.5 * pAνν / pA * σν^2 + pAxν / pA * σX * σν
  μpB = pBx / pB * μX + pBν / pB * μν + 0.5 * pBxx / pB * σX^2 + 0.5 * pBνν / pB * σν^2 + pBxν / pB * σX * σν
  μp = px / p * μX + pν / p * μν + 0.5 * pxx / p * σX^2 + 0.5 * pνν / p * σν^2 + pxν / p * σX * σν
  r = (1 - i) / p + g + μp + σ * σp - κ * (σ + σp) - γ / x * (ϕ * ν)^2

  # Market Pricing
  pAt = - pA * (1 / pA  + (ψ - 1) * τ / (1 - γ) * ((pA / pB)^((1 - γ) / (1 - ψ)) - 1) - ψ * ρ + (ψ - 1) * (r + κ * σA + κν * νA) + μpA - (ψ - 1) * γ / 2 * (σA^2 + νA^2) + (2 - ψ - γ) / (2 * (ψ - 1)) * σpA^2 + (1 - γ) * σpA * σA)
  pBt = - pB * (1 / pB - ψ * ρ + (ψ - 1) * (r + κ * σB) + μpB - (ψ - 1) * γ / 2 * σB^2 + (2 - ψ - γ) / (2 * (ψ - 1)) * σpB^2 + (1 - γ) * σpB * σB)
  # algebraic constraint
  pt = - p * ((1 - i) / p - x / pA - (1 - x) / pB)
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