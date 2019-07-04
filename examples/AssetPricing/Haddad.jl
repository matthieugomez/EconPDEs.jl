using EconPDEs, Distributions

struct HaddadModel
  # consumption process parameters
  μbar::Float64 
  vbar::Float64
  κμ::Float64 
  νμ::Float64 
  κv::Float64 
  νv::Float64 

  # active capital
  αbar::Float64
  λ::Float64

  # utility parameters
  ρ::Float64  
  γ::Float64 
  ψ::Float64
end

function HaddadModel(;μbar = 0.018, vbar = 0.00052, κμ = 0.3, νμ = 0.456, κv = 0.012, νv = 0.00472, αbar = 1.2, λ = 0.018, ρ = 0.0132, γ = 10.0, ψ = 1.5)
  HaddadModel(μbar, vbar, κμ, νμ, κv, νv, αbar, λ, ρ, γ, ψ)
end

function initialize_stategrid(m::HaddadModel; μn = 30, vn = 30)
  μbar = m.μbar; vbar = m.vbar ; κμ = m.κμ ; νμ = m.νμ ; κv = m.κv ; νv = m.νv ; αbar = m.αbar ; λ = m.λ ; ρ = m.ρ ; γ = m.γ ; ψ = m.ψ

  σ = sqrt(νμ^2 * vbar / (2 * κμ))
  μmin = quantile(Normal(μbar, σ), 0.025)
  μmax = quantile(Normal(μbar, σ), 0.975)
  μs = range(μmin, stop = μmax, length = μn)

  α = 2 * κv * vbar / νv^2
  β = νv^2 / (2 * κv)
  vmin = quantile(Gamma(α, β), 0.025)
  vmax = quantile(Gamma(α, β), 0.975)
  vs = range(vmin, stop = vmax, length = vn)

  OrderedDict(:μ => μs, :v => vs)
end
  
function initialize_y(m::HaddadModel, stategrid::OrderedDict)
  OrderedDict(:p => ones(length(stategrid[:μ]), length(stategrid[:v])))
end

function (m::HaddadModel)(state::NamedTuple, y::NamedTuple)
  μbar = m.μbar ; vbar = m.vbar ; κμ = m.κμ ; νμ = m.νμ ; κv = m.κv ; νv = m.νv ; αbar = m.αbar ; λ = m.λ ; ρ = m.ρ ; γ = m.γ ; ψ = m.ψ
  μ, v = state.μ, state.v
  p, pμ, pv, pμμ, pμv, pvv = y.p, y.pμ, y.pv, y.pμμ, y.pμv, y.pvv

  # drift and volatility of μ, ν, p
  μc = μ
  σc = sqrt(v)
  μμ = κμ * (μbar - μ)
  σμ = νμ * sqrt(v)
  μv = κv * (vbar - v)
  σv = νv * sqrt(v) 
  σp_Zμ = pμ / p * σμ
  σp_Zv = pv / p * σv
  μp = pμ / p * μμ + pv / p * μv + 0.5 * pμμ / p * σμ^2 + 0.5 * pvv / p * σv^2

  # Market Price of Risk κ
  αstar = αbar - sqrt(λ * αbar / (0.5 * γ * (σc^2 + σp_Zμ^2 + σp_Zv^2)))
  αstar = clamp(αstar, 0.0, 1.0)
  κ_Zc = αstar * γ * σc
  κ_Zμ = (αstar - (1 / γ - 1) / (ψ - 1)) * γ * σp_Zμ
  κ_Zv = (αstar - (1 / γ - 1) / (ψ - 1)) * γ * σp_Zv

  # Interest rate r
  κ2 = κ_Zc^2 + κ_Zμ^2 + κ_Zv^2
  σp2 = σp_Zμ^2 + σp_Zv^2
  r = ρ + μc / ψ - (1 - 1 / ψ) / (2 * γ) * κ2 - (1 / ψ - γ) / (2 * γ * (ψ - 1)) * σp2 - (1 - γ) / (γ * ψ) * (κμ * σμ + κv * σv) + 1 / ψ  * (κ_Zc * σc + κμ * σμ + κv * σv)

  # Market pricing
  out = p * (1 / p + μc + μp - r - κ_Zc * σc - κ_Zμ * σp_Zμ - κ_Zv * σp_Zv)

  #out = p * (1 / p - ρ + (1 - 1 / ψ) * (μc - 0.5 * γ * σc^2 * (1 - (αstar - 1)^2)) + μp + (0.5 * (1 / ψ - γ) / (1 - 1 / ψ) + 0.5 * γ * (1 - 1 / ψ) * (αstar - 1)^2) * σp2)

  return (out,) , (μμ, μv), (p = p, r = r, κ_Zc = κ_Zc, κ_Zμ = κ_Zμ, κ_Zv = κ_Zv, αstar = αstar)
end

m = HaddadModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
y, result, distance = pdesolve(m, stategrid, y0)