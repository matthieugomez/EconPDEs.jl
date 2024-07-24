# He Krishnamurthy "Intermediary Asset Pricing" AER (2013)

using EconPDEs

Base.@kwdef mutable struct HeKrishnamurthy
  # See Table 2—Parameters and Targets
  m::Float64 = 4  # intermediation multiplier
  λ::Float64 = 0.6 # Debt/assets of intermediary sector
  g::Float64 = 0.02 # dividend growth
  σ::Float64 = 0.09 # dividend volatility 
  ρ::Float64 = 0.04 # time discount rate
  γ::Float64 = 2.0 # risk aversion specialist
  l::Float64 = 1.84 # labor income ratio
end

# Description of the model
# Households cannot invest directly in the risky asset market: can only invest in bonds and intermediaries. Households unwilling to invest more than m * w_t in the equity of intermediary, where w_t is networth intermediary (constant on outside equity financing)
# Specialists have the knowledge to invest in the risky assets, and can invest in the risky asset on behalf of the households
# Specialists accept household funds and allocate between risky and riskless. networth of specialist is w and wealth allocated by households is H.
# Denote αI ratio of risky asset holdings to total capital w + H

# households have log utility. Assume that they die continuously so that we can treat financial wealth as the thigns to maximize in t + dt (instead of total wealth), so that we can treat labor income easily.
# We make assumptions so that the household sector chooses to keep a minimum of λ w^h in short-term deb issued by intermediary. Formally, we assume that a fraction λ can only ever invest in riskless bond while remaining can invest in intermediation market 

# Denote pS the wealth-to-consumption ratio of specialists
function (m::HeKrishnamurthy)(state::NamedTuple, y::NamedTuple)
  (; m, λ, g, σ, ρ, γ, l) = m
  (; pS, pSx_up, pSx_down, pSxx) = y
  (; x) = state
  iter = 0
  pSx = pSx_up
  @label start
  # market clearing for goods gives x / pS + (1 - x) * ρ = (1 + l) / p
  # This allows to obtain p (and its derivatives) in terms of pS (and its derivatives)
  p = (1 + l) / (x / pS + (1 - x) * ρ)
  px = - (1 + l) / (x / pS + (1 - x) * ρ)^2 * (1 / pS - ρ - x / pS^2 * pSx)
  pxx = - (1 + l) / (x / pS + (1 - x) * ρ)^2 * (- 1 / pS^2 * pSx - 1 / pS^2 * pSx - x / pS^2 * pSxx + 2 * x / pS^3 * pSx^2) + 2 * (1 + l) / (x / pS + (1 - x) * ρ)^3 * (1 / pS - ρ - x / pS^2 * pSx)^2
  @assert x / pS + (1 - x) * ρ ≈ (1 + l) / p


  αI = 1 / (1 - λ * (1 - x))
  if x * (1 + m) * αI < 1
    αI = 1 / (x * (1 + m))
  end
  # we have the usual feeback loop between asset prices and wealth
  # σx = x * (αI - 1) * (σ + σp)
  σx = x * (αI - 1) * σ / (1 -  x * (αI - 1) * px / p)
  σpS = pSx / pS * σx
  σp = px / p * σx
  @assert σx ≈ x * (αI - 1) * (σ + σp)
  σR = σ + σp
  κ = γ * (αI * σR - σpS)
  μx = x * (1 - x) * ((αI - (1 - αI * x) / (1-x)) * κ * σR - 1 / pS + ρ - (αI - (1 - αI * x) / (1 - x)) * σR^2 - l / ((1 - x) * p))
  if (iter == 0) && (μx < 0)
    iter += 1
    pSx = pSx_down
    @goto start
  end
  μpS = pSx / pS * μx + 0.5 * pSxx / pS * σx^2
  μp = px / p * μx + 0.5 * pxx / p * σx^2
  r = 1 / p + g + μp + σ * σp - κ * σR
  σCS = κ / γ
  μCS = (r - ρ) / γ + (1 + 1 / γ) / (2 * γ) * κ^2 
  pSt = - (1 / pS + μCS + μpS + σCS * σpS - r - κ * (σCS + σpS))
  μR = r + κ * σR
  σI = αI * σR
  (; pSt), (; pS, p, r, κ, σp, μR, σR, αI, μx, σx, x, px, pxx, σI)
end

m = HeKrishnamurthy()
xn = 100
stategrid =  OrderedDict(:x => range(0, 1, length = xn+2)[2:(end-1)].^1.5)
yend = OrderedDict(:pS => ones(length(stategrid[:x])))
y, residual_norm, a = pdesolve(m, stategrid, yend)
@assert residual_norm <= 1e-5