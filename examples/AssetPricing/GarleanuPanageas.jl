# Garleanu Panageas (2015) JPE "Young, Old, Conservative, and Bold: The Implications of Heterogeneity and Finite Lives for Asset Pricing"

using EconPDEs

Base.@kwdef struct GarleanuPanageasModel
  # utility function
  γA::Float64 = 1.5
  ψA::Float64 = 0.7
  γB::Float64 = 10.0
  ψB::Float64 = 0.05
  ρ::Float64 = 0.001
  δ::Float64 = 0.02

  # proportion a
  νA::Float64 = 0.01

  # consumption
  μ::Float64 = 0.02
  σ::Float64 = 0.041

  # earning function
  B1::Float64 = 30.72
  δ1::Float64 = 0.0525
  B2::Float64 = -30.29
  δ2::Float64 = 0.0611
  ω::Float64 = 0.92
end



function (m::GarleanuPanageasModel)(state::NamedTuple, y::NamedTuple)
  # pA is wealth / consumption ratio of agent A
  # pB is wealth / consumption ratio of agent B
  # ϕ1 is value of claim that promises 1 today and then grows at rate μ - δ- δ1
  # ϕ2 is value of claim that promises 1 today and then grows at rate μ - δ- δ2
  (; γA, ψA, γB, ψB, ρ, δ, νA, μ, σ, B1, δ1, B2, δ2, ω) = m  
  (; x) = state
  (; pA, pAx_up, pAx_down, pAxx, pB, pBx_up, pBx_down, pBxx, ϕ1, ϕ1x_up, ϕ1x_down, ϕ1xx, ϕ2, ϕ2x_up, ϕ2x_down, ϕ2xx) = y

  scale = δ / (δ + δ1) * B1 + δ / (δ + δ2) * B2
  B1 = B1 / scale
  B2 = B2 / scale

  # Market price of risk κ
  pAx, pBx, ϕ1x, ϕ2x = pAx_up, pBx_up, ϕ1x_up, ϕ2x_up
  iter = 0
  @label start
  Γ = 1 / (x / γA + (1 - x) / γB)
  p = x * pA + (1 - x) * pB
  σx = σ * x * (Γ / γA - 1) / (1 + Γ * x * (1 - x) / (γA * γB) * ((1 - γB * ψB) / (ψB - 1) * (pBx / pB) - (1 - γA * ψA) / (ψA - 1) * (pAx / pA)))
  σpA = pAx / pA * σx
  σpB = pBx / pB * σx
  σϕ1 = ϕ1x / ϕ1 * σx
  σϕ2 = ϕ2x / ϕ2 * σx
  κ = Γ * (σ - x * (1 - γA * ψA) / (γA * (ψA - 1)) * σpA - (1 - x) * (1 - γB * ψB) / (γB * (ψB - 1)) * σpB)
  σCA = κ / γA + (1 - γA * ψA) / (γA * (ψA - 1)) * σpA
  σCB = κ / γB + (1 - γB * ψB) / (γB * (ψB - 1)) * σpB
  # A.16 Equation in Garleanu Panageas has a typo
  mcA = κ^2 * (1 + ψA) / (2 * γA) + (1 - ψA * γA) / (γA * (ψA - 1)) * κ * σpA - (1 - γA * ψA) / (2 * γA * (ψA - 1)) * σpA^2
  mcB = κ^2 * (1 + ψB) / (2 * γB) + (1 - ψB * γB) / (γB * (ψB - 1)) * κ * σpB - (1 - γB * ψB) / (2 * γB * (ψB - 1)) * σpB^2
  r =  ρ + 1 / (ψA * x  + ψB * (1 - x))  * (μ - x * mcA - (1 - x) * mcB - δ * ((νA / pA + (1 - νA) / pB) * ω * (B1 * ϕ1 + B2 * ϕ2) - 1))
  μCA = ψA * (r - ρ) + mcA
  μCB = ψB * (r - ρ) + mcB
  μx = x * (μCA - μ) + δ * (νA / pA * ω * (B1 * ϕ1 + B2 * ϕ2) - x) - σ * σx  
  # upwinding
  if (iter == 0) & (μx <= 0)
    iter += 1
    pAx, pBx, ϕ1x, ϕ2x = pAx_down, pBx_down, ϕ1x_down, ϕ2x_down
    @goto start
  end

  μpA = pAx / pA * μx + 0.5 * pAxx / pA * σx^2
  μpB = pBx / pB * μx + 0.5 * pBxx / pB * σx^2
  μϕ1 = ϕ1x / ϕ1 * μx + 0.5 * ϕ1xx / ϕ1 * σx^2
  μϕ2 = ϕ2x / ϕ2 * μx + 0.5 * ϕ2xx / ϕ2 * σx^2
  
  # Market Pricing
  pAt = - pA * (1 / pA + μCA + μpA + σCA * σpA - r - δ - κ * (σpA + σCA))
  pBt = - pB * (1 / pB + μCB + μpB + σCB * σpB - r - δ - κ * (σpB + σCB))
  ϕ1t = - ϕ1 * (1 / ϕ1 + μ - δ - δ1 + μϕ1 + σ * σϕ1 - r - κ * (σϕ1 + σ))
  ϕ2t = - ϕ2 * (1 / ϕ2 + μ - δ - δ2 + μϕ2 + σ * σϕ2 - r - κ * (σϕ2 + σ))

  return (; pAt, pBt, ϕ1t, ϕ2t), (;r)
end

m = GarleanuPanageasModel()
stategrid = OrderedDict(:x => range(0.0, 1.0, length = 100))
yend = OrderedDict(:pA => ones(length(stategrid[:x])), :pB => ones(length(stategrid[:x])), :ϕ1 => ones(length(stategrid[:x])), :ϕ2 => ones(length(stategrid[:x])))
result = pdesolve(m, stategrid, yend)
@assert result.residual_norm <= 1e-5
