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
  (; γA, ψA, γB, ψB, ρ, δ, νA, μ, σ, B1, δ1, B2, δ2, ω) = m  
  (; x) = state
  (; pA, pAx_up, pAx_down, pAxx, pB, pBx_up, pBx_down, pBxx, ϕ1, ϕ1x_up, ϕ1x_down, ϕ1xx, ϕ2, ϕ2x_up, ϕ2x_down, ϕ2xx) = y

  scale = δ / (δ + δ1) * B1 + δ / (δ + δ2) * B2
  B1 = B1 / scale
  B2 = B2 / scale
  # pA is wealth / consumption ratio of agent A
  # pB is wealth / consumption ratio of agent B
  # ϕ1 is value of claim that promises B_1ωexp(-(δ+δ1)(s-t))C_s/C_t for s ≥ t
  # ϕ2 is value of claim that promises B_2ωexp(-(δ+δ2)(s-t))C_s/C_t for s ≥ t
  # Market price of risk κ
  Γ = 1 / (x / γA + (1 - x) / γB)
  p = x * pA + (1 - x) * pB
  σx_up = σ * x * (Γ / γA - 1) / (1 + Γ * x * (1 - x) / (γA * γB) * ((1 - γB * ψB) / (ψB - 1) * (pBx_up / pB) - (1 - γA * ψA) / (ψA - 1) * (pAx_up / pA)))
  σx_down = σ * x * (Γ / γA - 1) / (1 + Γ * x * (1 - x) / (γA * γB) * ((1 - γB * ψB) / (ψB - 1) * (pBx_down / pB) - (1 - γA * ψA) / (ψA - 1) * (pAx_down / pA)))
  σpA_up = pAx_up / pA * σx_up
  σpA_down = pAx_down / pA * σx_down
  σpB_up = pBx_up / pB * σx_up
  σpB_down = pBx_down / pB * σx_down
  σϕ1_up = ϕ1x_up / ϕ1 * σx_up
  σϕ1_down = ϕ1x_down / ϕ1 * σx_down
  σϕ2_up = ϕ2x_up / ϕ2 * σx_up
  σϕ2_down = ϕ2x_down / ϕ2 * σx_down
  κ_up = Γ * (σ - x * (1 - γA * ψA) / (γA * (ψA - 1)) * σpA_up - (1 - x) * (1 - γB * ψB) / (γB * (ψB - 1)) * σpB_up)
  κ_down = Γ * (σ - x * (1 - γA * ψA) / (γA * (ψA - 1)) * σpA_down - (1 - x) * (1 - γB * ψB) / (γB * (ψB - 1)) * σpB_down)

  σCA_up = κ_up / γA + (1 - γA * ψA) / (γA * (ψA - 1)) * σpA_up
  σCA_down = κ_down / γA + (1 - γA * ψA) / (γA * (ψA - 1)) * σpA_down
  σCB_up = κ_up / γB + (1 - γB * ψB) / (γB * (ψB - 1)) * σpB_up
  σCB_down = κ_down / γB + (1 - γB * ψB) / (γB * (ψB - 1)) * σpB_down

  # Interest rate r
  # A.16 Equation in Garleanu Panageas has a typo
  mcA_up = κ_up^2 * (1 + ψA) / (2 * γA) + (1 - ψA * γA) / (γA * (ψA - 1)) * κ_up * σpA_up - (1 - γA * ψA) / (2 * γA * (ψA - 1)) * σpA_up^2
  mcA_down = κ_down^2 * (1 + ψA) / (2 * γA) + (1 - ψA * γA) / (γA * (ψA - 1)) * κ_down * σpA_down - (1 - γA * ψA) / (2 * γA * (ψA - 1)) * σpA_down^2

  mcB_up = κ_up^2 * (1 + ψB) / (2 * γB) + (1 - ψB * γB) / (γB * (ψB - 1)) * κ_up * σpB_up - (1 - γB * ψB) / (2 * γB * (ψB - 1)) * σpB_up^2
  mcB_down = κ_down^2 * (1 + ψB) / (2 * γB) + (1 - ψB * γB) / (γB * (ψB - 1)) * κ_down * σpB_down - (1 - γB * ψB) / (2 * γB * (ψB - 1)) * σpB_down^2
  r_up =  ρ + 1 / (ψA * x  + ψB * (1 - x))  * (μ - x * mcA_up - (1 - x) * mcB_up - δ * ((νA / pA + (1 - νA) / pB) * (ϕ1 + ϕ2) - 1))
  r_down =  ρ + 1 / (ψA * x  + ψB * (1 - x))  * (μ - x * mcA_down - (1 - x) * mcB_down - δ * ((νA / pA + (1 - νA) / pB) * (ϕ1 + ϕ2) - 1))
  μCA_up = ψA * (r_up - ρ) + mcA_up
  μCA_down = ψA * (r_down - ρ) + mcA_down
  μCB_up = ψB * (r_up - ρ) + mcB_up
  μCB_down = ψB * (r_down - ρ) + mcB_down
  μx_up = x * (μCA_up - δ - μ) + δ * νA / pA * (ϕ1 + ϕ2) - σ * σx_up  
  μx_down = x * (μCA_down - δ - μ) + δ * νA / pA * (ϕ1 + ϕ2) - σ * σx_down  

  pAx = (μx_up >= 0) ? pAx_up : pAx_down
  pBx = (μx_up >= 0) ? pBx_up : pBx_down
  ϕ1x =  (μx_up >= 0) ? ϕ1x_up : ϕ1x_down
  ϕ2x =  (μx_up >= 0) ? ϕ2x_up : ϕ2x_down
  σpA = (μx_up >= 0) ? σpA_up : σpA_down
  σpB = (μx_up >= 0) ? σpB_up : σpB_down
  σϕ1 =  (μx_up >= 0) ? σϕ1_up : σϕ1_down
  σϕ2 =  (μx_up >= 0) ? σϕ2_up : σϕ2_down
  σCA = (μx_up >= 0) ? σCA_up : σCB_down
  σCB = (μx_up >= 0) ? σCB_up : σCB_down
  μCA = (μx_up >= 0) ? μCA_up : μCB_down
  μCB = (μx_up >= 0) ? μCB_up : μCB_down
  κ =   (μx_up >= 0) ? κ_up : κ_down
  r =   (μx_up >= 0) ? r_up : r_down
  σx =  (μx_up >= 0) ? σx_up : σx_down
  μx =  (μx_up >= 0) ? μx_up : μx_down


  μpA = pAx / pA * μx + 0.5 * pAxx / pA * σx^2
  μpB = pBx / pB * μx + 0.5 * pBxx / pB * σx^2
  μϕ1 = ϕ1x / ϕ1 * μx + 0.5 * ϕ1xx / ϕ1 * σx^2
  μϕ2 = ϕ2x / ϕ2 * μx + 0.5 * ϕ2xx / ϕ2 * σx^2
  
  # Market Pricing
  pAt = pA * (1 / pA + (μCA - δ) + μpA + σCA * σpA - r - κ * (σpA + σCA))
  pBt = pB * (1 / pB + (μCB - δ) + μpB + σCB * σpB - r - κ * (σpB + σCB))
  ϕ1t = ϕ1 * (B1 * ω / ϕ1 + (μ - δ - δ1) + μϕ1 + σ * σϕ1 - r - κ * (σϕ1 + σ))
  ϕ2t = ϕ2 * (B2 * ω / ϕ2 + (μ - δ - δ2) + μϕ2 + σ * σϕ2 - r - κ * (σϕ2 + σ))

  return (pAt, pBt, ϕ1t, ϕ2t)
end

m = GarleanuPanageasModel()
xn = 200
stategrid = OrderedDict(:x => range(0.0, 1.0, length = xn))
yend = OrderedDict(:pA => ones(xn), :pB => ones(xn), :ϕ1 => ones(xn), :ϕ2 => ones(xn))
y, residual_norm = pdesolve(m, stategrid, yend)