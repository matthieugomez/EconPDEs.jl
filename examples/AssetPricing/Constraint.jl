using EconPDEs

struct ConstraintModel

  # utility function
  γ::Float64 
  ψ::Float64
  ρ::Float64
  δ::Float64

  # FixedShare
  α::Float64
  # Proportion Fixed Share
  ν::Float64

  # consumption
  μ::Float64
  σ::Float64
end

function ConstraintModel(;γ  = 1.5, ψ = 0.7, ρ = 0.001, δ = 0.02, α = 1.0, ν = 0.01, μ = 0.02, σ = 0.041)
  ConstraintModel(γ , ψ, γ, ψ, ρ, δ, ν, μ, σ)
end

function initialize_stategrid(m::ConstraintModel; n = 200)
  OrderedDict(:x => range(0.0, stop = 1.0, length = n))
end

function initialize_y(m::ConstraintModel, stategrid::OrderedDict)
    x = ones(length(stategrid[:x]))
    OrderedDict(:pA => x, :pB => x)
end

# p = ρ / γ + (1- 1/γ) * (r + α * κ + α σ^2)

function (m::ConstraintModel)(state::NamedTuple, y::NamedTuple)
  γ = m.γ ; ψ = m.ψ ; ρ = m.ρ ; δ = m.δ ; ν = m.ν ; μ = m.μ ; σ = m.σ
  x = state.x
  # pA is wealth / consumption ratio of agent A
  # pB is wealth / consumption ratio of agent B
  pA = y.pA ; pAx = y.pAx ; pAxx = y.pAxx ; pB = y.pB ; pBx = y.pBx ; pBxx = y.pBxx

  # Market price of risk κ
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

  # Interest rate r
  # A.16 Equation in Garleanu Panageas has a typo
  mcA = κ^2 * (1 + ψA) / (2 * γA) + (1 - ψA * γA) / (γA * (ψA - 1)) * κ * σpA - (1 - γA * ψA) / (2 * γA * (ψA - 1)) * σpA^2
  mcB = κ^2 * (1 + ψB) / (2 * γB) + (1 - ψB * γB) / (γB * (ψB - 1)) * κ * σpB - (1 - γB * ψB) / (2 * γB * (ψB - 1)) * σpB^2
  r =  ρ + 1 / (ψA * x  + ψB * (1 - x))  * (μ - x * mcA - (1 - x) * mcB - δ * ((νA / pA + (1 - νA) / pB) * (ϕ1 + ϕ2) - 1))
  μCA = ψA * (r - ρ) + mcA
  μCB = ψB * (r - ρ) + mcB
  μx = x * (μCA - δ - μ) + δ * νA / pA * (ϕ1 + ϕ2) - σ * σx  
  μpA = pAx / pA * μx + 0.5 * pAxx / pA * σx^2
  μpB = pBx / pB * μx + 0.5 * pBxx / pB * σx^2
  μϕ1 = ϕ1x / ϕ1 * μx + 0.5 * ϕ1xx / ϕ1 * σx^2
  μϕ2 = ϕ2x / ϕ2 * μx + 0.5 * ϕ2xx / ϕ2 * σx^2
  
  # Market Pricing
  pAt = pA * (1 / pA + (μCA - δ) + μpA + σCA * σpA - r - κ * (σpA + σCA))
  pBt = pB * (1 / pB + (μCB - δ) + μpB + σCB * σpB - r - κ * (σpB + σCB))
  ϕ1t = ϕ1 * (B1 * ω / ϕ1 + (μ - δ - δ1) + μϕ1 + σ * σϕ1 - r - κ * (σϕ1 + σ))
  ϕ2t = ϕ2 * (B2 * ω / ϕ2 + (μ - δ - δ2) + μϕ2 + σ * σϕ2 - r - κ * (σϕ2 + σ))

  return (pAt, pBt, ϕ1t, ϕ2t), (μx, ), (μx = μx, p = p, pA = pA, pB = pB, κ = κ, r = r, σx = σx)
end

m = ConstraintModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
y, result, distance = pdesolve(m, stategrid, y0)