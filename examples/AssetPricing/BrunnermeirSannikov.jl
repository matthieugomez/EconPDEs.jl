# Brunnermeier Sannikov "A Macroeconomic Model with a Financial Sector" (2014)
# The idfference with Ditella is that Ψ appears in the last algebratic equation, making it much much harder to solve for it
# using non linear solver
using EconPDEs

Base.@kwdef mutable struct BrunnermeierSannikov
  # See Table H Numerical examples
  ρE::Float64 = 0.06
  ρH::Float64 = 0.06
  γ::Float64 = 2.0 
  ψ::Float64 = 0.5 
  δ::Float64 = 0.05
  σ::Float64 = 0.1
  κ::Float64 = 10.0
  a::Float64 = 0.11
  alow::Float64 = 0.1
  χlow::Float64 = 1.0
  firstregion::Bool = true
  q::Float64 = 1.0
  qx::Float64 = 0.0
end

function (m::BrunnermeierSannikov)(state::NamedTuple, y::NamedTuple, Δx, xmin)
  (; ρE, ρH, γ, ψ, δ, σ, κ, a, alow, χlow) = m
  (; x) = state
  (; pE, pEx_up, pEx_down, pExx, pH, pHx_up, pHx_down, pHxx) = y
  # χ is the fraction of risk held by experts.
  # ψ is the fraction of (real) capital held by experts
  # so χ * ψ / x corresponds to total risk held by experts  / their networth
  #   σE = κE / γ + (1 - γ) / (γ * (ψ - 1)) * σpE
  # implies κE = γ * σE -  (1 - γ) / ((ψ - 1)) * σpE
  # so κE - κH = γ * (σE -σH) -  (1 - γ) / ((ψ - 1)) * (σpE -σpH)
  #   = γ * (χ * ψ / x - (1-χ * ψ) / (1-x) - (1 - γ) / (γ * (ψ - 1)) * (σpE - σpH)
  # . Specifically, we suppose that experts must retain at least a fraction χbar ∈ (0, 1] of equity. (i
  # drift and volatility of state variable ν
  #@show x, p, ψ
  pE = max(pE, 1e-3)
  pH = max(pH, 1e-3)  



  χ = max(χlow, x)
  # market clear givges
  # x / pE + (1-x) / pH = (ψ * a + (1-ψ) * alow - i) / q
  # this implies ψ
  pEx, pHx = pEx_up, pHx_up
  iter = 0 
  @label start
  q = m.q
  if x ≈ xmin
    m.firstregion = true
    # poin toutside the gri, so zero
    Ψ = 0.0
    q = (Ψ * a + (1 - Ψ) * alow + 1 / κ) / ((x / pE + (1 - x) / pH) + 1 / κ)
  end
  if m.firstregion
    # in the second region (a - alow) / q = (κE - κH) * (σ + σq)
    # here q corresponds to value in the previous x grif
    i = (q - 1) / κ
    Ψ =  max(((x / pE + (1 - x) / pH) * q + i - alow) / (a - alow), (1e-10 + x) / χ)
    K = q / (a - alow) * χ * σ^2 * (γ / (x * (1 - x)) - (1 - γ) / (ψ - 1) * (pEx / pE - pHx / pH))
    X =χ * Ψ - x
    qx = (1 - sqrt(max(eps(), X * K))) / X * q
    qx = max(min(qx, 4.0), 0.0)
    q = q + qx * Δx
    #@assert abs((1 - qx / q * X)^2 - K * X) < 1e-8
    if Ψ > 1.0
      m.firstregion = false
      Ψ = 1
      q = (a + 1 / κ) / ((x / pE + (1 - x) / pH) + 1 / κ)
      qx = - (a + 1 / κ) / ((x / pE + (1 - x) / pH) + 1 / κ)^2 * (1 / pE - 1 / pH - x / pE^2 * pEx - (1-x) / pH^2 * pHx)
    end
  else
    # in the second region Ψ = 1
    Ψ = 1
    # We get q through market clearing equation after plugging i = (p - 1) / κ
    q = (a + 1 / κ) / ((x / pE + (1 - x) / pH) + 1 / κ)
    qx = - (a + 1 / κ) / ((x / pE + (1 - x) / pH) + 1 / κ)^2 * (1 / pE - 1 / pH - x / pE^2 * pEx - (1-x) / pH^2 * pHx)
    i = (q - 1) / κ
    @assert abs((a - i) / q - (x / pE + (1 - x) / pH)) < 1e-8
  end
  #@show m.firstregion, x, qx
  #@assert x isa Float64
  qxx = (qx - m.qx) / Δx
  m.q = q
  m.qx = qx
  q = max(1e-3, q)
  i = (q - 1) / κ
  Φ = log(q) / κ
  σx = (χ * Ψ - x) * σ / (1 - (χ * Ψ - x) * qx / q)
  σpE = pEx / pE * σx
  σpH = pHx / pH * σx
  σq = qx / q * σx
  σE = Ψ * χ / x * (σ + σq)
  σH = (1 - Ψ * χ) / (1 - x) * (σ + σq)
  κE = γ * σE - (1 - γ) / (ψ - 1) * σpE
  κH = γ * σH - (1 - γ) / (ψ - 1) * σpH
  μx = x * ((a - i) / q - 1 / pE) + σx * (κE - σ - σq) + x * (σ + σq) * (1 - χ) * (κE - κH)
  #μx = x * (1 - x) * (κE * σE - κH * σH + 1 / pH - 1 / pE - (σE - σH) * (σ + σq))
  if (iter == 0) && (μx < 0)
    iter += 1
    pEx, pHx = pEx_down, pHx_down
    @goto start
  end  
  μpE = pEx / pE * μx + 0.5 * pExx / pE * σx^2 
  μpH = pHx / pH * μx + 0.5 * pHxx / pH * σx^2 
  μq = qx / q * μx + 0.5 * qxx / q * σx^2

  r = min(max((a - i) / q + Φ - δ + μq + σ * σq - (χ * κE + (1 - χ) * κH) * (σ + σq), -0.1), 0.1)
  #if ψ < 1
    # @assert  r ≈ (alow - i) / p + Φ - δ + μp + σ * σp - κH * (σ + σp)
  #end
  # Market Pricing
  pEt = -  pE * (1 / pE - ψ * ρE + (ψ - 1) * (r + κE * σE) + μpE - (ψ - 1) * γ / 2 * σE^2 + (2 - ψ - γ) / (2 * (ψ - 1)) * σpE^2 + (1 - γ) * σpE * σE)
  pHt = -  pH * (1 / pH - ψ * ρH + (ψ - 1) * (r + κH * σH) + μpH - (ψ - 1) * γ / 2 * σH^2 + (2 - ψ - γ) / (2 * (ψ - 1)) * σpH^2 + (1 - γ) * σpH * σH)
  (; pEt, pHt), (; x, r, μx, σx, Ψ, κ, κE, κH, σE, σH, σq, pEt, pHt, q, qx, qxx, μq)
end



