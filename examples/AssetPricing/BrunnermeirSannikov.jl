using EconPDEs, Roots


# There are households and experts, with different productivity
# Experts can issue equity but must retain at least a fraction χlow of their equity. This is similar to HK: households can only invest up  to w^E * m in intermediaries; with m = 1 / χlow - 1
# Ψ is the fraction of capital held by experts. Compared to HK, can be lower than one
# Moreover, experts can issue equity at lower required returns than the return on capital (ie. they earn management fee)
Base.@kwdef mutable struct BrunnermeierSannikov
  # Calibration uses Section 3.6 of the textbook chapter (I think r is typo for ρH)
  ρE::Float64 = 0.06
  ρH::Float64 = 0.05
  γ::Float64 = 2.0 
  ψ::Float64 = 0.5 
  δ::Float64 = 0.05
  σ::Float64 = 0.1
  κ::Float64 = 10.0
  a::Float64 = 0.11
  alow::Float64 = 0.03
  χlow::Float64 = 0.5
  q_old::Float64 = 1.0
  qx_old::Float64 = 0.0
  Ψ_old::Float64 = 1.0
end

function (m::BrunnermeierSannikov)(state::NamedTuple, y::NamedTuple, Δx, xmin)
  (; ρE, ρH, γ, ψ, δ, σ, κ, a, alow, χlow) = m
  (; x) = state
  (; pE, pEx_up, pEx_down, pExx, pH, pHx_up, pHx_down, pHxx) = y
  pE = max(pE, 1e-3)
  pH = max(pH, 1e-3)  
  χ = max(χlow, x)
  pEx, pHx = pEx_up, pHx_up
  iter = 0 
  @label start
  q_old = m.q_old
  qx_old = m.qx_old
  Ψ_old = m.Ψ_old

  # First, solve for q, qx, and qxx
  if abs(x - xmin) < 1e-9
    # initial condition at lower boundery of x grid (point outside the grid technically)
    # We must have that expert consumption approaches zero as xt approaches zero (linked to transversality condition)
    Ψ_old = 0.0
    q_old = (alow + 1 / κ) / (1 / pH + 1 / κ)
    qx_old = 0.0
  end
  if Ψ_old < 1.0
    # crisis region. 
    # in this case, we do not have Ψ = 1.0 but we have
    # (a - alow) / q = (κE - κH) * (σ + σq)
    # note that 
    out = find_zero(Ψ_old) do Ψ
        local q = (Ψ * a + (1 - Ψ) * alow + 1 / κ) / (x / pE + (1 - x) / pH + 1 / κ)
        local qx = (q - q_old) / Δx
        local αE = χ * Ψ / x
        local σx = x * (αE - 1) * σ / (1 - x * (αE - 1) * qx / q)
        local σpE = pEx / pE * σx
        local σpH = pHx / pH * σx
        local σq = qx / q * σx
        local σE = αE * (σ + σq)
        local αH = (1 - αE * x) / (1 - x)
        local σH = αH * (σ + σq)
        local κE = γ * (σE - (1 / γ - 1) / (ψ - 1) * σpE)
        local κH = γ * (σH - (1 / γ - 1) / (ψ - 1) * σpH)
        return (a - alow) / q -  χ * (κE - κH) * (σ + σq)
      end
    Ψ = out[1]
    q = (Ψ * a + (1 - Ψ) * alow + 1 / κ) / (x / pE + (1 - x) / pH + 1 / κ)
    qx = (q - q_old) / Δx
    if Ψ > 1
      Ψ_old = 1.0
    end
  end
  if Ψ_old ≥ 1.0
    # normal region. 
    Ψ = 1.0
    q = (a + 1 / κ) / ((x / pE + (1 - x) / pH) + 1 / κ)
    qx = - (a + 1 / κ) / ((x / pE + (1 - x) / pH) + 1 / κ)^2 * (1 / pE - 1 / pH - x / pE^2 * pEx - (1-x) / pH^2 * pHx)
  end
  qxx = (qx - qx_old) / Δx
  m.q_old = q
  m.qx_old = qx
  m.Ψ_old = Ψ


  # having solved for q, do the time step
  q = max(q, sqrt(eps()))
  ι = (q - 1) / κ
  Φ = log(q) / κ
  αE = χ * Ψ / x
  σx = x * (αE - 1) * σ / (1 - x * (αE - 1) * qx / q)
  σpE = pEx / pE * σx
  σpH = pHx / pH * σx
  σq = qx / q * σx
  σE = αE * (σ + σq)
  αH = (1 - αE * x) / (1 - x)
  σH = αH * (σ + σq)
  κE = γ * (σE - (1 / γ - 1) / (ψ - 1) * σpE)
  κH = γ * (σH - (1 / γ - 1) / (ψ - 1) * σpH)
  μx = x * (1 - x) * (κE * σE - κH * σH + 1 / pH - 1 / pE - (σE - σH) * (σ + σq))
  μx2 = x * ((a - ι) / q - 1 / pE) + σx * (κE - σ - σq) + x * (1 - χ) * (κE - κH) * (σ + σq)
  if (iter == 0) && (μx < 0)
    # upwinding
    iter += 1
    pEx, pHx = pEx_down, pHx_down
    @goto start
  end  
  #@assert abs((a - ι) / q - (x / pE + (1 - x) / pH)) < 1e-6
  μpE = pEx / pE * μx + 0.5 * pExx / pE * σx^2 
  μpH = pHx / pH * μx + 0.5 * pHxx / pH * σx^2 
  μq = qx / q * μx + 0.5 * qxx / q * σx^2

  r = clamp((a - ι) / q + Φ - δ + μq + σ * σq - (χ * κE + (1 - χ) * κH) * (σ + σq), -3.0, 3.0)
  # Market Pricing
  pEt = -  pE * (1 / pE - ψ * ρE + (ψ - 1) * (r + κE * σE) + μpE - (ψ - 1) * γ / 2 * σE^2 + (2 - ψ - γ) / (2 * (ψ - 1)) * σpE^2 + (1 - γ) * σpE * σE)
  pHt = -  pH * (1 / pH - ψ * ρH + (ψ - 1) * (r + κH * σH) + μpH - (ψ - 1) * γ / 2 * σH^2 + (2 - ψ - γ) / (2 * (ψ - 1)) * σpH^2 + (1 - γ) * σpH * σH)
  d = Ψ * a + (1 - Ψ) * alow - ι
  μR = r + (χ * κE + (1 - χ) * κH) * (σ + σq)
  (; pEt, pHt), (; x, r, μx, σx, Ψ, κ, κE, κH, σE, σH, σq, pEt, pHt, q, qx, qxx, μq, μx2, Φ, χ, d, μR, pExx, pHxx)
end

m = BrunnermeierSannikov()
xn = 200
stategrid =  OrderedDict(:x => range(0, 1.0, length = xn+2)[2:(end-1)])
Δx = step(stategrid[:x])
xmin = minimum(stategrid[:x])
yend = OrderedDict(:pE =>  9 .* ones(xn), :pH =>   10 .* ones(xn))
y, residual_norm, a = pdesolve((state, y) -> m(state, y, Δx, xmin), stategrid, yend; autodiff = :finite, maxΔ = 1000.0)