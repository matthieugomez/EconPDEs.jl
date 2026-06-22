using EconPDEs

# Deterministic neoclassical (Ramsey–Cass–Koopmans) growth model.
#
# A representative household chooses consumption to solve
#     max ∫₀^∞ e^{-ρt} u(c) dt    s.t.    k̇ = A k^α - δ k - c
# with CRRA utility u(c) = c^{1-γ} / (1 - γ).
#
# The value function v(k) solves the HJB equation
#     ρ v(k) = max_c u(c) + v'(k) (A k^α - δ k - c)
# with the first-order condition u'(c) = v'(k), i.e. c = v'(k)^{-1/γ}.
#
# Since the drift k̇ can be of either sign, the first derivative is upwinded:
# the forward difference is used where the implied drift is positive, the
# backward difference where it is negative, and the consumption that sets the
# drift to zero at the steady state otherwise (Achdou–Han–Lasry–Lions–Moll).

Base.@kwdef struct NeoclassicalGrowthModel
    A::Float64 = 0.5    # total factor productivity
    α::Float64 = 0.3    # capital share
    δ::Float64 = 0.05   # depreciation rate
    ρ::Float64 = 0.05   # discount rate
    γ::Float64 = 2.0    # relative risk aversion
end

function (m::NeoclassicalGrowthModel)(state::NamedTuple, value::NamedTuple)
    (; A, α, δ, ρ, γ) = m
    (; k) = state
    (; v, vk_up, vk_down) = value
    # upwind the first derivative on the sign of the drift
    cF = max(vk_up, eps())^(-1 / γ)
    μF = A * k^α - δ * k - cF
    cB = max(vk_down, eps())^(-1 / γ)
    μB = A * k^α - δ * k - cB
    if μF > 0
        c, vk, μk = cF, vk_up, μF
    elseif μB < 0
        c, vk, μk = cB, vk_down, μB
    else
        # at the steady state the household consumes its net output
        c, vk, μk = A * k^α - δ * k, vk_up, 0.0
    end
    vt = - (c^(1 - γ) / (1 - γ) + μk * vk - ρ * v)
    return (; vt)
end

m = NeoclassicalGrowthModel()
(; A, α, δ, ρ, γ) = m

# Closed-form steady-state capital: f'(k) = ρ + δ ⟹ α A k^{α-1} = ρ + δ
kstar = (α * A / (ρ + δ))^(1 / (1 - α))

stategrid = OrderedDict(:k => range(0.1 * kstar, 5 * kstar, length = 1000))
# initial guess: value of consuming gross output forever (monotone increasing in k,
# so that v'(k) > 0 and the implied consumption v'(k)^{-1/γ} stays well-defined)
yend = OrderedDict(:v => [(A * k^α)^(1 - γ) / (1 - γ) / ρ for k in stategrid[:k]])
result = pdesolve(m, stategrid, yend)
@assert result.residual_norm <= 1e-5
