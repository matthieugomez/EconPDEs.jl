# # Leland (1994): optimal default
#
# Leland's model of a levered firm is a classic optimal-stopping problem. Cash flow (EBIT)
# ``\delta`` follows a geometric Brownian motion ``d\delta = \mu\delta\,dt + \sigma\delta\,dW``.
# The firm has issued perpetual debt with coupon ``C`` and faces a tax rate ``\tau``; equity
# holders receive the after-tax flow ``\delta-(1-\tau)C`` but may **default** (walk away)
# whenever they choose. Equity value ``E(\delta)`` is thus the value of the option to keep the
# firm alive, and solves the HJB variational inequality
#
# ```math
# \min\Bigl\{\, r E - \bigl(\delta-(1-\tau)C\bigr) - \mu\delta\,E'(\delta) - \tfrac12\sigma^2\delta^2 E''(\delta),\;\; E(\delta)\,\Bigr\} = 0,
# ```
#
# with default payoff ``0`` (limited liability). This delivers the value-matching
# ``E(\delta^*)=0`` and smooth-pasting ``E'(\delta^*)=0`` conditions at the endogenous default
# threshold ``\delta^*``.

# ## The model
#
# The parameters live in a `struct`:

using EconPDEs, Plots

Base.@kwdef struct LelandModel
    r::Float64 = 0.02    # risk-free rate
    μ::Float64 = 0.0     # drift of cash flow (EBIT growth)
    σ::Float64 = 0.1     # volatility of cash flow
    C::Float64 = 0.05    # coupon rate
    τ::Float64 = 0.2     # tax rate
end

# ## The state space
#
# We build the grid and the initial guess first, because they fix the names used everywhere
# else. The grid is a `NamedTuple` whose key is the state variable (`δ`, the cash flow); the
# guess is a `NamedTuple` whose key is the unknown function (`E`, equity value), holding one
# starting value at each grid point. These names reappear inside the equation below — e.g.
# `Eδ_up` will be the forward finite difference of `E` in `δ`. We start the guess from the value
# of the cash-flow claim net of the after-tax coupon, floored at zero.

m = LelandModel()
stategrid = (; δ = range(0, 0.2, step = 0.005))
yend = (; E = max.(stategrid[:δ] ./ (m.r - m.μ) .- (1 .- m.τ) .* m.C ./ m.r, 0.0))

# ## The equation
#
# We now write the function encoding the HJB equation. Following the package convention, it
# takes the current `state` (a grid point) and `u` (each unknown together with its
# finite-difference derivatives there) and returns the time derivative of each unknown.

function (m::LelandModel)(state::NamedTuple, u::NamedTuple)
    (; r, μ, σ, C, τ) = m
    (; δ) = state
    (; E, Eδ_up, Eδ_down, Eδδ) = u
    μδ, σδ = μ * δ, σ * δ
    Eδ = (μδ >= 0) ? Eδ_up : Eδ_down
    f = δ - (1 - τ) * C
    Et = -(f + Eδ * μδ + 0.5 * Eδδ * σδ^2 - r * E)
    return (; Et)
end

# The stopping (default) region is imposed with the lower bound `y̲ = 0`: equity can never be
# worth less than zero. We also pin the slope at the top boundary to that of a debt-free claim,
# ``1/(r-\mu)``. With these in hand, `pdesolve` solves the variational inequality:

y̲ = zeros(length(stategrid[:δ]))
bc = (; Eδ = (0.0, 1 / (m.r - m.μ)))
result = pdesolve(m, stategrid, yend; y̲ = y̲, bc = bc, reformulation = :smooth)

# ## The solution
#
# Equity is worthless below an endogenous default threshold ``\delta^*`` (red line) and rises
# smoothly above it. Shareholders keep servicing the debt only while cash flows are high
# enough to justify the coupon; below ``\delta^*`` they optimally default.

δs = stategrid[:δ]
E = result.zero[:E]
δstar = δs[findfirst(>(1e-8), E)]
plot(δs, E; xlabel = "cash flow δ", ylabel = "equity value E(δ)", legend = false)
vline!([δstar]; color = :red, linestyle = :dot)
