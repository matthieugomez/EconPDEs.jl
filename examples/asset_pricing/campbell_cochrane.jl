# # Campbell–Cochrane (1999): habit formation
#
# Campbell and Cochrane explain the equity premium with an external habit. The single state
# is the log surplus consumption ratio ``s = \log\frac{C-X}{C}`` — how far consumption sits
# above the habit ``X`` — which mean-reverts around ``\bar s``. The price–consumption ratio
# ``p(s)`` solves a linear second-order ODE (a special case of the HJB machinery):
#
# ```math
# 0 = 1 + \bigl(\mu + \mu_p(s) + \sigma_p(s)\,\sigma - r(s) - \kappa(s)\,(\sigma+\sigma_p(s))\bigr)\, p(s),
# ```
#
# where the drift ``\mu_p``, volatility ``\sigma_p``, price of risk ``\kappa``, and riskless
# rate ``r`` are all functions of ``s`` through the habit.

# ## The model
#
# The parameters live in a `struct`:

using EconPDEs, Plots

Base.@kwdef struct CampbellCochraneModel
    μ::Float64 = 0.0189    # mean consumption growth
    σ::Float64 = 0.015     # volatility of consumption growth
    γ::Float64 = 2.0       # utility curvature (relative risk aversion)
    ρ::Float64 = 0.116     # discount rate
    κs::Float64 = 0.138    # mean-reversion speed of log surplus consumption ratio
    b::Float64 = 0.0       # sensitivity of riskless rate to surplus ratio (b=0 ⇒ constant r)
end

# ## The state space
#
# We build the grid and the initial guess first, because they fix the names used everywhere
# else. The grid is a `NamedTuple` whose key is the state variable (`s`, the log surplus
# consumption ratio); the guess is a `NamedTuple` whose key is the unknown function (`p`, the
# price–consumption ratio). These names reappear inside the equation below — e.g. `ps_up` will
# be the forward finite difference of `p` in `s`. The grid concentrates points near the
# reflecting upper bound ``s_{\max}``, where the surplus ratio is highest, and stretches far
# into the left tail.

function initialize_stategrid(m::CampbellCochraneModel; sn = 1000)
    (; μ, σ, γ, ρ, κs, b) = m
    Sbar = σ * sqrt(γ / (κs - b / γ))
    sbar = log(Sbar)
    smax = sbar + 0.5 * (1 - Sbar^2)
    shigh = log.(range(0.0, exp(smax), length = div(sn, 10)))
    slow = range(-300.0, shigh[2], length = sn - div(sn, 10))
    (; s = vcat(slow[1:(end-1)], shigh[2:end]))
end

m = CampbellCochraneModel()
stategrid = initialize_stategrid(m)
guess = (; p = ones(length(stategrid[:s])))

# ## The equation
#
# We now write the function encoding the HJB equation. Following the package convention, it
# takes the current `state` (a grid point) and `u` (each unknown together with its
# finite-difference derivatives there) and returns the time derivative of each unknown.

function (m::CampbellCochraneModel)(state::NamedTuple, u::NamedTuple)
    (; μ, σ, γ, ρ, κs, b) = m
    (; s) = state
    (; p, ps_up, ps_down, pss) = u
    Sbar = σ * sqrt(γ / (κs - b / γ))
    sbar = log(Sbar)
    λ = 1 / Sbar * sqrt(1 - 2 * (s - sbar)) - 1
    μs = -κs * (s - sbar)
    σs = λ * σ
    ## pricing uses the risk-adjusted drift of surplus consumption, not the physical one
    μs_effective = μs + σs * (σ - γ * (σ + σs))
    ps = (μs_effective >= 0) ? ps_up : ps_down
    σp = ps / p * σs
    μp = ps / p * μs + 0.5 * pss / p * σs^2
    κ = γ * (σ + σs)
    r = ρ + γ * μ - (γ * κs - b) / 2 + b * (sbar - s)
    pt = -p * (1 / p + μ + μp + σp * σ - r - κ * (σ + σp))
    return (; pt)
end

# With the equation, grid, and guess in hand, `pdesolve` solves the stationary system:

result = pdesolve(m, stategrid, guess)

# ## The solution
#
# Plotted against the surplus consumption ratio ``S = e^s``, the price–consumption ratio rises
# steeply in good times: as consumption pulls away from the habit, effective risk aversion
# falls, discount rates drop, and valuations climb.

Sbar = m.σ * sqrt(m.γ / (m.κs - m.b / m.γ))
sbar = log(Sbar)
ss = stategrid[:s]
mask = ss .>= (sbar - 4)
plot(exp.(ss[mask]), result.solution.p[mask]; xlabel = "surplus consumption ratio S = exp(s)", ylabel = "price–consumption ratio p", legend = false)
