# # Tuckman–Vila (1992): a finite-horizon problem
#
# So far every model has been stationary. Tuckman–Vila is **time-dependent**: an arbitrageur
# faces holding costs over a finite horizon ``[0, T]``, so the value function ``F(z, \tau)``
# depends on calendar time. `pdesolve` handles this by taking a time grid `τs` as a fourth
# argument and solving **backward** from a terminal condition. The PDE function then takes a
# third argument `τ`:
#
# ```math
# \partial_\tau F + \max_i\Bigl\{ (\mu + \sigma^2 F_z)\, a i - \tfrac12 \sigma^2 (a i)^2 \Bigr\}
#   - \rho z F_z + \tfrac12 \sigma^2 (F_{zz} - F_z^2) = 0.
# ```

# ## The model

using EconPDEs, Distributions, Plots

Base.@kwdef struct TuckmanVilaModel
    c::Float64 = 0.06
    r::Float64 = 0.09
    ρ::Float64 = 5.42
    σ::Float64 = 26.72
    a::Float64 = 26.72
    T::Float64 = 100
end

function initialize_stategrid(m::TuckmanVilaModel; zn = 200)
    d = Normal(0, sqrt(m.σ^2 / (2 * m.ρ)))
    (; z = range(quantile(d, 0.00001), quantile(d, 0.99999), length = zn))
end

function initialize_y(m::TuckmanVilaModel, stategrid)
    zn = length(stategrid[:z])
    (; F = zeros(zn))
end

function (m::TuckmanVilaModel)(state::NamedTuple, u::NamedTuple, τ::Number)
    (; c, r, ρ, σ, a, T) = m
    (; z) = state
    (; F, Fz_up, Fz_down, Fzz) = u
    Fz = (z >= 0) ? Fz_up : Fz_down
    ϕ = z * (1 + z^2)^(-1 / 2)
    ϕz = (1 + z^2)^(-3 / 2)
    ϕzz = -3 * z * (1 + z^2)^(-5 / 2)
    sτ = c / r * (1 - exp(-r * (T - τ)))
    μL = ρ * z - 0.5 * σ^2 * ϕzz / ϕz + c / sτ * (ϕ - 1) / ϕz
    iL = 1 / a * (μL / σ^2 + Fz)
    μS = ρ * z - 0.5 * σ^2 * ϕzz / ϕz + c / sτ * (ϕ + 1) / ϕz
    iS = 1 / a * (μS / σ^2 + Fz)
    μ, i = 0.0, 0.0
    if iL > 0
        i, μ = iL, μL
    elseif iS < 0
        i, μ = iS, μS
    end
    Ft = -((μ + σ^2 * Fz) * a * i - 0.5 * σ^2 * (a * i)^2 - ρ * z * Fz + 0.5 * σ^2 * (Fzz - Fz^2))
    return (; Ft)
end

# ## Solving it
#
# `yend` is the terminal value at `τs[end]`; `pdesolve` marches backward over `τs`.

m = TuckmanVilaModel()
stategrid = initialize_stategrid(m)
yend = initialize_y(m, stategrid)
τs = range(0, m.T, length = 10)
result = pdesolve(m, stategrid, yend, τs)

# ## The solution
#
# For a time-dependent problem `result.zero[i]` is the solution at time `τs[i]` (so
# `result.zero[end]` is the terminal condition). Plotting a few time slices shows the value
# building up as the remaining horizon lengthens — from the flat terminal ``F = 0`` toward the
# long-horizon profile.

zs = stategrid[:z]
plt = plot(; xlabel = "state z", ylabel = "value F(z)", legend = :topright)
for i in (1, 4, 7, 10)
    plot!(plt, zs, result.zero[i][:F]; label = "τ = $(round(τs[i], digits = 0))")
end
plt
