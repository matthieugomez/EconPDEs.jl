# # Achdou–Han–Lasry–Lions–Moll: consumption–saving with two income states
#
# A household saves in a riskless asset ``a`` at rate ``r`` and earns a labor income that
# switches between a low value ``y_l`` and a high value ``y_h`` as a two-state Poisson process
# (intensities ``\lambda_{lh}``, ``\lambda_{hl}``). There is **one** continuous state (assets
# ``a``) but **two coupled value functions** ``v_l(a)`` and ``v_h(a)``, one per income state,
# linked by the income transitions:
#
# ```math
# \rho\, v_j(a) = \max_c \frac{c^{1-\gamma}}{1-\gamma} + v_j'(a)\,(y_j + r a - c) + \lambda_{jk}\,\bigl(v_k(a)-v_j(a)\bigr), \qquad j \ne k.
# ```
#
# A **borrowing constraint** ``a \ge \underline a`` is enforced by preventing the asset drift
# from turning negative at the lower bound.

# ## The model

using EconPDEs, Plots

mutable struct AchdouHanLasryLionsMoll_TwoStatesModel
    yl::Float64
    yh::Float64
    λlh::Float64
    λhl::Float64
    r::Float64
    ρ::Float64
    γ::Float64
    amin::Float64
    amax::Float64
end

function AchdouHanLasryLionsMoll_TwoStatesModel(; yl = 0.5, yh = 1.5, λlh = 0.2, λhl = 0.2, r = 0.03, ρ = 0.04, γ = 2.0, amin = -yl / r, amax = 50.0)
    AchdouHanLasryLionsMoll_TwoStatesModel(yl, yh, λlh, λhl, r, ρ, γ, amin, amax)
end

# In each income state we upwind the asset drift on its sign, capping the implied consumption
# rather than flooring the marginal value (Newton may try negative marginal values). At the
# borrowing constraint the drift is set to zero. We save consumption `cl`, `ch` to plot.

function (m::AchdouHanLasryLionsMoll_TwoStatesModel)(state::NamedTuple, u::NamedTuple)
    (; yl, yh, λlh, λhl, r, ρ, γ, amin, amax) = m
    (; a) = state
    (; vl, vla_up, vla_down, vh, vha_up, vha_down) = u
    clmax = 100.0 * (yl + r * max(a, 0.0))
    chmax = 100.0 * (yh + r * max(a, 0.0))

    ## upwind the low-income value function
    cl_up = vla_up > 0 ? min(vla_up^(-1 / γ), clmax) : clmax
    μla_up = yl + r * a - cl_up
    if μla_up >= 0.0
        vla, cl, μla = vla_up, cl_up, μla_up
    else
        cl_down = vla_down > 0 ? min(vla_down^(-1 / γ), clmax) : clmax
        μla_down = yl + r * a - cl_down
        if μla_down <= 0.0 && a > amin
            vla, cl, μla = vla_down, cl_down, μla_down
        else
            cl = yl + r * a          # borrowing constraint binds: drift is zero
            μla = 0.0
            vla = cl^(-γ)
        end
    end
    vlt = -(cl^(1 - γ) / (1 - γ) + μla * vla + λlh * (vh - vl) - ρ * vl)

    ## upwind the high-income value function
    ch_up = vha_up > 0 ? min(vha_up^(-1 / γ), chmax) : chmax
    μha_up = yh + r * a - ch_up
    if μha_up >= 0.0
        vha, ch, μha = vha_up, ch_up, μha_up
    else
        ch_down = vha_down > 0 ? min(vha_down^(-1 / γ), chmax) : chmax
        μha_down = yh + r * a - ch_down
        if μha_down <= 0.0 && a > amin
            vha, ch, μha = vha_down, ch_down, μha_down
        else
            ch = yh + r * a
            μha = 0.0
            vha = ch^(-γ)
        end
    end
    vht = -(ch^(1 - γ) / (1 - γ) + μha * vha + λhl * (vl - vh) - ρ * vh)

    return (; vlt, vht), (; cl, ch, μla, μha)
end

# ## Solving it
#
# Build the asset grid (finer near the borrowing limit), start from an autarky-style guess,
# and solve.

m = AchdouHanLasryLionsMoll_TwoStatesModel()
m.amin += 0.001
stategrid = (; a = m.amin .+ range(0, (m.amax - m.amin)^(1 / 2), length = 200) .^ 2)
yend = (;
    vl = (m.ρ ./ m.γ .+ (1 .- 1 / m.γ) .* m.r)^(-m.γ) .* (stategrid[:a] .+ m.yl ./ m.r) .^ (1 - m.γ) ./ (1 - m.γ),
    vh = (m.ρ ./ m.γ .+ (1 .- m.γ) .* m.r)^(-m.γ) .* (stategrid[:a] .+ m.yh ./ m.r) .^ (1 - m.γ) ./ (1 - m.γ),
)
result = pdesolve(m, stategrid, yend)

# ## The solution
#
# The value function is higher in the high-income state (left). Consumption rises with wealth
# in both states and is higher when income is high (right). Near the borrowing limit the
# low-income household is forced to consume its income ``y_l + r a`` — the constraint binds and
# precautionary saving disappears.

as = stategrid[:a]
mask = as .<= 5.0
p1 = plot(as[mask], [result.zero[:vl][mask] result.zero[:vh][mask]]; label = ["low income" "high income"], xlabel = "assets a", ylabel = "value v(a)", legend = :bottomright)
p2 = plot(as[mask], [result.optional[:cl][mask] result.optional[:ch][mask]]; label = ["low income" "high income"], xlabel = "assets a", ylabel = "consumption c(a)", legend = :bottomright)
plot(p1, p2; layout = (1, 2), size = (800, 300))
