# # Gârleanu–Panageas with long-run risk: alternating values and prices
#
# This example solves a heterogeneous-agent Gârleanu–Panageas economy in which the
# aggregate endowment process follows long-run-risk dynamics. There are three state
# variables: `x`, the type-A consumption share; `μ`, expected endowment growth; and `v`,
# the conditional variance of endowment growth. The extra long-run-risk states make
# equilibrium prices and prices of risk vary over the state grid rather than reduce to
# constants.
#
# To my knowledge, this three-state GP long-run-risk economy has not previously been
# presented as a worked numerical solution. To help with convergence, the example below
# presents a helpful solution method: formulate the model in EconPDEs by treating the
# interest rate and prices of risk as state-dependent unknown functions, and alternate
# between valuation equations and market-clearing price updates.
#
# This example illustrates a useful way to write general-equilibrium asset-pricing
# systems. Rather than first solving the equilibrium pricing equations for the interest
# rate and prices of risk as functions of the value functions, we keep those prices as
# unknown functions on the grid. The unknown is
#
#     (pA, pB, ϕ1, ϕ2, r, κC, κM, κV)
#
# where the first four entries are valuation functions and the last four are
# state-dependent equilibrium price functions. Their equations are algebraic
# restrictions: the supplied `r`, `κC`, `κM`, and `κV` must equal the values implied by
# market clearing at that state.
#
# This formulation is often better conditioned than eliminating prices. In a simple
# one-dimensional model, substituting the equilibrium `r(V)` and `κ(V)` formulas directly
# into the HJB can simplify the system. In higher-dimensional models with several prices
# of risk, the same substitution makes the PDE operator much more nonlinear: prices affect
# consumption volatilities, the wealth-share drift, upwind directions, and cross
# derivatives. Keeping prices as algebraic unknowns makes the system larger but cleaner:
# optimization equations and market-clearing equations are solved simultaneously.
#
# The same trick can be useful in other asset-pricing models whenever endogenous prices,
# risk premia, or volatilities feed back strongly into the law of motion. It is not needed
# in every model: when the equilibrium prices have simple closed forms and substituting
# them leaves a well-conditioned HJB, the shorter reduced-form PDE is preferable.
# The computation alternates between fixed-price valuation equations and price updates
# from market clearing.

using EconPDEs
using Distributions
using Logging
using Plots
using Printf

struct GarleanuPanageasLongRunRiskModel
    γA::Float64
    ψA::Float64
    γB::Float64
    ψB::Float64
    ρ::Float64
    δ::Float64
    νA::Float64
    μbar::Float64
    vbar::Float64
    κμ::Float64
    νμ::Float64
    κv::Float64
    νv::Float64
    B1::Float64
    δ1::Float64
    B2::Float64
    δ2::Float64
    ω::Float64
end

function GarleanuPanageasLongRunRiskModel(;
        γA = 1.5, ψA = 0.7, γB = 10.0, ψB = 0.05, ρ = 0.001,
        δ = 0.02, νA = 0.01, μbar = 0.018, vbar = 0.00073,
        κμ = 0.252, νμ = 0.528, κv = 0.156, νv = 0.00354,
        B1 = 30.72, δ1 = 0.0525, B2 = -30.29, δ2 = 0.0611,
        ω = 0.92)

    scale = δ / (δ + δ1) * B1 + δ / (δ + δ2) * B2
    return GarleanuPanageasLongRunRiskModel(γA, ψA, γB, ψB, ρ, δ,
        νA, μbar, vbar, κμ, νμ, κv, νv, B1 / scale, δ1, B2 / scale, δ2, ω)
end

function (model::GarleanuPanageasLongRunRiskModel)(state::NamedTuple, u::NamedTuple)
    (; γA, ψA, γB, ψB, ρ, δ, νA, μbar, vbar, κμ, νμ, κv, νv, B1, δ1, B2, δ2, ω) = model
    (; x, μ, v) = state

    pA, pB, ϕ1, ϕ2 = u.pA, u.pB, u.ϕ1, u.ϕ2
    r, κC, κM, κV = u.r, u.κC, u.κM, u.κV

    σc = sqrt(v)
    μμ = κμ * (μbar - μ)
    σμ = νμ * sqrt(v)
    μv = κv * (vbar - v)
    σv = νv * sqrt(v)

    pAμ = (μμ >= 0) ? u.pAμ_up : u.pAμ_down
    pBμ = (μμ >= 0) ? u.pBμ_up : u.pBμ_down
    ϕ1μ = (μμ >= 0) ? u.ϕ1μ_up : u.ϕ1μ_down
    ϕ2μ = (μμ >= 0) ? u.ϕ2μ_up : u.ϕ2μ_down
    pAv = (μv >= 0) ? u.pAv_up : u.pAv_down
    pBv = (μv >= 0) ? u.pBv_up : u.pBv_down
    ϕ1v = (μv >= 0) ? u.ϕ1v_up : u.ϕ1v_down
    ϕ2v = (μv >= 0) ? u.ϕ2v_up : u.ϕ2v_down

    pAx, pBx, ϕ1x, ϕ2x = u.pAx_up, u.pBx_up, u.ϕ1x_up, u.ϕ2x_up
    iter = 0
    μx = 0.0
    σxC = 0.0
    σxM = 0.0
    σxV = 0.0
    σCA_C = 0.0
    σCA_M = 0.0
    σCA_V = 0.0
    σCB_C = 0.0
    σCB_M = 0.0
    σCB_V = 0.0
    κ2 = 0.0
    mcA = 0.0
    mcB = 0.0
    @label start

    Γ = 1 / (x / γA + (1 - x) / γB)
    aA = (1 - γA * ψA) / (γA * (ψA - 1))
    aB = (1 - γB * ψB) / (γB * (ψB - 1))

    pAx_pA = pAx / pA
    pBx_pB = pBx / pB
    pAμ_pA = pAμ / pA
    pBμ_pB = pBμ / pB
    pAv_pA = pAv / pA
    pBv_pB = pBv / pB

    denom = 1 - x * aA * pAx_pA
    σxC = x * (κC / γA - σc) / denom
    σxM = x * (κM / γA + aA * pAμ_pA * σμ) / denom
    σxV = x * (κV / γA + aA * pAv_pA * σv) / denom
    σx2 = σxC^2 + σxM^2 + σxV^2

    σpA_C = pAx_pA * σxC
    σpA_M = pAx_pA * σxM + pAμ_pA * σμ
    σpA_V = pAx_pA * σxV + pAv_pA * σv
    σpB_C = pBx_pB * σxC
    σpB_M = pBx_pB * σxM + pBμ_pB * σμ
    σpB_V = pBx_pB * σxV + pBv_pB * σv

    σCA_C = κC / γA + aA * σpA_C
    σCA_M = κM / γA + aA * σpA_M
    σCA_V = κV / γA + aA * σpA_V
    σCB_C = κC / γB + aB * σpB_C
    σCB_M = κM / γB + aB * σpB_M
    σCB_V = κV / γB + aB * σpB_V

    κ2 = κC^2 + κM^2 + κV^2
    σpA2 = σpA_C^2 + σpA_M^2 + σpA_V^2
    σpB2 = σpB_C^2 + σpB_M^2 + σpB_V^2
    κσpA = κC * σpA_C + κM * σpA_M + κV * σpA_V
    κσpB = κC * σpB_C + κM * σpB_M + κV * σpB_V
    mcA = (1 + ψA) / (2 * γA) * κ2 + aA * κσpA - 0.5 * aA * σpA2
    mcB = (1 + ψB) / (2 * γB) * κ2 + aB * κσpB - 0.5 * aB * σpB2

    human_capital = ω * (B1 * ϕ1 + B2 * ϕ2)
    μCA = ψA * (r - ρ) + mcA
    μCB = ψB * (r - ρ) + mcB
    μx = x * (μCA - μ) + δ * (νA / pA * human_capital - x) - σc * σxC

    if (iter == 0) && (μx <= 0)
        iter += 1
        pAx, pBx, ϕ1x, ϕ2x = u.pAx_down, u.pBx_down, u.ϕ1x_down, u.ϕ2x_down
        @goto start
    end

    pAxμ = (σxM * σμ >= 0) ? u.pAxμ_up : u.pAxμ_down
    pBxμ = (σxM * σμ >= 0) ? u.pBxμ_up : u.pBxμ_down
    ϕ1xμ = (σxM * σμ >= 0) ? u.ϕ1xμ_up : u.ϕ1xμ_down
    ϕ2xμ = (σxM * σμ >= 0) ? u.ϕ2xμ_up : u.ϕ2xμ_down
    pAxv = (σxV * σv >= 0) ? u.pAxv_up : u.pAxv_down
    pBxv = (σxV * σv >= 0) ? u.pBxv_up : u.pBxv_down
    ϕ1xv = (σxV * σv >= 0) ? u.ϕ1xv_up : u.ϕ1xv_down
    ϕ2xv = (σxV * σv >= 0) ? u.ϕ2xv_up : u.ϕ2xv_down

    ϕ1x_ϕ1 = ϕ1x / ϕ1
    ϕ2x_ϕ2 = ϕ2x / ϕ2
    ϕ1μ_ϕ1 = ϕ1μ / ϕ1
    ϕ2μ_ϕ2 = ϕ2μ / ϕ2
    ϕ1v_ϕ1 = ϕ1v / ϕ1
    ϕ2v_ϕ2 = ϕ2v / ϕ2

    σϕ1_C = ϕ1x_ϕ1 * σxC
    σϕ1_M = ϕ1x_ϕ1 * σxM + ϕ1μ_ϕ1 * σμ
    σϕ1_V = ϕ1x_ϕ1 * σxV + ϕ1v_ϕ1 * σv
    σϕ2_C = ϕ2x_ϕ2 * σxC
    σϕ2_M = ϕ2x_ϕ2 * σxM + ϕ2μ_ϕ2 * σμ
    σϕ2_V = ϕ2x_ϕ2 * σxV + ϕ2v_ϕ2 * σv

    μpA = pAx_pA * μx + pAμ_pA * μμ + pAv_pA * μv +
          0.5 * u.pAxx / pA * σx2 + 0.5 * u.pAμμ / pA * σμ^2 +
          0.5 * u.pAvv / pA * σv^2 + pAxμ / pA * σxM * σμ +
          pAxv / pA * σxV * σv
    μpB = pBx_pB * μx + pBμ_pB * μμ + pBv_pB * μv +
          0.5 * u.pBxx / pB * σx2 + 0.5 * u.pBμμ / pB * σμ^2 +
          0.5 * u.pBvv / pB * σv^2 + pBxμ / pB * σxM * σμ +
          pBxv / pB * σxV * σv
    μϕ1 = ϕ1x_ϕ1 * μx + ϕ1μ_ϕ1 * μμ + ϕ1v_ϕ1 * μv +
          0.5 * u.ϕ1xx / ϕ1 * σx2 + 0.5 * u.ϕ1μμ / ϕ1 * σμ^2 +
          0.5 * u.ϕ1vv / ϕ1 * σv^2 + ϕ1xμ / ϕ1 * σxM * σμ +
          ϕ1xv / ϕ1 * σxV * σv
    μϕ2 = ϕ2x_ϕ2 * μx + ϕ2μ_ϕ2 * μμ + ϕ2v_ϕ2 * μv +
          0.5 * u.ϕ2xx / ϕ2 * σx2 + 0.5 * u.ϕ2μμ / ϕ2 * σμ^2 +
          0.5 * u.ϕ2vv / ϕ2 * σv^2 + ϕ2xμ / ϕ2 * σxM * σμ +
          ϕ2xv / ϕ2 * σxV * σv

    σCAσpA = σCA_C * σpA_C + σCA_M * σpA_M + σCA_V * σpA_V
    σCBσpB = σCB_C * σpB_C + σCB_M * σpB_M + σCB_V * σpB_V
    κσpA_CA = κC * (σpA_C + σCA_C) + κM * (σpA_M + σCA_M) +
              κV * (σpA_V + σCA_V)
    κσpB_CB = κC * (σpB_C + σCB_C) + κM * (σpB_M + σCB_M) +
              κV * (σpB_V + σCB_V)

    pAt = -pA * (1 / pA + μCA + μpA + σCAσpA - r - δ - κσpA_CA)
    pBt = -pB * (1 / pB + μCB + μpB + σCBσpB - r - δ - κσpB_CB)

    κσϕ1_Y = κC * (σϕ1_C + σc) + κM * σϕ1_M + κV * σϕ1_V
    κσϕ2_Y = κC * (σϕ2_C + σc) + κM * σϕ2_M + κV * σϕ2_V
    ϕ1t = -ϕ1 * (1 / ϕ1 + μ - δ - δ1 + μϕ1 + σc * σϕ1_C - r - κσϕ1_Y)
    ϕ2t = -ϕ2 * (1 / ϕ2 + μ - δ - δ2 + μϕ2 + σc * σϕ2_C - r - κσϕ2_Y)

    κC_implied = Γ * (σc - x * aA * σpA_C - (1 - x) * aB * σpB_C)
    κM_implied = Γ * (0.0 - x * aA * σpA_M - (1 - x) * aB * σpB_M)
    κV_implied = Γ * (0.0 - x * aA * σpA_V - (1 - x) * aB * σpB_V)
    r_implied = ρ + (μ - x * mcA - (1 - x) * mcB -
                δ * ((νA / pA + (1 - νA) / pB) * human_capital - 1)) /
                (ψA * x + ψB * (1 - x))

    rt = r_implied - r
    κCt = κC_implied - κC
    κMt = κM_implied - κM
    κVt = κV_implied - κV

    return (; pAt, pBt, ϕ1t, ϕ2t, rt, κCt, κMt, κVt),
           (; κ = sqrt(κ2), r_implied, κC_implied, κM_implied, κV_implied,
              r_gap = rt, κC_gap = κCt, κM_gap = κMt, κV_gap = κVt,
              pA_residual = pAt, pB_residual = pBt,
              ϕ1_residual = ϕ1t, ϕ2_residual = ϕ2t,
              μx, σx = sqrt(σx2), σxC, σxM, σxV)
end




# Define the calibration and state grid.
m = GarleanuPanageasLongRunRiskModel()
xn, μn, vn = 10, 5, 5
μsd = sqrt(m.νμ^2 * m.vbar / (2 * m.κμ))
μdistribution = Normal(m.μbar, μsd)
vdistribution = Gamma(2 * m.κv * m.vbar / m.νv^2, m.νv^2 / (2 * m.κv))
xs = collect(range(0.0, 1.0, length = xn).^2)
μs = collect(range(quantile(μdistribution, 0.10),
                   quantile(μdistribution, 0.90), length = μn))
vs = collect(range(quantile(vdistribution, 0.10),
                   quantile(vdistribution, 0.90), length = vn))
stategrid = (; x = xs, μ = μs, v = vs)

# Define the initial guess for value functions and equilibrium price functions.
gridsize = map(length, values(stategrid))
pA = ones(gridsize...)
pB = ones(gridsize...)
ϕ1 = ones(gridsize...)
ϕ2 = ones(gridsize...)
r = zeros(gridsize...)
κC = zeros(gridsize...)
κM = zeros(gridsize...)
κV = zeros(gridsize...)
human_capital = m.ω * (m.B1 + m.B2)
for (ix, x) in pairs(stategrid.x), (iμ, μ) in pairs(stategrid.μ), (iv, v) in pairs(stategrid.v)
    σc = sqrt(v)
    Γ = 1 / (x / m.γA + (1 - x) / m.γB)
    κC[ix, iμ, iv] = Γ * σc
    κ2 = κC[ix, iμ, iv]^2
    mcA = (1 + m.ψA) / (2 * m.γA) * κ2
    mcB = (1 + m.ψB) / (2 * m.γB) * κ2
    r[ix, iμ, iv] = m.ρ + (μ - x * mcA - (1 - x) * mcB -
                      m.δ * (human_capital - 1)) / (m.ψA * x + m.ψB * (1 - x))
end

guess = (; pA, pB, ϕ1, ϕ2, r, κC, κM, κV)

# Alternate between prices and value functions. Holding `r`, `κC`, `κM`, and `κV` fixed,
# solve the valuation block `(pA, pB, ϕ1, ϕ2)`. Then update prices toward the values
# implied by market clearing.

struct GarleanuPanageasLongRunRiskFrozenPriceModel
    model::GarleanuPanageasLongRunRiskModel
    stategrid::NamedTuple
    r::Array{Float64, 3}
    κC::Array{Float64, 3}
    κM::Array{Float64, 3}
    κV::Array{Float64, 3}
end

function (model::GarleanuPanageasLongRunRiskFrozenPriceModel)(
        state::NamedTuple, u::NamedTuple)
    ix = searchsortedfirst(model.stategrid.x, state.x)
    iμ = searchsortedfirst(model.stategrid.μ, state.μ)
    iv = searchsortedfirst(model.stategrid.v, state.v)
    u_with_prices = merge(u, (; r = model.r[ix, iμ, iv],
                              κC = model.κC[ix, iμ, iv],
                              κM = model.κM[ix, iμ, iv],
                              κV = model.κV[ix, iμ, iv]))
    out, saved = model.model(state, u_with_prices)
    return (; pAt = out.pAt, pBt = out.pBt, ϕ1t = out.ϕ1t, ϕ2t = out.ϕ2t), saved
end

damping = 0.5
warmup_maxiters = 20
warmup_inner_maxiters = 10
warmups_between_global_attempts = 10
reuse_failed_global_tol = 1e-4

initial_result = with_logger(NullLogger()) do
    pdesolve(m, stategrid, guess; alg = NonlinearSolve.TrustRegion(), Δ = Inf,
             verbose = false)
end

@printf "%4s %-8s %12s %12s %12s %12s %12s %12s\n" "iter" "step" "r_gap" "κC_gap" "κM_gap" "κV_gap" "value_change" "residual"
@printf "%4d %-8s %12s %12s %12s %12s %12s %12.3e\n" 0 "global" "-" "-" "-" "-" "-" initial_result.residual_norm

if initial_result.converged
    new_guess = initial_result.solution
    result = initial_result
else
    new_guess, result = let
    if initial_result.residual_norm <= reuse_failed_global_tol
        value_guess = (; pA = copy(initial_result.solution.pA),
                       pB = copy(initial_result.solution.pB),
                       ϕ1 = copy(initial_result.solution.ϕ1),
                       ϕ2 = copy(initial_result.solution.ϕ2))
        r = copy(initial_result.solution.r)
        κC = copy(initial_result.solution.κC)
        κM = copy(initial_result.solution.κM)
        κV = copy(initial_result.solution.κV)
        candidate_guess = initial_result.solution
    else
        value_guess = (; pA = copy(guess.pA), pB = copy(guess.pB),
                       ϕ1 = copy(guess.ϕ1), ϕ2 = copy(guess.ϕ2))
        r = copy(guess.r)
        κC = copy(guess.κC)
        κM = copy(guess.κM)
        κV = copy(guess.κV)
        candidate_guess = guess
    end
    price_iteration = 0
    trial_result = nothing
    warmups_since_global_attempt = 0

    while trial_result === nothing
        price_iteration += 1

        frozen_model = GarleanuPanageasLongRunRiskFrozenPriceModel(
            m, stategrid, r, κC, κM, κV)
        frozen_result = with_logger(NullLogger()) do
            pdesolve(frozen_model, stategrid, value_guess;
                     is_algebraic = (; pA = false, pB = false,
                                     ϕ1 = false, ϕ2 = false),
                     alg = NonlinearSolve.TrustRegion(),
                     maxiters = warmup_maxiters,
                     inner_maxiters = warmup_inner_maxiters, verbose = false)
        end

        r_gap = maximum(abs, frozen_result.saved.r_gap)
        κC_gap = maximum(abs, frozen_result.saved.κC_gap)
        κM_gap = maximum(abs, frozen_result.saved.κM_gap)
        κV_gap = maximum(abs, frozen_result.saved.κV_gap)
        pA_change = maximum(abs, frozen_result.solution.pA .- value_guess.pA) /
                    max(1.0, maximum(abs, value_guess.pA))
        pB_change = maximum(abs, frozen_result.solution.pB .- value_guess.pB) /
                    max(1.0, maximum(abs, value_guess.pB))
        ϕ1_change = maximum(abs, frozen_result.solution.ϕ1 .- value_guess.ϕ1) /
                    max(1.0, maximum(abs, value_guess.ϕ1))
        ϕ2_change = maximum(abs, frozen_result.solution.ϕ2 .- value_guess.ϕ2) /
                    max(1.0, maximum(abs, value_guess.ϕ2))
        warmup_change = max(pA_change, pB_change, ϕ1_change, ϕ2_change)

        value_guess = frozen_result.solution
        r .= (1 - damping) .* r .+ damping .* frozen_result.saved.r_implied
        κC .= (1 - damping) .* κC .+ damping .* frozen_result.saved.κC_implied
        κM .= (1 - damping) .* κM .+ damping .* frozen_result.saved.κM_implied
        κV .= (1 - damping) .* κV .+ damping .* frozen_result.saved.κV_implied
        warmups_since_global_attempt += 1
        candidate_guess = (; value_guess.pA, value_guess.pB, value_guess.ϕ1,
                           value_guess.ϕ2, r, κC, κM, κV)
        @printf "%4d %-8s %12.3e %12.3e %12.3e %12.3e %12.3e %12s\n" price_iteration "values" r_gap κC_gap κM_gap κV_gap warmup_change "-"

        if warmups_since_global_attempt >= warmups_between_global_attempts
            trial_result = with_logger(NullLogger()) do
                pdesolve(m, stategrid, candidate_guess;
                         alg = NonlinearSolve.TrustRegion(), Δ = Inf,
                         verbose = false)
            end
            @printf "%4d %-8s %12s %12s %12s %12s %12s %12.3e\n" price_iteration "global" "-" "-" "-" "-" "-" trial_result.residual_norm
            warmups_since_global_attempt = 0
        end

        if trial_result !== nothing
            if trial_result.converged
                candidate_guess = trial_result.solution
                break
            end
            if trial_result.residual_norm <= reuse_failed_global_tol
                value_guess = (; pA = trial_result.solution.pA,
                               pB = trial_result.solution.pB,
                               ϕ1 = trial_result.solution.ϕ1,
                               ϕ2 = trial_result.solution.ϕ2)
                r .= (1 - damping) .* r .+ damping .* trial_result.solution.r
                κC .= (1 - damping) .* κC .+ damping .* trial_result.solution.κC
                κM .= (1 - damping) .* κM .+ damping .* trial_result.solution.κM
                κV .= (1 - damping) .* κV .+ damping .* trial_result.solution.κV
            end
            trial_result = nothing
            continue
        end
    end

    (; value_guess.pA, value_guess.pB, value_guess.ϕ1, value_guess.ϕ2,
     r, κC, κM, κV), trial_result
end
end

# ## The solution
#
# The model has three state variables, so we inspect one-dimensional slices of the
# equilibrium price schedules. First vary expected growth ``\mu`` while holding variance
# fixed at the middle grid point. Then vary variance ``v`` while holding expected growth
# fixed at the middle grid point. The total price of risk is
# ``\kappa = \sqrt{\kappa_C^2 + \kappa_\mu^2 + \kappa_v^2}``.

xgrid = stategrid.x
μ_mid_index = cld(length(stategrid.μ), 2)
v_mid_index = cld(length(stategrid.v), 2)
μ_indices = unique([1, μ_mid_index, length(stategrid.μ)])
v_indices = unique([1, v_mid_index, length(stategrid.v)])

r_by_μ = plot(; xlabel = "type-A consumption share x", ylabel = "interest rate r")
κ_by_μ = plot(; xlabel = "type-A consumption share x", ylabel = "total price of risk κ")
for iμ in μ_indices
    label = "μ = $(round(stategrid.μ[iμ]; digits = 3))"
    plot!(r_by_μ, xgrid, result.solution.r[:, iμ, v_mid_index]; label)
    plot!(κ_by_μ, xgrid, result.saved.κ[:, iμ, v_mid_index]; label)
end
plot(r_by_μ, κ_by_μ; layout = (1, 2), size = (820, 300))

r_by_v = plot(; xlabel = "type-A consumption share x", ylabel = "interest rate r")
κ_by_v = plot(; xlabel = "type-A consumption share x", ylabel = "total price of risk κ")
for iv in v_indices
    label = "v = $(round(stategrid.v[iv]; sigdigits = 3))"
    plot!(r_by_v, xgrid, result.solution.r[:, μ_mid_index, iv]; label)
    plot!(κ_by_v, xgrid, result.saved.κ[:, μ_mid_index, iv]; label)
end
plot(r_by_v, κ_by_v; layout = (1, 2), size = (820, 300))
