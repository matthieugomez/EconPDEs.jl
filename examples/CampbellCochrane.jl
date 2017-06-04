mutable struct CampbellCochraneModel
    # consumption process parameters
    μ::Float64 
    σ::Float64

    # utility
    γ::Float64
    ρ::Float64

    # habit
    κs::Float64
    b::Float64
end

function CampbellCochraneModel(;μ = 0.0189, σ = 0.015, γ = 2.0, ρ = 0.116, κs = 0.138, b = 0.0)
    # I choose persistence so that monthly simulation of the model matches processes in CC (1999)
    # ρ = 12 * (1 - 0.89^(1/12))
    # κs = 12 * (1 - 0.87^(1/12))
    CampbellCochraneModel(μ, σ, γ, ρ, κs, b)
end

function initialize_state(m::CampbellCochraneModel; smin = -300.0, n = 1000)
    μ = m.μ ; σ = m.σ ; γ = m.γ ; ρ = m.ρ ; κs = m.κs ; b = m.b
    Sbar = σ * sqrt(γ / (κs - b / γ))
    sbar = log.(Sbar)
    smax =  sbar + 0.5 * (1 - Sbar^2)
    # corresponds to Grid 3 in Wachter (2005)
    shigh = log.(linspace(0.0, exp(smax), div(n, 10)))
    slow = linspace(smin, shigh[2], n - div(n, 10))
    OrderedDict(:s => vcat(slow[1:(end-1)], shigh[2:end]))
end

function initialize_y(m::CampbellCochraneModel, state)
    OrderedDict(:p => ones(length(state[:s])))
end
	
function (m::CampbellCochraneModel)(state, y)
    μ = m.μ ; σ = m.σ ; γ = m.γ ; ρ = m.ρ ; κs = m.κs ; b = m.b
    s = state.s
    p, ps, pss = y.p, y.ps, y.pss
    
    # drift and volatility of state variable s
    Sbar = σ * sqrt(γ / (κs - b / γ))
    sbar = log(Sbar)
    λ = 1 / Sbar * sqrt(1 - 2 * (s - sbar)) - 1
    μs = - κs * (s - sbar)
    σs = λ * σ

    # market price of risk κ
    κ = γ * (σ + σs)

    # risk free rate  r
    r = ρ + γ * μ - (γ * κs - b) / 2 + b * (sbar - s)

    # drift and volatility of p
    σp = ps / p * σs
    μp = ps / p * μs + 0.5 * pss / p * σs^2

    # PDE
    pt = p * (1 / p + μ + μp + σp * σ - r - κ * (σ + σp))
    return pt, μs, tuple(:p => p, :κ => κ, :λ => λ, :r => r, :σp => σp, :μs => μs, :σs => σs)
end


# # Habit Model
# ## Campbell Cochrane (1999)
# m = CampbellCochrane()
# state = initialize_state(m)
# y0 = initialize_y(m, state)
# result, distance = pdesolve(m, state, y0)

# ## Wachter (2005) calibration
# m = CampbellCochraneModel(μ = 0.022, σ = 0.0086, γ = 2.0, ρ = 0.073, κs = 0.116, b = 0.011 * 4)
# state = initialize_state(m)
# y0 = initialize_y(m, state)
# result, distance = pdesolve(m, state, y0)
