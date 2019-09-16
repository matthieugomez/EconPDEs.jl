using EconPDEs

Base.@kwdef  struct CampbellCochraneModel
    # consumption process parameters
    μ::Float64 = 0.0189
    σ::Float64 = 0.015

    # utility
    γ::Float64 = 2.0
    ρ::Float64 = 0.116

    # habit
    κs::Float64 = 0.138
    b::Float64 = 0.0
end
# I choose persistence so that monthly simulation of the model matches processes in CC (1999)
# ρ = 12 * (1 - 0.89^(1/12))
# κs = 12 * (1 - 0.87^(1/12))

function initialize_stategrid(m::CampbellCochraneModel; smin = -300.0, n = 1000)
    μ = m.μ ; σ = m.σ ; γ = m.γ ; ρ = m.ρ ; κs = m.κs ; b = m.b
    Sbar = σ * sqrt(γ / (κs - b / γ))
    sbar = log.(Sbar)
    smax =  sbar + 0.5 * (1 - Sbar^2)
    # corresponds to Grid 3 in Wachter (2005)
    shigh = log.(range(0.0, stop = exp(smax), length = div(n, 10)))
    slow = range(smin, stop = shigh[2], length = n - div(n, 10))
    OrderedDict(:s => vcat(slow[1:(end-1)], shigh[2:end]))
end

function initialize_y(m::CampbellCochraneModel, stategrid)
    OrderedDict(:p => ones(length(stategrid[:s])))
end
	
function (m::CampbellCochraneModel)(state::NamedTuple, y::NamedTuple)
    μ = m.μ ; σ = m.σ ; γ = m.γ ; ρ = m.ρ ; κs = m.κs ; b = m.b
    s = state.s
    p, ps, pss = y.p, y.ps, y.pss
    
    # drift and volatility of  s and p
    Sbar = σ * sqrt(γ / (κs - b / γ))
    sbar = log(Sbar)
    λ = 1 / Sbar * sqrt(1 - 2 * (s - sbar)) - 1
    μs = - κs * (s - sbar)
    σs = λ * σ
    σp = ps / p * σs
    μp = ps / p * μs + 0.5 * pss / p * σs^2

    # market price of risk κ
    κ = γ * (σ + σs)

    # risk free rate  r
    r = ρ + γ * μ - (γ * κs - b) / 2 + b * (sbar - s)

    # PDE
    pt = p * (1 / p + μ + μp + σp * σ - r - κ * (σ + σp))
    return (pt,), (μs, ), (p = p, κ = κ, λ = λ, r = r, σp = σp, μs = μs, σs = σs)
end


# Campbell Cochrane (1999)
m = CampbellCochraneModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
y, result, distance = pdesolve(m, stategrid, y0)



# Wachter (2005) calibration
# m = CampbellCochraneModel(μ = 0.022, σ = 0.0086, γ = 2.0, ρ = 0.073, κs = 0.116, b = 0.011 * 4)
# stategrid = initialize_stategrid(m)
# y0 = initialize_y(m, stategrid)
# y, result, distance = pdesolve(m, stategrid, y0)

