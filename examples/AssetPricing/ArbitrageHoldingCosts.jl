using EconPDEs, Distributions

# Arbitrage With Holding Costs: A Utility-Based Approach
# Author(s): Bruce Tuckman and Jean-Luc Vila

struct ArbitrageHoldingCosts
    c::Float64
    r::Float64
    ρ::Float64 
    σ::Float64 
    a::Float64 
    T::Float64
end

function ArbitrageHoldingCosts(;c = 0.06, r = 0.09,  ρ = 5.42, σ = 0.88, a = 26.72, T = 100)
    ArbitrageHoldingCosts(c, r, ρ, σ, a, T)
end

function initialize_stategrid(m::ArbitrageHoldingCosts; n = 200)
    d = Normal(0, sqrt(m.σ^2 / (2 * m.ρ)))
    zmin = quantile(d, 0.00001)
    zmax = quantile(d, 0.99999)
    OrderedDict(:z => range(zmin, stop = zmax, length = n))
end

function initialize_y(m::ArbitrageHoldingCosts, stategrid::OrderedDict)
    OrderedDict(:F => zeros(length(stategrid[:z])))
end

function (m::ArbitrageHoldingCosts)(state::NamedTuple, y::NamedTuple, τ::Number)
    (; c, r, ρ, σ, a, T) = m
    (; z) = state
    (; F, Fz, Fzz) = y
    ϕ = z * (1 + z^2)^(-1/2)
    ϕz = (1 + z^2)^(-3/2)
    ϕzz = - 3 * z * (1 + z^2)^(-5/2)
    sτ = c / r * (1 - exp(- r * (T - τ)))
    μL = ρ * z - 0.5 * σ^2 * ϕzz / ϕz + c / sτ * (ϕ - 1) / ϕz
    iL = 1 / a * (μL / σ^2 + Fz)
    μS = ρ * z - 0.5 * σ^2 * ϕzz / ϕz + c / sτ * (ϕ + 1) / ϕz
    iS = 1 / a * (μS / σ^2 + Fz)
    μ, i  = 0.0, 0.0
    if iL > 0
        i = iL
        μ = μL
    elseif iS < 0
        i = iS
        μ = μS
    end
    if z >= 0
        i_myopic = max(0, i - Fz / a)
    else
        i_myopic = min(0, i - Fz / a)
    end
    # otherwise i = 0.0 and value of μ does not matter
    Ft = (μ + σ^2 * Fz) * a * i - 0.5 * σ^2 * (a * i)^2 - ρ * z * Fz + 0.5 * σ^2 * (Fzz - Fz^2) 
    return (Ft,), (-ρ * z,)
end

m = ArbitrageHoldingCosts()
τs = range(m.T, stop = 0, length = 100)
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
y, result, distance = pdesolve(m, stategrid, y0, τs)


### reproduce Fig 2
#d = Normal(0, sqrt(m.σ^2 / (2 * m.ρ)))
#zmin = quantile(d, 0.025)
#zmax = quantile(d, 0.975)
#idx = (state[:z] .>= zmin) .& (state[:z] .<= zmax)
#
#using Plots
#plot(result[:x][idx, 20], [result[:I_myopic][idx, 2] result[:I][idx, 2]], label = ["myopic" "all"])
#plot(result[:x][idx, 20], result[:I][idx, 2])
#plot!(result[:x][idx, 60], result[:I][idx, 60])
#plot!(result[:x][idx, 100], result[:I][idx, 100])
#

