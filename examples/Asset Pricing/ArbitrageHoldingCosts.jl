using Distributions

# Arbitrage With Holding Costs: A Utility-Based Approach
# Author(s): Bruce Tuckman and Jean-Luc Vila

struct ArbitrageHoldingCosts
    c::Float64
    r::Float64
    ρ::Float64 
    σ::Float64 
    a::Float64 
end

function ArbitrageHoldingCosts(;c = 0.06, r = 0.09,  ρ = 5.42, σ = 0.88, a = 26.72)
    ArbitrageHoldingCosts(c, r, ρ, σ, a)
end


function initialize_state(m::ArbitrageHoldingCosts; n = 200)
    d = Normal(0, sqrt(m.σ^2 / (2 * m.ρ)))
    zmin = quantile(d, 0.00001)
    zmax = quantile(d, 0.99999)
    OrderedDict(:z => collect(range(zmin, stop = zmax, length = n)))
end

function initialize_y(m::ArbitrageHoldingCosts, state)
    OrderedDict(:F => zeros(length(state[:z])))
end

function (m::ArbitrageHoldingCosts)(state, y, τ)
    c = m.c ; r = m.r ; ρ = m.ρ ; σ = m.σ ; a = m.a
    z = state.z
    F, Fz, Fzz = y.F, y.Fz, y.Fzz
    ϕ = z * (1 + z^2)^(-1/2)
    ϕz = (1 + z^2)^(-3/2)
    ϕzz = - 3 * z * (1 + z^2)^(-5/2)
    sτ = c / r * (1 - exp(- r * τ))
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
    return (Ft,), (-ρ * z,), (F = F, i = i, I = i / (sτ * ϕz * exp(r * τ)), x  = sτ * ϕ, I_myopic = i_myopic / (sτ * ϕz * exp(r * τ)), I_hedging = Fz / a / (sτ * ϕz * exp(r * τ)))
end


using EconPDEs
m = ArbitrageHoldingCosts()
state = initialize_state(m)
y0 = initialize_y(m, state)
τs = range(0, stop = 100, length = 100)
y, result, distance = pdesolve(m, state, y0, τs)


## reproduce Fig 2
d = Normal(0, sqrt(m.σ^2 / (2 * m.ρ)))
zmin = quantile(d, 0.025)
zmax = quantile(d, 0.975)
idx = (state[:z] .>= zmin) .& (state[:z] .<= zmax)

using Plots
plot(result[:x][idx, 20], [result[:I_myopic][idx, 2] result[:I][idx, 2]], label = ["myopic" "all"])



plot(result[:x][idx, 20], result[:I][idx, 2])
plot!(result[:x][idx, 60], result[:I][idx, 60])
plot!(result[:x][idx, 100], result[:I][idx, 100])



m = ArbitrageHoldingCosts(; ρ = 0.00001)
state = initialize_state(m)
y0 = initialize_y(m, state)
τs = range(0, stop = 100, length = 100)
y, result, distance = pdesolve(m, state, y0, τs)

