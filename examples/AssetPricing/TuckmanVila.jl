using EconPDEs, Distributions

# Tuckman Vila (1992) JF Arbitrage With Holding Costs: A Utility-Based Approach
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
    OrderedDict(:z => range(quantile(d, 0.00001), quantile(d, 0.99999), length = zn))
end

function initialize_y(m::TuckmanVilaModel, stategrid)
    zn = length(stategrid[:z])
    OrderedDict(:F => zeros(zn))
end

function (m::TuckmanVilaModel)(state::NamedTuple, y::NamedTuple, τ::Number)
    (; c, r, ρ, σ, a, T) = m
    (; z) = state
    (; F, Fz_up, Fz_down, Fzz) = y
    Fz = (z >= 0) ? Fz_up : Fz_down
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
    Ft = - ((μ + σ^2 * Fz) * a * i - 0.5 * σ^2 * (a * i)^2 - ρ * z * Fz + 0.5 * σ^2 * (Fzz - Fz^2))
    return (; Ft)
end

m = TuckmanVilaModel()
stategrid = initialize_stategrid(m)
yend = initialize_y(m, stategrid)
τs = range(0, m.T, length = 100)
result = pdesolve(m, stategrid, yend, τs)
residual_norm = maximum(result.residual_norm)

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

