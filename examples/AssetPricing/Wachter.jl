using EconPDEs, Distributions

struct WachterModel{T<: Distribution}
    # consumption process parameters
    μ::Float64 
    σ::Float64

    # idiosyncratic
    λbar::Float64
    κλ::Float64
    νλ::Float64
    ZDistribution::T

    # utility parameters
    ρ::Float64  
    γ::Float64 
    ψ::Float64
    ϕ::Float64
end

# distirbution of jump comes from Ian Martin higher order cumulants paper
function WachterModel(;μ = 0.025, σ = 0.02, λbar = 0.0355, κλ = 0.08, νλ = 0.067, ZDistribution = Normal(-0.4, 0.25), ρ = 0.012, γ = 3.0, ψ = 1.1, ϕ = 2.6)
    WachterModel(μ, σ, λbar, κλ, νλ, ZDistribution, ρ, γ, ψ, ϕ)
end

function initialize_stategrid(m::WachterModel; n = 30)
    (; μ, σ, λbar, κλ, νλ, ZDistribution, ρ, γ, ψ, ϕ) = m
    OrderedDict(:λ => range(0.0, stop = 0.1, length = n))
end

function initialize_y(m::WachterModel, stategrid::OrderedDict)
    OrderedDict(:p => ones(length(stategrid[:λ])))
end

function (m::WachterModel)(state::NamedTuple, y::NamedTuple)
    (; μ, σ, λbar, κλ, νλ, ZDistribution, ρ, γ, ψ, ϕ) = m
    (; λ) = state
    (; p, pλ, pλλ) = y

    # Drift and volatility of λ, p
    μλ = κλ * (λbar - λ)
    σλ = νλ * sqrt(λ)
    μp = pλ / p * μλ + 0.5 * pλλ / p * σλ^2 
    σp_Zλ = pλ / p * σλ

    # Market price of risk 
    κ_Zc = γ * σ 
    κ_Zλ = (γ * ψ - 1) / (ψ - 1) * σp_Zλ
    η =  λ * (mgf(ZDistribution, 1) - 1 + mgf(ZDistribution, -γ) - 1 - (mgf(ZDistribution, 1 - γ) - 1))

    # Interest rate r
    r = ρ + μ / ψ - (1 + 1 / ψ) / 2 * γ * σ^2 - λ * (mgf(ZDistribution, -γ) - 1) + (1 / ψ - γ) / (1 - γ) * λ * (mgf(ZDistribution, 1 - γ) - 1) - (γ * ψ - 1) / (2 * (ψ - 1)) * σp_Zλ^2

    # Market Pricing
    pt = p * (1 / p  + μ + μp + λ * (mgf(ZDistribution, 1) - 1) - r - κ_Zc * σ - κ_Zλ * σp_Zλ - η)
    
    return (pt,), (μλ,)
end

m = WachterModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
y, result, distance = pdesolve(m, stategrid, y0)

#========================================================================================

Solve for levered equity claim

========================================================================================#
 
function pde_levered(m, state, y, r, κ_Zc, κ_Zλ)
    μ = m.μ ; σ = m.σ ; λbar = m.λbar ; κλ = m.κλ ; νλ = m.νλ ; ZDistribution = m.ZDistribution ; ρ = m.ρ ; γ = m.γ ; ψ = m.ψ ; ϕ = m.ϕ
    λ = state.λ

    pe, peλ, peλλ = y.pe, y.peλ, y.peλλ
    μλ = κλ * (λbar - λ)
    σλ = νλ * sqrt(λ)
    μpe = peλ / pe * μλ + 0.5 * peλλ / pe * σλ^2 
    σpe_Zλ = peλ / pe * σλ
    η =  λ * (mgf(ZDistribution, ϕ) - 1 + mgf(ZDistribution, -γ) - 1 - (mgf(ZDistribution, ϕ - γ) - 1))
    pet = pe * (1 / pe  + μ + μpe + λ * (mgf(ZDistribution, ϕ) - 1) - r(λ) - κ_Zc(λ) * ϕ * σ - κ_Zλ(λ) * σpe_Zλ - η)
    return (pet,), (μλ,)
end
#using Interpolations
#r = interpolate((state[:λ],), result[:r], Gridded(Linear()))
#κ_Zc = interpolate((state[:λ],), result[:κ_Zc], Gridded(Linear()))
#κ_Zλ = interpolate((state[:λ],), result[:κ_Zλ], Gridded(Linear()))
#y2, _, distance =  pdesolve((state, y) -> pde_levered(m, state, y, r, κ_Zc, κ_Zλ), stategrid, OrderedDict(:pe => y[:p]))