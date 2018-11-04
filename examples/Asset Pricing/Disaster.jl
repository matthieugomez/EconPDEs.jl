##############################################################################
##
## Type
##
##############################################################################
using Distributions

struct DisasterModel{T<: Distribution}
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
function DisasterModel(;μ = 0.025, σ = 0.02, λbar = 0.0355, κλ = 0.08, νλ = 0.067, ZDistribution = Normal(-0.4, 0.25), ρ = 0.012, γ = 3.0, ψ = 1.1, ϕ = 2.6)
    DisasterModel(μ, σ, λbar, κλ, νλ, ZDistribution, ρ, γ, ψ, ϕ)
end

function initialize_state(m::DisasterModel; n = 30)
    μ = m.μ ; σ = m.σ ; λbar = m.λbar ; κλ = m.κλ ; νλ = m.νλ ; ZDistribution = m.ZDistribution ; ρ = m.ρ ; γ = m.γ ; ψ = m.ψ
    λs = collect(range(0.0, stop = 0.1, length = n))
    OrderedDict(:λ => λs)
end

function initialize_y(m::DisasterModel, state)
    OrderedDict(:p => fill(1.0, length(state[:λ])))
end

function (m::DisasterModel)(state, y)
    μ = m.μ ; σ = m.σ ; λbar = m.λbar ; κλ = m.κλ ; νλ = m.νλ ; ZDistribution = m.ZDistribution ; ρ = m.ρ ; γ = m.γ ; ψ = m.ψ ; ϕ = m.ϕ
    λ = state.λ
    p, pλ, pλλ = y.p, y.pλ, y.pλλ

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
    
    return (pt,), (μλ,), (p = p, r = r, κ_Zc = κ_Zc, κ_Zλ = κ_Zλ, η = η)
end


# solve for levered equity claim
function f(m, state, y, r, κ_Zc, κ_Zλ)
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

#using EconPDEs
#m = DisasterModel()
#state = initialize_state(m)
#y0 = initialize_y(m, state)
#y, result, distance = pdesolve(m, state, y0)
#using Interpolations
#r = interpolate((state[:λ],), result[:r], Gridded(Linear()))
#κ_Zc = interpolate((state[:λ],), result[:κ_Zc], Gridded(Linear()))
#κ_Zλ = interpolate((state[:λ],), result[:κ_Zλ], Gridded(Linear()))
#y2, _, distance =  pdesolve((state, y) -> f(m, state, y, r, κ_Zc, κ_Zλ), state, OrderedDict(:pe => y[:p]))

