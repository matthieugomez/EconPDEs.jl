using Distributions

mutable struct LongRunRiskModel
    # consumption process parameters
    μbar::Float64 
    vbar::Float64
    κμ::Float64 
    νμ::Float64 
    κv::Float64 
    νv::Float64 

    # utility parameters
    ρ::Float64  
    γ::Float64 
    ψ::Float64
end

function LongRunRiskModel(;μbar = 0.018, vbar = 0.00073, κμ = 0.252, νμ = 0.528, κv = 0.156, νv = 0.00354, ρ = 0.024, γ = 7.5, ψ = 1.5)
    LongRunRiskModel(μbar, vbar, κμ, νμ, κv, νv, ρ, γ, ψ)
end


function initialize_state(m::LongRunRiskModel; μn = 30, vn = 30)
    μbar = m.μbar ; vbar = m.vbar ; κμ = m.κμ ; νμ = m.νμ ; κv = m.κv ; νv = m.νv ; ρ = m.ρ ; γ = m.γ ; ψ = m.ψ

    σ = sqrt(νμ^2 * vbar / (2 * κμ))
    μmin = quantile(Normal(μbar, σ), 0.001)
    μmax = quantile(Normal(μbar, σ), 0.999)
    μs = collect(range(μmin, stop = μmax, length = μn))

    α = 2 * κv * vbar / νv^2
    β = νv^2 / (2 * κv)
    vmin = quantile(Gamma(α, β), 0.001)
    vmax = quantile(Gamma(α, β), 0.999)
    vs = collect(range(vmin, stop = vmax, length = vn))
    OrderedDict(:μ => μs, :v => vs)
end

function initialize_y(m::LongRunRiskModel, state)
    OrderedDict(:p => fill(1.0, length(state[:μ]), length(state[:v])))
end


function (m::LongRunRiskModel)(state, y)
    μbar = m.μbar ; vbar = m.vbar ; κμ = m.κμ ; νμ = m.νμ ; κv = m.κv ; νv = m.νv ; ρ = m.ρ ; γ = m.γ ; ψ = m.ψ
    μ, v = state.μ, state.v
    p, pμ, pv, pμμ, pμv, pvv = y.p, y.pμ, y.pv, y.pμμ, y.pμv, y.pvv
    
    # drift and volatility of c, μ, σ, p
    μc = μ
    σc = sqrt(v)
    μμ = κμ * (μbar - μ)
    σμ = νμ * sqrt(v)
    μv = κv * (vbar - v)
    σv = νv * sqrt(v) 
    σp_Zμ = pμ / p * σμ
    σp_Zv = pv / p * σv
    σp2 = σp_Zμ^2 + σp_Zv^2
    μp = pμ / p * μμ + pv / p * μv + 0.5 * pμμ / p * σμ^2 + 0.5 * pvv / p * σv^2

    # Market price of risk κ
    κ_Zc = γ * σc
    κ_Zμ = - (1 - γ * ψ) / (ψ - 1) * σp_Zμ
    κ_Zv = - (1 - γ * ψ) / (ψ - 1) * σp_Zv
    κ2 = κ_Zc^2 + κ_Zμ^2 + κ_Zv^2

    # Risk free rate r
    r = ρ + μc / ψ - (1 + 1 / ψ) / 2 * γ * σc^2 - (γ * ψ - 1) / (2 * (ψ - 1)) * σp2
    # PDE
    #pt = p * (1 / p + μc + μp - r - κ_Zc * σc - κ_Zμ * σp_Zμ - κ_Zv * σp_Zv)
    pt = p * (1 / p - ρ + (1 - 1 / ψ) * (μc - 0.5 * γ * σc^2) + μp + 0.5 * (1 / ψ - γ) / (1 - 1 / ψ) * σp2)

    return (pt,), (μμ, μv), (p = p, r = r, κ_Zc = κ_Zc, κ_Zμ = κ_Zμ, κ_Zv = κ_Zv, σμ = σμ, σμ_Zv = 0.0, σv_Zμ = 0.0, σv = σv, μ = μ, v = v, σμ2 = σμ^2, σv2 = σv^2, σμσv = 0.0, μμ = μμ, μv = μv)
end

## Long Run Risk Models
### Bansal Yaron (2004)
# m = LongRunRiskModel()
# state = initialize_state(m)
# y0 = initialize_y(m, state)
# y, result, distance = pdesolve(m, state, y0)
