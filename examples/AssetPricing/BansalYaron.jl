using EconPDEs, Distributions

Base.@kwdef struct BansalYaronModel
    # consumption process parameters
    μbar::Float64 = 0.018
    vbar::Float64 = 0.00073
    κμ::Float64 = 0.252
    νμ::Float64 = 0.528 
    κv::Float64 = 0.156 
    νv::Float64 = 0.00354

    # utility parameters
    ρ::Float64 = 0.024  
    γ::Float64 = 7.5
    ψ::Float64 = 1.5
end


function (m::BansalYaronModel)(state::NamedTuple, y::NamedTuple)
    (; μbar, vbar, κμ, νμ, κv, νv, ρ, γ, ψ) = m
    (; μ, v) = state
    (; p, pμ, pv, pμμ, pμv, pvv) = y
    
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
    
    # Market Pricing
    pt = p * (1 / p + μc + μp - r - κ_Zc * σc - κ_Zμ * σp_Zμ - κ_Zv * σp_Zv)
    #pt = p * (1 / p - ρ + (1 - 1 / ψ) * (μc - 0.5 * γ * σc^2) + μp + 0.5 * (1 / ψ - γ) / (1 - 1 / ψ) * σp2)

    return (pt,), (μμ, μv)
end

# Bansal Yaron (2004)
m = BansalYaronModel()
μn = 30
νn = 30
μdistribution = Normal(m.μbar,  sqrt(m.νμ^2 * m.vbar / (2 * m.κμ)))
μs = range(quantile(μdistribution, 0.025), quantile(μdistribution, 0.975), length = μn)
νdistribution = Gamma(2 * m.κv * m.vbar / m.νv^2, m.νv^2 / (2 * m.κv))
vs = range(quantile(νdistribution, 0.025), νquantile(distribution, 0.975), length = νn)
stategrid = OrderedDict(:μ => μs, :v => vs)
yend = OrderedDict(:p => ones(μn, νn))
y, result, distance = pdesolve(m, stategrid, yend)

# Bansal Yaron (2009)
## m = BansalYaronModel(μbar = 0.018, vbar = 0.00062, κμ = 0.3, νμ = 0.456, κv = 0.012, νv = 0.00472, ρ = 0.0132, γ = 10, ψ = 1.5)
## stategrid = initialize_stategrid(m)
## yend = OrderedDict(:p => ones(length(stategrid[:μ]), length(stategrid[:v])))
## y, result, distance = pdesolve(m, stategrid, yend)