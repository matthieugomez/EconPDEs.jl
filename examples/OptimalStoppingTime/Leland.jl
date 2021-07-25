using EconPDEs

Base.@kwdef struct LelandModel
    # consumption process parameters
    r::Float64 = 0.02
    μ::Float64 = 0.0
    σ::Float64 = 0.1
    C::Float64 = 0.05 # coupon rate
    τ::Float64 = 0.2 # tax rate
end

function (m::LelandModel)(state::NamedTuple, y::NamedTuple)
    (; r, μ, σ, C, τ) = m
    (; δ) = state
    (; E, Eδ, Eδδ) = y
    μδ, σδ = μ * δ, σ * δ
    f = δ - (1 - τ) * C
    Et = f + Eδ * μδ + 0.5 * Eδδ * σδ^2 - r * E
    return (Et,), (μδ, )
end

m = LelandModel()
stategrid = OrderedDict(:δ => [0:0.005:0.2;])
yend = OrderedDict(:E => [max(δ / (m.r - m.μ) - (1 - m.τ) * m.C / m.r, 0.0) for δ in stategrid[:δ]])
y̲ = zeros(length(stategrid[:δ]))
bc = OrderedDict(:Eδ => (0.0, 1 / (m.r - m.μ),))
y, result, distance = pdesolve(m, stategrid, yend; y̲ = y̲, bc = bc, reformulation = :smooth)


"""
    Analytical solution to Leland model. Given current cashflow δ, return the equity value.
"""
function analytical_solution(m::LelandModel, δ)
    (; r, μ, σ, C, τ) = m
    Δ = (1/2 * σ^2 - μ)^2 + 2 * σ^2 * r
    γ = (μ - 1/2 * σ^2 + sqrt(Δ))/σ^2
    η = -(μ - 1/2 * σ^2 - sqrt(Δ))/σ^2
    δb = (1-τ) * C * (r-μ)/r * γ/(1+γ)
    if δ < δb
        return zero(δ)
    else
        return δ / (r - μ) - (1 - τ) * C / r + ((1 - τ) * C / r - δb / (r - μ)) * (δ / δb)^(-γ)
    end
end
