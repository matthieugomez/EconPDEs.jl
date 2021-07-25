using EconPDEs

struct LelandModel{T<:Real}
    # consumption process parameters
    r::T
    μ::T
    σ::T
    C::T # coupon rate
    τ::T # tax rate
    # parameters for analytical solution
    Δ::T
    γ::T
    η::T
    δb::T
end

function LelandModel(; 
    r = 0.02, 
    μ = 0.0, 
    σ = 0.1, 
    C = 0.05, 
    τ = 0.2, 
    Δ = (1/2 * σ^2 - μ)^2 + 2 * σ^2 * r,
    γ = (μ - 1/2 * σ^2 + sqrt(Δ))/σ^2,
    η = -(μ - 1/2 * σ^2 - sqrt(Δ))/σ^2,
    δb = (1-τ) * C * (r-μ)/r * γ/(1+γ)
    )
    return LelandModel(r, μ, σ, C, τ, Δ, γ, η, δb)
end
"""
Analytical solution to Leland model. Given current cashflow δ, return the equity value.
"""
function (leland::LelandModel)(δ)
    r, μ, C, τ, γ, δb = leland.r, leland.μ, leland.C, leland.τ, leland.γ, leland.δb
    if δ < δb
        return zero(δ)
    else
        return δ/(r-μ) - (1-τ) * C/r + ((1-τ)*C/r - δb/(r-μ))*(δ/δb)^-γ
    end
end


function initialize_y(leland, grid)
    r, μ, τ, C = leland.r, leland.μ, leland.τ, leland.C
    Efunc = δ -> δ/(r-μ) - (1-τ) * C/r
    return OrderedDict(:E => max.(Efunc.(grid[:δ]), 0.0))
end


function (leland::LelandModel)(state::NamedTuple, y::NamedTuple)
    r, μ, σ, τ, C = leland.r, leland.μ, leland.σ, leland.τ, leland.C
    δ = state.δ
    E, Eδ, Eδδ = y.E, y.Eδ, y.Eδδ
    μδ, σδ = μ * δ, σ * δ
    f = δ - (1-τ) * C
    Et = f + Eδ * μδ + 1/2 * Eδδ * σδ^2 - r * E
    return (Et,), (μδ, )
end

le = LelandModel()
δgrid = [0:0.005:0.2;]
grid = OrderedDict(:δ => δgrid)

y0 = initialize_y(le, grid)

y̲ = zeros(length(grid[:δ]))
y, result, distance = pdesolve(le, grid, y0; y̲ = y̲, bc = OrderedDict(:Eδ => (0.0, 1/(le.r - le.μ),)), reformulation = :smooth)



