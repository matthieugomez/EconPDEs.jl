using EconPDEs

struct LelandModel
    # consumption process parameters
    r::Float64
    μ::Float64
    σ::Float64
    C::Float64 # coupon rate
    τ::Float64 # tax rate
end

function LelandModel(; r = 0.02, μ = 0.0, σ = 0.1, C = 0.05, τ = 0.2)
    return LelandModel(r, μ, σ, C, τ)
end

function initialize_y(m::LelandModel, grid)
    r, μ, σ, C, τ = m.r, m.μ, m.σ, m.C, m.τ
    Efunc = δ -> δ/(r-μ) - (1-τ) * C/r
    return OrderedDict(:E => max.(Efunc.(grid[:δ]), 0.0))
end


function (m::LelandModel)(state::NamedTuple, y::NamedTuple)
    r, μ, σ, C, τ = m.r, m.μ, m.σ, m.C, m.τ
    δ = state.δ
    E, Eδ, Eδδ = y.E, y.Eδ, y.Eδδ
    μδ, σδ = μ * δ, σ * δ
    f = δ - (1-τ) * C
    Et = f + Eδ * μδ + 1/2 * Eδδ * σδ^2 - r * E
    return (Et,), (μδ, )
end

le = LelandModel()
grid = OrderedDict(:δ => [0:0.005:0.2;])
y0 = initialize_y(le, grid)
y̲ = zeros(length(grid[:δ]))
y, result, distance = pdesolve(le, grid, y0; y̲ = y̲, bc = OrderedDict(:Eδ => (0.0, 1/(le.r - le.μ),)), reformulation = :smooth)


"""
    Analytical solution to Leland model. Given current cashflow δ, return the equity value.
"""
function analytical_solution(m::LelandModel, δ)
    r, μ, σ, C, τ = m.r, m.μ, m.σ, m.C, m.τ
    Δ = (1/2 * σ^2 - μ)^2 + 2 * σ^2 * r
    γ = (μ - 1/2 * σ^2 + sqrt(Δ))/σ^2
    η = -(μ - 1/2 * σ^2 - sqrt(Δ))/σ^2
    δb = (1-τ) * C * (r-μ)/r * γ/(1+γ)
    if δ < δb
        return zero(δ)
    else
        return δ/(r-μ) - (1-τ) * C/r + ((1-τ)*C/r - δb/(r-μ))*(δ/δb)^-γ
    end
end
