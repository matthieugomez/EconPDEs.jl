using EconPDEs
using Parameters

@with_kw struct LelandModel
    @deftype Float64
    # consumption process parameters
    r = 0.02
    μ = 0.0
    σ = 0.1
    C = 0.05 # coupon rate
    τ = 0.2 # tax rate

    Δ = (1/2 * σ^2 - μ)^2 + 2 * σ^2 * r
    γ = (μ - 1/2 * σ^2 + sqrt(Δ))/σ^2
    η = -(μ - 1/2 * σ^2 - sqrt(Δ))/σ^2
    δb = (1-τ) * C * (r-μ)/r * γ/(1+γ)
end

"""
Analytical solution to Leland model. Given current cashflow δ, return the equity value.
"""
function (leland::LelandModel)(δ)
    @unpack_LelandModel leland
    if δ < δb
        return zero(δ)
    else
        return δ/(r-μ) - (1-τ) * C/r + ((1-τ)*C/r - δb/(r-μ))*(δ/δb)^-γ
    end
end


function initialize_y(leland, grid)
    @unpack_LelandModel leland
    Efunc = δ -> δ/(r-μ) - (1-τ) * C/r
    return OrderedDict(:E => max.(Efunc.(grid[:δ]), 0.0))
end


function (leland::LelandModel)(state::NamedTuple, y::NamedTuple)
    @unpack_LelandModel leland
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
@unpack_LelandModel le

y̲, ȳ = zeros(length(grid[:δ])), Inf * ones(length(grid[:δ]))
y, result, distance = pdesolve(le, grid, y0; y̲ = y̲, ȳ = ȳ, bc = OrderedDict(:Eδ => (0.0, 1/(r - μ),)), reformulation = :smooth)



