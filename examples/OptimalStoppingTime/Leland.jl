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
    (; E, Eδ_up, Eδ_down, Eδδ) = y
    μδ, σδ = μ * δ, σ * δ
    Eδ = (μδ >= 0) ? Eδ_up : Eδ_down
    f = δ - (1 - τ) * C
    Et = - (f + Eδ * μδ + 0.5 * Eδδ * σδ^2 - r * E)
    return (; Et)
end

m = LelandModel()
stategrid = OrderedDict(:δ => range(0, 0.2, step = 0.005))
yend = OrderedDict(:E => max.(stategrid[:δ] ./ (m.r - m.μ) .- (1 .- m.τ) .* m.C ./ m.r, 0.0))
y̲ = zeros(length(stategrid[:δ]))
bc = OrderedDict(:Eδ => (0.0, 1 / (m.r - m.μ),))
result = pdesolve(m, stategrid, yend; y̲ = y̲, bc = bc, reformulation = :smooth)
@assert result.residual_norm <= 1e-5