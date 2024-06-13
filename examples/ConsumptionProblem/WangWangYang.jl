using EconPDEs

Base.@kwdef mutable struct WangWangYangModel
    μ::Float64 = 0.015
    σ::Float64 = 0.1
    r::Float64 = 0.035
    ρ::Float64 = 0.04
    γ::Float64 = 3.0
    ψ::Float64 = 1.1
    wmin::Float64 = 0.0
    wmax::Float64 = 1000.0
end

    
function (m::WangWangYangModel)(state::NamedTuple, y::NamedTuple)
    (; μ, σ, r, ρ, γ, ψ, wmin, wmax) = m
    (; w) = state
    (; p, pw_up, pw_down, pww) = y
    pw = pw_up
    iter = 0
    @label start
    pw = max(pw, sqrt(eps()))
    c = (r + ψ * (ρ - r)) * p * pw^(-ψ)
    μw = (r - μ + σ^2) * w + 1 - c
    if (iter == 0) & (μw <= 0)
        iter += 1
        pw = pw_down
        @goto start
    end
   #  One only needs a ghost node if μw <= 0 (since w^2p_ww = 0). In this case, we obtain a formula for pw so that c <= 1
    if w ≈ wmin && μw <= 0.0
        μw = 0.0
        c = 1.0
        pw = (c / ((r + ψ * (ρ - r))))^(-1 / ψ)
    end
    # At the top, I use the solution of the unconstrainted, i.e. pw = 1 (I could also do reflecting boundary but less elegant)
    pt = - ((((r + ψ * (ρ - r)) * pw^(1 - ψ) - ψ * ρ) / (ψ - 1) + μ - γ * σ^2 / 2) * p + ((r - μ + γ * σ^2) * w + 1) * pw + σ^2 * w^2 / 2  * (pww - γ * pw^2 / p))
    return (; pt)
end

m = WangWangYangModel()
stategrid = OrderedDict(:w => range(m.wmin^(1/2), m.wmax^(1/2), length = 100).^2)
yend = OrderedDict(:p => 1 .+ stategrid[:w])
result = pdesolve(m, stategrid, yend, bc = OrderedDict(:pw => (1.0, 1.0)))
@assert result.residual_norm <= 1e-5


# Alternative solution bypassing pdesolve
# just encode the PDE has a vector equation
using InfinitesimalGenerators
function solve!(pts, m, ws, ps)
    (; μ, σ, r, ρ, γ, ψ, wmin, wmax) = m
    pw_ups = FirstDerivative(ws, ps; direction = :upward, bc = (0.0, 1.0))
    pw_downs = FirstDerivative(ws, ps; direction = :downward, bc = (0.0, 1.0))
    pwws = SecondDerivative(ws, ps, bc = (0.0, 1.0))
    for i in eachindex(ws)
        w = ws[i]
        p, pw_up, pw_down, pww = ps[i], pw_ups[i], pw_downs[i], pwws[i]
        pw = pw_up
        iter = 0
        @label start
        pw = max(pw, sqrt(eps()))
        c = (r + ψ * (ρ - r)) * p * pw^(-ψ)
        μw = (r - μ + σ^2) * w + 1 - c
        if (iter == 0) & (μw <= 0)
            iter += 1
            pw = pw_down
            @goto start
        end
        #  One only needs a ghost node if μw <= 0 (since w^2p_ww = 0). In this case, we obtain a formula for pw so that c <= 1
        if w ≈ wmin && μw <= 0.0
            μw = 0.0
            c = 1.0
            pw = (c / ((r + ψ * (ρ - r))))^(-1 / ψ)
        end
        pts[i] = - ((((r + ψ * (ρ - r)) * pw^(1 - ψ) - ψ * ρ) / (ψ - 1) + μ - γ * σ^2 / 2) * p + ((r - μ + γ * σ^2) * w + 1) * pw + σ^2 * w^2 / 2  * (pww - γ * pw^2 / p))
    end
    return pts
end
m = WangWangYangModel()
ws = range(m.wmin^(1/2), m.wmax^(1/2), length = 100).^2
ps = 1 .+ stategrid[:w]
finiteschemesolve((ydot, y) -> solve!(ydot, m, ws, y), ps)