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
    A = r + ψ * (ρ - r)

    # Upwind on the physical wealth drift of the controlled state process.
    pw_up = max(pw_up, sqrt(eps()))
    c_up = A * p * pw_up^(-ψ)
    μw_up = (r - μ + σ^2) * w + 1 - c_up
    if μw_up >= 0
        pw = pw_up
        c = c_up
        μw = μw_up
    else
        pw_down = max(pw_down, sqrt(eps()))
        c_down = A * p * pw_down^(-ψ)
        μw_down = (r - μ + σ^2) * w + 1 - c_down
        if (μw_down <= 0) && (w > wmin)
            pw = pw_down
            c = c_down
            μw = μw_down
        else
            # If the two candidates straddle zero OR drift is negative at minimum asset threshold  
            # we impose drift μw = 0.
            μw = 0.0
            c = 1 + (r - μ + σ^2) * w
            pw = (c / (A * p))^(-1 / ψ)
        end
    end
    # At the top, I use the solution of the unconstrainted, i.e. pw = 1 (I could also do reflecting boundary but less elegant)
    pt = - ((((r + ψ * (ρ - r)) * pw^(1 - ψ) - ψ * ρ) / (ψ - 1) + μ - γ * σ^2 / 2) * p + ((r - μ + γ * σ^2) * w + 1) * pw + σ^2 * w^2 / 2  * (pww - γ * pw^2 / p))
    return (; pt)
end

m = WangWangYangModel()
stategrid = OrderedDict(:w => range(m.wmin^(1/2), m.wmax^(1/2), length = 100).^2)
yend = OrderedDict(:p => 1 .+ stategrid[:w])
@time result = pdesolve(m, stategrid, yend, bc = OrderedDict(:pw => (1.0, 1.0)))
@assert result.residual_norm <= 1e-5


# Alternative solution bypassing pdesolve
# just encode the PDE has a vector equation
using InfinitesimalGenerators

function solve!(pts, m, ws, ps)
    (; μ, σ, r, ρ, γ, ψ, wmin, wmax) = m
    pw_ups = FirstDerivative(ws, ps; direction = :forward, bc = (0.0, 1.0))
    pw_downs = FirstDerivative(ws, ps; direction = :backward, bc = (0.0, 1.0))
    pwws = SecondDerivative(ws, ps, bc = (0.0, 1.0))
    for i in eachindex(ws)
        w = ws[i]
        p, pw_up, pw_down, pww = ps[i], pw_ups[i], pw_downs[i], pwws[i]
        A = r + ψ * (ρ - r)

        pw_up = max(pw_up, sqrt(eps()))
        c_up = A * p * pw_up^(-ψ)
        μw_up = (r - μ + σ^2) * w + 1 - c_up
        if μw_up >= 0
            pw = pw_up
            c = c_up
            μw = μw_up
        else
            pw = pw_down
            pw = max(pw, sqrt(eps()))
            c = A * p * pw^(-ψ)
            μw = (r - μ + σ^2) * w + 1 - c
        end
        # One only needs a ghost node if μw <= 0 (since w^2p_ww = 0). In this case, we obtain a formula for pw so that c <= 1.
        if w ≈ wmin && μw <= 0.0
            μw = 0.0
            c = 1.0
            pw = (c / (A * p))^(-1 / ψ)
        end
        pts[i] = - ((((r + ψ * (ρ - r)) * pw^(1 - ψ) - ψ * ρ) / (ψ - 1) + μ - γ * σ^2 / 2) * p + ((r - μ + γ * σ^2) * w + 1) * pw + σ^2 * w^2 / 2  * (pww - γ * pw^2 / p))
    end
    return pts
end

m = WangWangYangModel()
ws = range(m.wmin^(1/2), m.wmax^(1/2), length = 100).^2
ps = 1 .+ ws
out = finiteschemesolve((ydot, y) -> solve!(ydot, m, ws, y), ps)

@assert sum((out[1] .- result.zero[:p]).^2) / length(out[1]) < 1e-5
