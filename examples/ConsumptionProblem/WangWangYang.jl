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
    # Newton can try negative marginal values, so cap implied consumption instead of flooring derivatives.
    cmax = 100.0 * (1 + max((r - μ + σ^2) * w, 0.0))

    # Upwind on the physical wealth drift of the controlled state process.
    c_up = pw_up > 0 ? min((r + ψ * (ρ - r)) * p * pw_up^(-ψ), cmax) : cmax
    μw_up = (r - μ + σ^2) * w + 1 - c_up
    if μw_up >= 0
        pw = pw_up
        c = c_up
        μw = μw_up
    else
        c_down = pw_down > 0 ? min((r + ψ * (ρ - r)) * p * pw_down^(-ψ), cmax) : cmax
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
            pw = (c / ((r + ψ * (ρ - r)) * p))^(-1 / ψ)
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
        cmax = 100.0 * (1 + max((r - μ + σ^2) * w, 0.0))
        c_up = pw_up > 0 ? min((r + ψ * (ρ - r)) * p * pw_up^(-ψ), cmax) : cmax
        μw_up = (r - μ + σ^2) * w + 1 - c_up
        if μw_up >= 0
            pw = pw_up
            c = c_up
            μw = μw_up
        else
            c_down = pw_down > 0 ? min((r + ψ * (ρ - r)) * p * pw_down^(-ψ), cmax) : cmax
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
                pw = (c / ((r + ψ * (ρ - r)) * p))^(-1 / ψ)
            end
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
