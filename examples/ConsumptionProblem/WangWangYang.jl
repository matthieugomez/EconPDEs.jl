using EconPDEs

struct WangWangYangModel
    μ::Float64 
    σ::Float64
    r::Float64
    ρ::Float64  
    γ::Float64 
    ψ::Float64
    wmax::Float64
end

function WangWangYangModel(;μ = 0.015, σ = 0.1, r = 0.035, ρ = 0.04, γ = 3, ψ = 1.1, wmax = 5000.0)
    WangWangYangModel(μ, σ, r, ρ, γ, ψ, wmax)
end

function initialize_stategrid(m::WangWangYangModel; n = 500)
    OrderedDict(:w => range(0.0, stop = sqrt(m.wmax), length = n).^2)
end

function initialize_y(m::WangWangYangModel, stategrid)
    OrderedDict(:p => 1 .+ stategrid[:w])
end
    
function (m::WangWangYangModel)(state::NamedTuple, y::NamedTuple)
    μ = m.μ ;  σ = m.σ ;  r = m.r ;  ρ = m.ρ ;  γ = m.γ ;  ψ = m.ψ ; wmax = m.wmax
    w = state.w
    p, pw, pww = y.p, y.pw, y.pww
    c = (r + ψ * (ρ - r)) * p * pw^(-ψ)
    μw = (r - μ + σ^2) * w + 1 - c
    # One only needs a ghost node if μw <= 0 (since w^2p_ww = 0). In this case, we obtain a formula for pw so that c <= 1
    if w ≈ 0.0 && μw <= 0.0
       pw = ((r + ψ * (ρ - r)) * p)^(1 / ψ)
       c = (r + ψ * (ρ - r)) * p * pw^(-ψ)
       μw = (r - μ + σ^2) * w + 1 - c
    end
    # At the top, I use the solution of the unconstrainted, i.e. pw = 1 (I could also do reflecting boundary but less elegant)
    pt = (((r + ψ * (ρ - r)) * pw^(1 - ψ) - ψ * ρ) / (ψ - 1) + μ - γ * σ^2 / 2) * p + ((r - μ + γ * σ^2) * w + 1) * pw + σ^2 * w^2 / 2  * (pww - γ * pw^2 / p)
    μw = 1 + (r - μ + σ^2) * w - c
    return (pt,), (μw,), (w = w, p = p, pw = pw, pww = pww, μw = μw, c = c)
end

m = WangWangYangModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
y, result, distance = pdesolve(m, stategrid, y0; bc = OrderedDict(:pw => (1.0, 1.0)))