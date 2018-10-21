using EconPDEs

mutable struct WangWangYangModel
    μ::Float64 
    σ::Float64
    r::Float64
    ρ::Float64  
    γ::Float64 
    ψ::Float64
end

function WangWangYangModel(;μ = 0.015, σ = 0.1, r = 0.035, ρ = 0.04, γ = 3, ψ = 1.1)
    WangWangYangModel(μ, σ, r, ρ, γ, ψ)
end

function initialize_state(m::WangWangYangModel; n = 100)
    OrderedDict(:w => collect(range(0.0, stop = 1_000.0, length = n)))
end

function initialize_y(m::WangWangYangModel, state)
    OrderedDict(:p => 1 .+ state[:w])
end
    
function (m::WangWangYangModel)(state, y)
    μ = m.μ ;  σ = m.σ ;  r = m.r ;  ρ = m.ρ ;  γ = m.γ ;  ψ = m.ψ 
    w = state.w
    p, pw, pww = y.p, y.pw, y.pww
    m = r + ψ * (ρ - r)
    c = m * p * pw^(-ψ)
    pt = ((m * pw^(1 - ψ) - ψ * ρ) / (ψ - 1) + μ - γ * σ^2 / 2) * p + ((r - μ + γ * σ^2) * w + 1) * pw + σ^2 * w^2 / 2  * (pww - γ * pw^2 / p)
    μw = (r - μ + σ^2) * w + 1 - c
    return (pt,), (μw,), (w = w, p = p, pw = pw, pww = pww, μw = μw, c = c)
end

m = WangWangYangModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
y, result, distance = pdesolve(m, state, y0, bc = OrderedDict(:pw => (3.0, 1.0)))


# becuase second derivative drops out at the bottom, no need for boundary condition at the bottom, only at the top
# that's weird cause I don't even have to use the fact that firctions at the bottom?
# also I don't have that derivative equals what I wnat
result[:pw][1] - ((m.r + m.ψ * (m.ρ - m.r)) * result[:p][1])^(1 / m.ψ)
