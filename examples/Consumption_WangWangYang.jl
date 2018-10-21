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
    # financial friction: check consumption < 1 when w = 0. Turns out that does not matter.
    # I think that one way to undertstand this is that, since second derivative drops out at the bottom, there is already a boundary condition given by relationship between first and second derivative.
    if w == 0.0
       m = r + ψ * (ρ - r)
       c = m * p * pw^(-ψ)
       if c >= 1.0
           pw = (m * p)^(1 / ψ)
       end
    end
    pt = ((m * pw^(1 - ψ) - ψ * ρ) / (ψ - 1) + μ - γ * σ^2 / 2) * p + ((r - μ + γ * σ^2) * w + 1) * pw + σ^2 * w^2 / 2  * (pww - γ * pw^2 / p)
    μw = (r - μ + σ^2) * w + 1 - c
    return (pt,), (μw,), (w = w, p = p, pw = pw, pww = pww, μw = μw, c = c)
end

m = WangWangYangModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
y, result, distance = pdesolve(m, state, y0, bc = OrderedDict(:pw => (3.0, 1.0)))