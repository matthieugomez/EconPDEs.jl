using EconPDEs

mutable struct WangWangYangModel
    μ::Float64 
    σ::Float64
    r::Float64
    ρ::Float64  
    γ::Float64 
    ψ::Float64
    wmax::Float64
end

function WangWangYangModel(;μ = 0.015, σ = 0.1, r = 0.035, ρ = 0.04, γ = 3, ψ = 1.1, wmax = 1000.0)
    WangWangYangModel(μ, σ, r, ρ, γ, ψ, wmax)
end

function initialize_state(m::WangWangYangModel; n = 500)
    OrderedDict(:w => collect(range(0.0, stop = m.wmax, length = n)))
end

function initialize_y(m::WangWangYangModel, state)
    OrderedDict(:p => 1 .+ state[:w])
end
    
function (m::WangWangYangModel)(state, y)
    μ = m.μ ;  σ = m.σ ;  r = m.r ;  ρ = m.ρ ;  γ = m.γ ;  ψ = m.ψ  ; wmax = m.wmax
    w = state.w
    p, pw, pww = y.p, y.pw, y.pww
    c = (r + ψ * (ρ - r)) * p * pw^(-ψ)
    μw = (r - μ + σ^2) * w + 1 - c
    # Turns out that this branch is unecessary. One way to understand this may be that, since second derivative drops out at the bottom, there is already a boundary condition, which is given by relationship between first and second derivative.
    if w ≈ 0.0 && μw <= 0.0
       pw = ((r + ψ * (ρ - r)) * p)^(1 / ψ)
       c = (r + ψ * (ρ - r)) * p * pw^(-ψ)
       μw = (r - μ + σ^2) * w + 1 - c
    end
    # this branch is unnecessary when individuals dissave at the top (default)
    if w ≈ wmax && μw >= 0.0
        pw = 1.0
        pww = 0.0
        c = (r + ψ * (ρ - r)) * p * pw^(-ψ)
        μw = (r - μ + σ^2) * w + 1 - c
    end
    # Since first derivative is upwinded, I can directly impose value of second derivative
    if w ≈ wmax
        pww = 0.0
    end
    pt = (((r + ψ * (ρ - r)) * pw^(1 - ψ) - ψ * ρ) / (ψ - 1) + μ - γ * σ^2 / 2) * p + ((r - μ + γ * σ^2) * w + 1) * pw + σ^2 * w^2 / 2  * (pww - γ * pw^2 / p)
    μw = 1 + (r - μ + σ^2) * w - c
    return (pt,), (μw,), (w = w, p = p, pw = pw, pww = pww, μw = μw, c = c)
end

m = WangWangYangModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
y0, result0, distance0 = pdesolve(m, state, y0)
# Important: check marginal value of wealth converges to 1.0
result0[:pw]