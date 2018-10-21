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

function initialize_state(m::WangWangYangModel; n = 300)
    OrderedDict(:w => collect(range(0.0, stop = 500.0, length = n)))
end

function initialize_y(m::WangWangYangModel, state)
    OrderedDict(:p => 1 .+ state[:w])
end
	
function (m::WangWangYangModel)(state, y)
    μ = m.μ ;  σ = m.σ ;  r = m.r ;  ρ = m.ρ ;  γ = m.γ ;  ψ = m.ψ 
    w = state.w
    p, pw, pww = y.p, y.pw, y.pww
    # financial friction: check consumption < 1 when w = 0
    pt = 0.0
    if w == 0.0
        m = r + ψ * (ρ - r)
        c = m * p * pw^(-ψ)
        if c >= 1.0
            pw =  (m * p)^(1 / ψ)
        end
    end
    m = r + ψ * (ρ - r)
    c = m * p * pw^(-ψ)
    pt = ((m * pw^(1 - ψ) - ψ * ρ) / (ψ - 1) + μ - γ * σ^2 / 2) * p + ((r - μ + γ * σ^2) * w + 1) * pw + σ^2 * w^2 / 2  * (pww - γ * pw^2 / p)
    μw = (r - μ + σ^2) * w + 1 - c
    return (pt,), (μw,), (w = w, p = p, pw = pw, pww = pww, μw = μw, c = c)
end

m = WangWangYangModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
y2, result2, distance = pdesolve(m, state, y0, bc = OrderedDict(:pw => (3.0, 1.0)))

# boundary condition  for derivative value at lower boundary is not used. Note that at the frontier, we have w = 0 so second derivative drops out. 
# what I am imposing right now is a condition about value of p at boundary tijme = 0 (since pww disappears and pw is a function of p, the PDE at w = 0 is just a condition on p) and vlaue of pw at upward boundary.
# Note that I do not really satisfy  y2[:p][end] - state[:w][end] - 1 / (m.r - m.μ) but I  think it is true in the limit