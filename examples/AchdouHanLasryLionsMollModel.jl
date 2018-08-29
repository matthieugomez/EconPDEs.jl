using Distributions

mutable struct AchdouHanLasryLionsMollModel
    # income process parameters
    κy::Float64 
    ybar::Float64
    σy::Float64

    abar::Float64

    r::Float64

    # utility parameters
    ρ::Float64  
    γ::Float64 
end

function AchdouHanLasryLionsMollModel(;κy = 0.018, ybar = 1.0, σy = 0.05, abar = 0.00, r = 0.01, ρ = 0.024, γ = 7.5)
    AchdouHanLasryLionsMollModel(κy, ybar, σy, abar, r, ρ, γ)
end


function initialize_state(m::AchdouHanLasryLionsMollModel; yn = 10, an = 50)
    κy = m.κy ; ybar = m.ybar ; σy = m.σy ; abar = m.abar ; ρ = m.ρ ; γ = m.γ

    distribution = Gamma(2 * κy * ybar / σy^2, σy^2 / (2 * κy))
    ymin = quantile(distribution, 0.001)
    ymax = quantile(distribution, 0.999)
    ys = collect(range(ymin, stop = ymax, length = yn))

    as = collect(range(abar, stop = 100.0, length = an))

    OrderedDict(:y => ys, :a => as)
end

function initialize_y(m::AchdouHanLasryLionsMollModel, state)
    OrderedDict(:v => fill(1.0, length(state[:y]), length(state[:a])))
end

function (m::AchdouHanLasryLionsMollModel)(state, value)
    κy = m.κy ; σy = m.σy ; ybar = m.ybar ; abar = m.abar ; r = m.r ; ρ = m.ρ ; γ = m.γ
    y, a = state.y, state.a
    v, vy, va, vyy, vya, vaa = value.v, value.vy, value.va, value.vyy, value.vya, value.vaa
    va = max(1e-10, va)
    if a == abar
        c = va^(-1 / γ)
        if c >= y + r * abar
            va = (y + r * abar)^(-γ)
        end
    end
    c = va^(-1 / γ)
    μy = κy * (ybar - y)
    μa = y + r * a - c
    vt = c^(1 - γ) / (1 - γ) + va * μa + vy * μy + 0.5 * vyy * σy^2 - ρ * v
    return (vt,), (μy, μa), (c = c, va = va, vy = vy, y = y, a = a, μa = μa)
end



m = AchdouHanLasryLionsMollModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
result, a, distance = pdesolve(m, state, y0)
using Plots
surface(state[:a], state[:y], a[:μa])
