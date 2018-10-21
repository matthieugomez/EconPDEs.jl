using EconPDEs, Distributions

mutable struct AchdouHanLasryLionsMollModel
    # income process parameters
    κy::Float64 
    ybar::Float64
    σy::Float64

    r::Float64

    # utility parameters
    ρ::Float64  
    γ::Float64

    amin::Float64
    amax::Float64 
end

function AchdouHanLasryLionsMollModel(;κy = 0.1, ybar = 1.0, σy = 0.07, r = 0.03, ρ = 0.05, γ = 2.0, amin = 0.0, amax = 50.0)
    AchdouHanLasryLionsMollModel(κy, ybar, σy, r, ρ, γ, amin, amax)
end


function initialize_state(m::AchdouHanLasryLionsMollModel; yn = 20, an = 50)
    κy = m.κy ; ybar = m.ybar ; σy = m.σy  ; ρ = m.ρ ; γ = m.γ ; amin = m.amin ; amax = m.amax

    distribution = Gamma(2 * κy * ybar / σy^2, σy^2 / (2 * κy))
    ymin = quantile(distribution, 0.001)
    ymax = quantile(distribution, 0.999)
    ys = collect(range(ymin, stop = ymax, length = yn))
    as = collect(range(amin, stop = amax, length = an))
    OrderedDict(:y => ys, :a => as)
end

function initialize_y(m::AchdouHanLasryLionsMollModel, state)
    OrderedDict(:v => [(y + m.r * a)^(1-m.γ)/(1-m.γ)/m.ρ for y in state[:y], a in state[:a]])
end

function (m::AchdouHanLasryLionsMollModel)(state, value)
    κy = m.κy ; σy = m.σy ; ybar = m.ybar ; r = m.r ; ρ = m.ρ ; γ = m.γ ; amin = m.amin ; amax = m.amax
    y, a = state.y, state.a
    v, vy, va, vyy, vya, vaa = value.v, value.vy, value.va, value.vyy, value.vya, value.vaa
    μy = κy * (ybar - y)
    va = max(va, eps())
    
    # There is no second derivative with respect to a. Therefore, a simple way to impose boundary conditions is to specify them directly
    c = va^(-1 / γ)
    μa = y + r * a - c
    if (a ≈ amin) && (μa <= 0.0)
        va = (y + r * amin)^(-γ)
        c = y + r * amin
        μa = 0.0
    end
    # does not matter if individuals dissave at the top, which is true for the set of parameters
    if (a ≈ amax) && (μa >= 0.0)
        va = (((ρ - r) / γ + r) * a)^(-γ)
        c = ((ρ - r) / γ + r) * a
        μa = y + (r - ρ) / γ * a
    end

    vt = c^(1 - γ) / (1 - γ) + va * μa + vy * μy + 0.5 * vyy * σy^2 - ρ * v
    return (vt,), (μy, μa), (c = c, va = va, vy = vy, y = y, a = a, μa = μa)
end



m = AchdouHanLasryLionsMollModel()
state = initialize_state(m)
y0 = initialize_y(m, state)
y, result, distance = pdesolve(m, state, y0)
using Plots
surface(state[:a], state[:y], result[:μa])
